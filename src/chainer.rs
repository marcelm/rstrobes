use std::cmp::Reverse;
use std::time::Instant;
use crate::details::NamDetails;
use crate::hit::{Hit, find_hits};
use crate::index::StrobemerIndex;
use crate::mapper::{MappingParameters, QueryRandstrobe};
use crate::mcsstrategy::McsStrategy;
use crate::nam::Nam;

const N_PRECOMPUTED: usize = 1024;

#[derive(Debug, Ord, PartialOrd, Eq, PartialEq)]
struct Anchor {
    ref_id: usize,
    ref_start: usize,
    query_start: usize,
}

#[derive(Debug)]
struct ChainingParameters {
    max_lookback: usize,
    diag_diff_penalty: f32,
    gap_length_penalty: f32,
    valid_score_threshold: f32,
    max_ref_gap: usize,
    matches_weight: f32,
}

#[derive(Debug)]
struct Chainer {
    k: usize,
    parameters: ChainingParameters,
    precomputed_scores: [f32; N_PRECOMPUTED],
}

impl Chainer {
    pub fn new(k: usize, parameters: ChainingParameters) -> Self {
        let mut precomputed_scores: [f32; N_PRECOMPUTED];
        for i in 1..N_PRECOMPUTED {
            precomputed_scores[i] = compute_score(i, i, k, &parameters);
        }
        Chainer { k, parameters, precomputed_scores }
    }

    fn compute_score(&self, dq: i64, dr: i64) -> f32 {
        // dq == dr is usually the most common case
        if dq == dr && (dq as usize) < N_PRECOMPUTED {
            self.precomputed_scores[dq]
        } else {
            compute_score(dq, dr, self.k, &self.parameters)
        }
    }

    fn collinear_chaining(
        &self,
        anchors: &[Anchor],
        // dp: &[f32],
        // predecessors: &[usize],
    ) -> (f32, Vec<f32>, Vec<usize>) {
        let n = anchors.len();
        if n == 0 {
            return (0.0, vec![], vec![]);
        }

        let dp = vec![self.k as f32; n];
        let predecessors = vec![usize::MAX; n];
        let mut best_score = 0.0;

        for i in 0..n {
            let lookup_end = i.saturating_sub(self.parameters.max_lookback);
            for j in (lookup_end..i).rev() {
                let ai = &anchors[i];
                let aj = &anchors[j];

                if ai.ref_id != aj.ref_id {
                    break;
                }

                let dq = ai.query_start - aj.query_start;
                let dr = ai.ref_start - aj.ref_start;

                if dr >= self.parameters.max_ref_gap {
                    break;
                }
                if dq <= 0 || dr <= 0 {
                    continue;
                }

                let score = self.compute_score(dq, dr);

                let new_score = dp[j] + score;
                if new_score > dp[i] {
                    dp[i] = new_score;
                    predecessors[i] = j;
                    if dq == dr {
                        break;
                    }
                }
            }
            if dp[i] > best_score {
                best_score = dp[i];
            }
        }

        (best_score, dp, predecessors)
    }

    fn get_chains(
        &self,
        query_randstrobes: &[Vec<QueryRandstrobe>; 2],
        index: &StrobemerIndex,
        map_param: &MappingParameters,
    ) -> (NamDetails, Vec<Nam>) {
        let hits_timer = Instant::now();

        let mut total_hits = 0;
        let mut partial_hits = 0;
        let mut hits = [vec![], vec![]];
        for is_revcomp in 0..2 {
            let total_hits1;
            let partial_hits1;
            let sorting_needed1;

            (total_hits1, partial_hits1, sorting_needed1, hits[is_revcomp]) =
                find_hits(&query_randstrobes[is_revcomp], index, index.filter_cutoff, map_param.mcs_strategy);
            total_hits += total_hits1;
            partial_hits += partial_hits1;
        }
        let nonrepetitive_hits = hits[0].len() + hits[1].len();
        let nonrepetitive_fraction = if total_hits > 0 { (nonrepetitive_hits as f32) / (total_hits as f32) } else { 1.0 };
        let mut time_find_hits = hits_timer.elapsed().as_secs_f64();

        let mut n_rescue_hits = 0;
        let mut n_rescue_nams = 0;
        let mut n_hits = nonrepetitive_hits;
        let mut n_partial_hits = partial_hits;
        let mut rescue_done = false;
        let mut time_rescue;
        let mut time_chaining;
        let mut chains = vec![];
        for is_revcomp in 0..2 {
            let mut best_score = 0.0;
            let mut anchors = vec![];
            let mut dp = vec![];
            let mut predecessors = vec![];

            // Rescue if requested and needed
            if map_param.rescue_level > 1 && (nonrepetitive_hits == 0 || nonrepetitive_fraction < 0.7) {
                let rescue_timer = Instant::now();
                rescue_done = true;
                let (n_rescue_hits1, n_partial_hits1) = find_anchors_rescue(
                    &query_randstrobes[is_revcomp], index, index.rescue_cutoff, map_param.mcs_strategy, &mut anchors
                );
                n_rescue_hits += n_rescue_hits1;
                n_partial_hits += n_partial_hits1;
                n_rescue_nams += chains.len();
                rescue_done = true;
                time_rescue += rescue_timer.elapsed().as_secs_f64();
            } else {
                let hits_timer = Instant::now();
                add_hits_to_anchors(&hits[is_revcomp], index, &mut anchors);
                time_find_hits += hits_timer.elapsed().as_secs_f64();
            }
            let chaining_timer = Instant::now();
            pdqsort_branchless(anchors.begin(), anchors.end());
            anchors.erase(
                std::unique(anchors.begin(), anchors.end()), anchors.end()
            );
            let (score, dp, predecessors) = self.collinear_chaining(&anchors);
            best_score = score;

            extract_chains_from_dp(
                &anchors, &dp, &predecessors, best_score,
                index.k(), is_revcomp == 1, &mut chains, &self.parameters
            );
            time_chaining += chaining_timer.elapsed().as_secs_f64();
        }
        let details = NamDetails {
            n_reads: 1,
            n_randstrobes,
            n_nams: chains.len(),
            n_rescue_nams,
            nam_rescue: rescue_done as usize,
            n_hits: nonrepetitive_hits,
            n_partial_hits: partial_hits,
            n_rescue_hits,
            time_randstrobes,
            time_find_hits,
            time_chaining,
            time_rescue,
            time_sort_nams,
        };

        (details, chains)
    }
}

/**
* @brief Compute the chaining score between two anchors based on their distance.
*
* This function calculates the score contribution when chaining two anchors,
* considering their relative positions on the query and reference sequences.
* It penalizes large gaps and diagonal differences to encourage collinear chains.
*
* @param dq Difference in query start positions between anchor i and anchor j
*           (i.e., ai.query_start - aj.query_start). Must be > 0.
* @param dr Difference in reference start positions between anchor i and anchor j
*           (i.e., ai.ref_start - aj.ref_start). Must be > 0.
* @param k  Length of the k-mer used to form the anchor.
* @param chaining_params Parameters controlling penalties:
*        - diag_diff_penalty: Multiplier for the absolute difference |dr - dq|, penalizing non-diagonal moves.
*        - gap_length_penalty: Multiplier for min(dq, dr), penalizing longer gaps.
*
* @return A float representing the score for chaining the two anchors.
*/
fn compute_score(dq: i32, dr: i32, k: i32, parameters: &ChainingParameters) -> f32 {
    let dd = dr.abs_diff(dq);
    let dg = dq.min(dr);
    let mut score = k.min(dg) as f32;

    let lin_penalty = parameters.diag_diff_penalty * dd as f32 + parameters.gap_length_penalty * dg as f32;
    let log_penalty = if dd >= 1 { (dd as f32 + 1f32).log2() } else { 0.0 };
    score -= lin_penalty + 0.5 * log_penalty;

    score
}

fn add_to_anchors_full(
    anchors: &mut Vec<Anchor>,
    query_start: usize,
    query_end: usize,
    index: &StrobemerIndex,
    mut position: usize,
) {
    let mut min_length_diff = usize::MAX;
    for randstrobe in &index.randstrobes[position..] {
        if randstrobe.hash() != index.randstrobes[position].hash() {
            break;
        }
        let ref_start = randstrobe.position();
        let ref_end = ref_start + randstrobe.strobe2_offset() + index.k();
        let length_diff = (query_end - query_start).abs_diff(ref_end - ref_start);
        if length_diff <= min_length_diff {
            let ref_id = randstrobe.reference_index();
            anchors.push(Anchor {ref_id, ref_start, query_start});
            anchors.push(Anchor {ref_id, ref_start: ref_end - index.k(), query_start: query_end - index.k()});
            min_length_diff = length_diff;
        }
        position += 1;
    }
}

fn add_to_anchors_partial(
    anchors: &mut Vec<Anchor>,
    query_start: usize,
    index: &StrobemerIndex,
    position: usize,
) {
    let hash = index.get_hash_partial(position);
    for pos in position..index.randstrobes.len() {
        if index.get_hash_partial(pos) != hash {
            break;
        }
        let randstrobe = &index.randstrobes[pos];
        let ref_id = randstrobe.reference_index();
        let (ref_start, ref_end) = index.strobe_extent_partial(pos);

        anchors.push(Anchor { ref_id, ref_start, query_start });
    }
}

fn add_hits_to_anchors(
    hits: &Vec<Hit>,
    index: &StrobemerIndex,
    anchors: &mut Vec<Anchor>,
) {
    for hit in hits {
        if hit.is_partial {
            add_to_anchors_partial(anchors, hit.query_start, index, hit.position);
        } else {
            add_to_anchors_full(anchors, hit.query_start, hit.query_end, index, hit.position);
        }
    }
}

fn find_anchors_rescue(
    query_randstrobes: &Vec<QueryRandstrobe>,
    index: &StrobemerIndex,
    rescue_cutoff: usize,
    mcs_strategy: McsStrategy,
    anchors: &mut Vec<Anchor>,
) -> (usize, usize) {
    struct RescueHit {
        count: usize,
        query_start: usize,
        query_end: usize,
        position: usize,
        is_partial: bool,
    }

    let mut n_hits = 0;
    let mut partial_hits = 0;
    let mut rescue_hits = vec![];
    for qr in query_randstrobes {
        if let Some(position) = index.get_full(qr.hash) {
            let mut count = index.get_count_full(position);
            if let Some(position_revcomp) = index.get_full(qr.hash_revcomp) {
                count += index.get_count_full(position_revcomp);
            }
            rescue_hits.push(RescueHit {position, count, query_start: qr.start, query_end: qr.end, is_partial: false});
        } else if mcs_strategy == McsStrategy::Always {
            if let Some(partial_pos) = index.get_partial(qr.hash) {
                let mut partial_count = index.get_count_partial(partial_pos);
                if let Some(position_revcomp) = index.get_partial(qr.hash_revcomp) {
                    partial_count += index.get_count_partial(position_revcomp);
                }
                rescue_hits.push(RescueHit {position: partial_pos, count: partial_count, query_start: qr.start, query_end: qr.start + index.k(), is_partial: true});
                partial_hits += 1;
            }
        }
    }

    let cmp = |a: &RescueHit, b: &RescueHit| (a.count, a.query_start, a.query_end).cmp(&(b.count, b.query_start, b.query_end));
    rescue_hits.sort_by(cmp);

    for (i, rh) in rescue_hits.iter().enumerate() {
        if (rh.count > rescue_cutoff && i >= 5) || rh.count > 1000 {
            break;
        }
        if rh.is_partial {
            partial_hits += 1;
            add_to_anchors_partial(anchors, rh.query_start, index, rh.position);
        } else {
            add_to_anchors_full(anchors, rh.query_start, rh.query_end, index, rh.position);
        }
        n_hits += 1;
    }

    (n_hits, partial_hits)
}

fn extract_chains_from_dp(
    anchors: &[Anchor],
    dp: &[f32],
    predecessors: &[usize],
    best_score: f32,
    k: usize,
    is_revcomp: bool,
    chains: &mut Vec<Nam>,
    chaining_parameters: &ChainingParameters,
) {
    let n = anchors.len();
    let valid_score = best_score * chaining_parameters.valid_score_threshold;
    let mut used = vec![false; n];

    let mut candidates = vec![];//std::vector<std::pair<int, float>> candidates;
    for i in 0..n {
        if dp[i] >= valid_score {
            candidates.push((i, dp[i]));
        }
    }

    candidates.sort_by(|a, b| b.1.total_cmp(&a.1));

    for (i, score) in candidates {
        if used[i] {
            continue;
        }

        let mut j = i;
        let mut c = 1;
        let mut overlaps = false;

        while predecessors[j] >= 0 {
            j = predecessors[j];
            if used[j] {
                overlaps = true;
                break;
            }
            used[j] = true;
            c += 1;
        }

        if overlaps {
            continue;
        }

        let first = &anchors[j];
        let last = &anchors[j];
        chains.push(
            Nam {
                nam_id: chains.len(),
                query_start: first.query_start,
                query_end: last.query_start + k,
                query_prev_match_startpos: usize::MAX,
                ref_start: first.ref_start,
                ref_end: last.ref_start + k,
                ref_prev_match_startpos: usize::MAX,
                n_matches: c,
                ref_id: last.ref_id,
                score: score + c as f32 * chaining_parameters.matches_weight,
                is_revcomp,
            }
        );
    }
}
