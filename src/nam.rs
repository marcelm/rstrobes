use std::cmp::{max, min};
use std::collections::HashMap;
use crate::fasta::RefSequence;
use crate::index::StrobemerIndex;
use crate::mapper::QueryRandstrobe;
use crate::read::Read;

/// Non-overlapping approximate match
#[derive(Clone,Copy,Debug)]
pub struct Nam {
    nam_id: usize,
    query_start: usize,
    query_end: usize,
    query_prev_hit_startpos: usize,
    ref_start: usize,
    ref_end: usize,
    ref_prev_hit_startpos: usize,
    pub n_hits: usize,
    ref_id: usize,
    pub score: u32,
    is_revcomp: bool,
}

impl Nam {
    fn ref_span(&self) -> usize {
        self.ref_end - self.ref_start
    }

    fn query_span(&self) -> usize {
        self.query_end - self.query_start
    }
}

struct Hit {
    query_start: usize,
    query_end: usize,
    ref_start: usize,
    ref_end: usize,
    is_revcomp: bool,
}

/// Find a query’s NAMs, ignoring randstrobes that occur too often in the
/// reference (have a count above filter_cutoff).
///
/// Return the fraction of nonrepetitive hits (those not above the filter_cutoff threshold)
///
pub fn find_nams(query_randstrobes: &Vec<QueryRandstrobe>, index: &StrobemerIndex) -> (f32, Vec<Nam>) {
    let mut hits_per_ref = HashMap::new();
    hits_per_ref.reserve(100);
    let mut nr_good_hits = 0;
    let mut total_hits = 0;
    for randstrobe in query_randstrobes {
        if let Some(position) = index.get(randstrobe.hash) {
            total_hits += 1;
            if index.is_filtered(position) {
                continue;
            }
            nr_good_hits += 1;
            add_to_hits_per_ref(&mut hits_per_ref, randstrobe.start, randstrobe.end, randstrobe.is_reverse, index, position, 100_000);
        }
    }
    let nonrepetitive_fraction = if total_hits > 0 { (nr_good_hits as f32) / (total_hits as f32) } else { 1.0 };
    let nams = merge_hits_into_nams(&mut hits_per_ref, index.parameters.syncmer.k, false);

    (nonrepetitive_fraction, nams)
}

/// Find a query’s NAMs, using also some of the randstrobes that occur more often
/// than filter_cutoff.
fn find_nams_rescue(
    query_randstrobes: &Vec<QueryRandstrobe>,
    index: &StrobemerIndex,
    filter_cutoff: u32
) -> Vec<Nam> {
/*
    struct RescueHit {
        unsigned int count;
        size_t position;
        unsigned int query_s;
        unsigned int query_e;
        bool is_revcomp;

        bool operator< (const RescueHit& rhs) const {
            return std::tie(count, query_s, query_e, is_revcomp)
                < std::tie(rhs.count, rhs.query_s, rhs.query_e, rhs.is_revcomp);
        }
    };

    robin_hood::unordered_map<unsigned int, std::vector<Hit>> hits_per_ref;
    std::vector<RescueHit> hits_fw;
    std::vector<RescueHit> hits_rc;
    hits_per_ref.reserve(100);
    hits_fw.reserve(5000);
    hits_rc.reserve(5000);

    for (auto &qr : query_randstrobes) {
        size_t position = index.find(qr.hash);
        if (position != index.end()) {
            unsigned int count = index.get_count(position);
            RescueHit rh{count, position, qr.start, qr.end, qr.is_reverse};
            if (qr.is_reverse){
                hits_rc.push_back(rh);
            } else {
                hits_fw.push_back(rh);
            }
        }
    }

    std::sort(hits_fw.begin(), hits_fw.end());
    std::sort(hits_rc.begin(), hits_rc.end());
    for (auto& rescue_hits : {hits_fw, hits_rc}) {
        int cnt = 0;
        for (auto &rh : rescue_hits) {
            if ((rh.count > filter_cutoff && cnt >= 5) || rh.count > 1000) {
                break;
            }
            add_to_hits_per_ref(hits_per_ref, rh.query_s, rh.query_e, rh.is_revcomp, index, rh.position, 1000);
            cnt++;
        }
    }

    return merge_hits_into_nams(hits_per_ref, index.k(), true);

 */

    todo!()
}

fn add_to_hits_per_ref(
    hits_per_ref: &mut HashMap<usize, Vec<Hit>>,
    query_start: usize,
    query_end: usize,
    is_revcomp: bool,
    index: &StrobemerIndex,
    position: usize,
    max_length_diff: usize,
) {
    let query_length = query_end - query_start;
    let mut max_length_diff = max_length_diff;
    for randstrobe in &index.randstrobes[position..] {
        if randstrobe.hash() != index.randstrobes[position].hash() {
            break;
        }
        let ref_start = randstrobe.position();
        let ref_end = ref_start + randstrobe.strobe2_offset() + index.parameters.syncmer.k;
        let ref_length = ref_end - ref_start;
        let length_diff = (query_length as isize - ref_length as isize).abs() as usize;
        if length_diff <= max_length_diff {
            let hit = Hit{query_start, query_end, ref_start, ref_end, is_revcomp};
            let ref_id = randstrobe.reference_index();
            if hits_per_ref.contains_key(&ref_id) {
                hits_per_ref.get_mut(&ref_id).unwrap().push(hit);
            } else {
                hits_per_ref.insert(ref_id, vec![hit]);
            }

            max_length_diff = length_diff;
        }
    }
}

// TODO should not be mut
fn merge_hits_into_nams(hits_per_ref: &mut HashMap<usize, Vec<Hit>>, k: usize, sort: bool) -> Vec<Nam> {
    let mut nams: Vec<Nam> = Vec::new();
    let mut nam_id_cnt = 0;
    for (ref_id, hits) in hits_per_ref.iter_mut() {
        if sort {
            hits.sort_by_key(|k| (k.query_start, k.ref_start));
        }

        let mut open_nams: Vec<Nam> = Vec::new();
        let mut prev_q_start = 0;
        for h in hits {
            let mut is_added = false;
            for o in &mut open_nams {
                // Extend NAM
                if (o.is_revcomp == h.is_revcomp)
                    && (o.query_prev_hit_startpos < h.query_start)
                    && (h.query_start <= o.query_end)
                    && (o.ref_prev_hit_startpos < h.ref_start)
                    && (h.ref_start <= o.ref_end)
                {
                    if (h.query_end > o.query_end) && (h.ref_end > o.ref_end) {
                        o.query_end = h.query_end;
                        o.ref_end = h.ref_end;
                        o.query_prev_hit_startpos = h.query_start; // log the last strobemer hit in case of outputting paf
                        o.ref_prev_hit_startpos = h.ref_start; // log the last strobemer hit in case of outputting paf
                        o.n_hits += 1;
                        is_added = true;
                        break;
                    }
                    else if (h.query_end <= o.query_end) && (h.ref_end <= o.ref_end) {
                        o.query_prev_hit_startpos = h.query_start; // log the last strobemer hit in case of outputting paf
                        o.ref_prev_hit_startpos = h.ref_start; // log the last strobemer hit in case of outputting paf
                        o.n_hits += 1;
                        is_added = true;
                        break;
                    }
                }

            }
            // Add the hit to open matches
            if !is_added {
                open_nams.push(Nam {
                    nam_id: nam_id_cnt,
                    query_start: h.query_start,
                    query_end: h.query_end,
                    ref_start: h.ref_start,
                    ref_end: h.ref_end,
                    ref_id: *ref_id,
                    query_prev_hit_startpos: h.query_start,
                    ref_prev_hit_startpos: h.ref_start,
                    n_hits: 1,
                    is_revcomp: h.is_revcomp,
                    score: 0,
                });
                nam_id_cnt += 1;
            }

            // Only filter if we have advanced at least k nucleotides
            if h.query_start > prev_q_start + k {

                // Output all NAMs from open_matches to final_nams that the current hit have passed
                for n in &open_nams {
                    if n.query_end < h.query_start {
                        let n_max_span = max(n.query_span(), n.ref_span());
                        let n_min_span = min(n.query_span(), n.ref_span());
                        let n_score =
                            if 2 * n_min_span - n_max_span > 0 {
                                // this is really just n_hits * (min_span - (offset_in_span) ) );
                                n.n_hits * (2 * n_min_span - n_max_span)
                            } else {
                                1
                            } as u32;
//                        n_score = n.n_hits * n.query_span();
                        let mut nam = n.clone();
                        nam.score = n_score;
                        nams.push(nam);
                    }
                }

                // Remove all NAMs from open_matches that the current hit have passed
                let c = h.query_start;

                open_nams = open_nams.into_iter().filter(|&x| x.query_end >= c).collect();

                prev_q_start = h.query_start;
            }
        }

        // Add all current open_matches to final NAMs
        for mut n in open_nams {
            let n_max_span = max(n.query_span(), n.ref_span());
            let n_min_span = min(n.query_span(), n.ref_span());
            let n_score =
                if 2 * n_min_span - n_max_span > 0 {
                    // this is really just n_hits * (min_span - (offset_in_span) ) );
                    n.n_hits * (2 * n_min_span - n_max_span)
                } else {
                    1
                } as u32;
            n.score = n_score;
            nams.push(n);
        }
    }

    nams
}

/// Determine whether the NAM represents a match to the forward or
/// reverse-complemented sequence by checking in which orientation the
/// first and last strobe in the NAM match
///
/// - If first and last strobe match in forward orientation, return true.
/// - If first and last strobe match in reverse orientation, update the NAM
///   in place and return true.
/// - If first and last strobe do not match consistently, return false.
pub fn reverse_nam_if_needed(n: &mut Nam, read: &Read, references: &Vec<RefSequence>, k: usize) -> bool {
    // TODO rename n to nam
    let ref_start_kmer = &references[n.ref_id].sequence[n.ref_start..k];
    let ref_end_kmer = &references[n.ref_id].sequence[n.ref_end - k..k];

    let (seq, seq_rc) = if n.is_revcomp {
        (read.rc(), read.seq())
    } else {
        (read.seq(), read.rc())
    };
    let read_start_kmer = &seq[n.query_start..n.query_start + k];
    let read_end_kmer = &seq[n.query_end - k.. n.query_end];
    if ref_start_kmer == read_start_kmer && ref_end_kmer == read_end_kmer {
        return true;
    }

    // False forward or false reverse (possible due to symmetrical hash values)
    // we need two extra checks for this - hopefully this will remove all the false hits we see
    // (true hash collisions should be very few)
    let read_len = read.len();
    let q_start_tmp = read_len - n.query_end;
    let q_end_tmp = read_len - n.query_start;
    // false reverse hit, change coordinates in nam to forward
    let read_start_kmer = &seq_rc[q_start_tmp..q_start_tmp + k];
    let read_end_kmer = &seq_rc[q_end_tmp - k..q_end_tmp];
    if ref_start_kmer == read_start_kmer && ref_end_kmer == read_end_kmer {
        n.is_revcomp = !n.is_revcomp;
        n.query_start = q_start_tmp;
        n.query_end = q_end_tmp;
        return true;
    }

    false
}
