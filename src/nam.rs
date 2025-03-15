use std::cmp::{max, min};
use std::collections::hash_map::Entry;
use std::collections::HashMap;
use std::fmt::{Display, Formatter};
use std::hash::{BuildHasherDefault, DefaultHasher};
use std::time::Instant;
use fastrand::Rng;
use log::Level::Trace;
use log::trace;
use crate::details::NamDetails;
use crate::fasta::RefSequence;
use crate::index::StrobemerIndex;
use crate::mapper;
use crate::mapper::QueryRandstrobe;
use crate::read::Read;

/// Non-overlapping approximate match
#[derive(Clone, Debug)]
pub struct Nam {
    pub nam_id: usize,
    pub ref_start: usize,
    pub ref_end: usize,
    pub query_start: usize,
    pub query_end: usize,
    query_prev_match_startpos: usize,
    ref_prev_match_startpos: usize,
    pub n_matches: usize,
    pub ref_id: usize,
    pub score: u32,
    pub is_revcomp: bool,
}

impl Nam {
    pub fn ref_span(&self) -> usize {
        self.ref_end - self.ref_start
    }

    pub fn query_span(&self) -> usize {
        self.query_end - self.query_start
    }

    pub fn projected_ref_start(&self) -> usize {
        self.ref_start.saturating_sub(self.query_start)
    }
}

impl Display for Nam {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "Nam(ref_id={}, query: {}..{}, ref: {}..{}, rc={}, score={})",
               self.ref_id, self.query_start, self.query_end, self.ref_start, self.ref_end, self.is_revcomp as u8, self.score
        )?;
        Ok(())
    }
}

struct Match {
    query_start: usize,
    query_end: usize,
    ref_start: usize,
    ref_end: usize,
}

// The regular HashMap uses a randomly seeded hash function, but we want reproducible results
type DefaultHashMap<K, V> = HashMap<K, V, BuildHasherDefault<DefaultHasher>>;

/// Find a query’s NAMs, ignoring randstrobes that occur too often in the
/// reference (have a count above filter_cutoff).
///
/// Return the fraction of nonrepetitive hits (those not above the filter_cutoff threshold)
///
pub fn find_nams(query_randstrobes: &Vec<QueryRandstrobe>, index: &StrobemerIndex, filter_cutoff: usize) -> (f32, usize, Vec<Nam>) {
   
    let mut matches_map = [DefaultHashMap::default(), DefaultHashMap::default()];
    matches_map[0].reserve(100);
    matches_map[1].reserve(100);
    let mut nr_good_hits = 0;
    let mut total_hits = 0;
    for randstrobe in query_randstrobes {
        if let Some(position) = index.get(randstrobe.hash) {
            total_hits += 1;
            if index.is_too_frequent(position, filter_cutoff) {
                continue;
            }
            nr_good_hits += 1;
            add_to_matches_map(&mut matches_map[randstrobe.is_revcomp as usize], randstrobe.start, randstrobe.end, index, position);
        }
    }
    let nonrepetitive_fraction = if total_hits > 0 { (nr_good_hits as f32) / (total_hits as f32) } else { 1.0 };
    let nams = merge_matches_into_nams_forward_and_reverse(&mut matches_map, index.parameters.syncmer.k, false);

    (nonrepetitive_fraction, nr_good_hits, nams)
}

/// Find a query’s NAMs, using also some of the randstrobes that occur more often
/// than the normal (non-rescue) filter cutoff.
///
/// Return the number of hits and the vector of NAMs
pub fn find_nams_rescue(
    query_randstrobes: &[QueryRandstrobe],
    index: &StrobemerIndex,
    rescue_cutoff: usize,
) -> (usize, Vec<Nam>) {

    struct RescueHit {
        count: usize,
        position: usize,
        query_start: usize,
        query_end: usize,
    }
/*
        bool operator< (const RescueHit& rhs) const {
            return std::tie(count, query_s, query_e, is_revcomp)
                < std::tie(rhs.count, rhs.query_s, rhs.query_e, rhs.is_revcomp);
        }
    }*/

    let mut hits_fw = Vec::with_capacity(5000);
    let mut hits_rc = Vec::with_capacity(5000);

    for randstrobe in query_randstrobes {

        if let Some(position) = index.get(randstrobe.hash) {
            let count = index.get_count(position);
            let rh = RescueHit {
                count,
                position,
                query_start: randstrobe.start,
                query_end: randstrobe.end,
            };
            if randstrobe.is_revcomp {
                hits_rc.push(rh);
            } else {
                hits_fw.push(rh);
            }
        }
    }

    let cmp = |a: &RescueHit, b: &RescueHit| (a.count, a.query_start, a.query_end).cmp(&(b.count, b.query_start, b.query_end));
    hits_fw.sort_by(cmp);
    hits_rc.sort_by(cmp);

    let mut matches_map = [DefaultHashMap::default(), DefaultHashMap::default()];
    matches_map[0].reserve(100);
    matches_map[1].reserve(100);

    let mut n_hits = 0;
    for (is_revcomp, rescue_hits) in [(false, hits_fw), (true, hits_rc)] {
        for (i, rh) in rescue_hits.iter().enumerate() {
            if (rh.count > rescue_cutoff && i >= 5) || rh.count > 1000 {
                break;
            }
            add_to_matches_map(&mut matches_map[is_revcomp as usize], rh.query_start, rh.query_end, index, rh.position);
            n_hits += 1;
        }
    }

    (n_hits, merge_matches_into_nams_forward_and_reverse(&mut matches_map, index.parameters.syncmer.k, true))
}

fn add_to_matches_map(
    matches_map: &mut DefaultHashMap<usize, Vec<Match>>,
    query_start: usize,
    query_end: usize,
    index: &StrobemerIndex,
    position: usize,
) {
    let query_length = query_end - query_start;
    let mut max_length_diff = usize::MAX;
    for randstrobe in &index.randstrobes[position..] {
        if randstrobe.hash() != index.randstrobes[position].hash() {
            break;
        }
        let ref_start = randstrobe.position();
        let ref_end = ref_start + randstrobe.strobe2_offset() + index.parameters.syncmer.k;
        let ref_length = ref_end - ref_start;
        let length_diff = (query_length as isize - ref_length as isize).unsigned_abs();
        if length_diff <= max_length_diff {
            let match_ = Match {query_start, query_end, ref_start, ref_end};
            let ref_id = randstrobe.reference_index();

            if let Entry::Vacant(e) = matches_map.entry(ref_id) {
                e.insert(vec![match_]);
            } else {
                matches_map.get_mut(&ref_id).unwrap().push(match_);
            }
            max_length_diff = length_diff;
        }
    }
}

// TODO should not be mut
fn merge_matches_into_nams(matches_map: &mut DefaultHashMap<usize, Vec<Match>>, k: usize, sort: bool, is_revcomp: bool, nams: &mut Vec<Nam>) {
    for (ref_id, matches) in matches_map.iter_mut() {
        if sort {
            matches.sort_by_key(|k| (k.query_start, k.ref_start));
        }

        let mut open_nams: Vec<Nam> = Vec::new();
        let mut prev_q_start = 0;
        for m in matches {
            let mut is_added = false;
            for o in &mut open_nams {
                // Extend NAM
                if (o.query_prev_match_startpos < m.query_start)
                    && (m.query_start <= o.query_end)
                    && (o.ref_prev_match_startpos < m.ref_start)
                    && (m.ref_start <= o.ref_end)
                {
                    if (m.query_end > o.query_end) && (m.ref_end > o.ref_end) {
                        o.query_end = m.query_end;
                        o.ref_end = m.ref_end;
                        o.query_prev_match_startpos = m.query_start; // log the last strobemer hit in case of outputting paf
                        o.ref_prev_match_startpos = m.ref_start; // log the last strobemer hit in case of outputting paf
                        o.n_matches += 1;
                        is_added = true;
                        break;
                    } else if (m.query_end <= o.query_end) && (m.ref_end <= o.ref_end) {
                        o.query_prev_match_startpos = m.query_start; // log the last strobemer hit in case of outputting paf
                        o.ref_prev_match_startpos = m.ref_start; // log the last strobemer hit in case of outputting paf
                        o.n_matches += 1;
                        is_added = true;
                        break;
                    }
                }
            }
            // Add the match to open matches
            if !is_added {
                open_nams.push(Nam {
                    nam_id: nams.len() + open_nams.len(),
                    query_start: m.query_start,
                    query_end: m.query_end,
                    ref_start: m.ref_start,
                    ref_end: m.ref_end,
                    ref_id: *ref_id,
                    query_prev_match_startpos: m.query_start,
                    ref_prev_match_startpos: m.ref_start,
                    n_matches: 1,
                    is_revcomp,
                    score: 0,
                });
            }

            // Only filter if we have advanced at least k nucleotides
            if m.query_start > prev_q_start + k {

                // Output all NAMs from open_matches to final_nams that the current match have passed
                for n in &open_nams {
                    if n.query_end < m.query_start {
                        let n_max_span = max(n.query_span(), n.ref_span());
                        let n_min_span = min(n.query_span(), n.ref_span());
                        let n_score =
                            if 2 * n_min_span > n_max_span {
                                // this is really just n_matches * (min_span - (offset_in_span) ) );
                                n.n_matches * (2 * n_min_span - n_max_span)
                            } else {
                                1
                            } as u32;
//                        n_score = n.n_matches * n.query_span();
                        let mut nam = n.clone();
                        nam.score = n_score;
                        nams.push(nam);
                    }
                }

                // Remove all NAMs from open_matches that the current match have passed
                let c = m.query_start;

                open_nams.retain(|x| x.query_end >= c);

                prev_q_start = m.query_start;
            }
        }

        // Add all current open_matches to final NAMs
        for mut n in open_nams {
            let n_max_span = max(n.query_span(), n.ref_span());
            let n_min_span = min(n.query_span(), n.ref_span());
            n.score =
                if 2 * n_min_span > n_max_span {
                    n.n_matches * (2 * n_min_span - n_max_span)
                } else {
                    1
                } as u32;
            nams.push(n);
        }
    }
}

fn merge_matches_into_nams_forward_and_reverse(
    matches_map: &mut [DefaultHashMap<usize, Vec<Match>>; 2],
    k: usize,
    sort: bool
) -> Vec<Nam> {
    let mut nams = Vec::new();
    for is_revcomp in [false, true] {
        merge_matches_into_nams(&mut matches_map[is_revcomp as usize], k, sort, is_revcomp, &mut nams);
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
pub fn reverse_nam_if_needed(nam: &mut Nam, read: &Read, references: &[RefSequence], k: usize) -> bool {
    let ref_start_kmer = &references[nam.ref_id].sequence[nam.ref_start..nam.ref_start + k];
    let ref_end_kmer = &references[nam.ref_id].sequence[nam.ref_end - k..nam.ref_end];

    let (seq, seq_rc) = if nam.is_revcomp {
        (read.rc(), read.seq())
    } else {
        (read.seq(), read.rc())
    };
    let read_start_kmer = &seq[nam.query_start..nam.query_start + k];
    let read_end_kmer = &seq[nam.query_end - k.. nam.query_end];
    if ref_start_kmer == read_start_kmer && ref_end_kmer == read_end_kmer {
        return true;
    }

    // False forward or false reverse (possible due to symmetrical hash values)
    // we need two extra checks for this - hopefully this will remove all the false matches we see
    // (true hash collisions should be very few)
    let read_len = read.len();
    let q_start_tmp = read_len - nam.query_end;
    let q_end_tmp = read_len - nam.query_start;
    // false reverse match, change coordinates in nam to forward
    let read_start_kmer = &seq_rc[q_start_tmp..q_start_tmp + k];
    let read_end_kmer = &seq_rc[q_end_tmp - k..q_end_tmp];
    if ref_start_kmer == read_start_kmer && ref_end_kmer == read_end_kmer {
        nam.is_revcomp = !nam.is_revcomp;
        nam.query_start = q_start_tmp;
        nam.query_end = q_end_tmp;
        true
    } else {
        false
    }
}

// TODO rename
pub fn get_nams(sequence: &[u8], index: &StrobemerIndex, rescue_level: usize, rng: &mut Rng) -> (NamDetails, Vec<Nam>) {
    let timer = Instant::now();
    let query_randstrobes = mapper::randstrobes_query(sequence, &index.parameters);
    let n_randstrobes = query_randstrobes.len();
    let time_randstrobes = timer.elapsed().as_secs_f64();

    let timer = Instant::now();
    let (nonrepetitive_fraction, n_hits, mut nams) = find_nams(&query_randstrobes, index, index.filter_cutoff);
    let n_nams = nams.len();
    let time_find_nams = timer.elapsed().as_secs_f64();

    let mut n_rescue_hits = 0;
    let mut n_rescue_nams = 0;
    let mut nam_rescue = false;
    let time_rescue;
    if rescue_level > 1 {
        let timer = Instant::now();
        if nams.is_empty() || nonrepetitive_fraction < 0.7 {
            nam_rescue = true;
            (n_rescue_hits, nams) = find_nams_rescue(&query_randstrobes, index, index.rescue_cutoff);
            n_rescue_nams = nams.len();
        }
        time_rescue = timer.elapsed().as_secs_f64();
    } else {
        time_rescue = 0f64;
    }

    let timer = Instant::now();
    nams.sort_by_key(|k| -(k.score as i32));
    shuffle_top_nams(&mut nams, rng);
    let time_sort_nams = timer.elapsed().as_secs_f64();

    if log::log_enabled!(Trace) {
        trace!("Found {} NAMs (rescue done: {})", nams.len(), nam_rescue);
        for nam in &nams {
            trace!("- {}", nam);
        }
    }

    let nam_details = NamDetails {
        n_reads: 1,
        n_randstrobes,
        n_nams,
        n_rescue_nams,
        nam_rescue: nam_rescue as usize,
        n_hits,
        n_rescue_hits,
        time_randstrobes,
        time_find_nams,
        time_rescue,
        time_sort_nams,
    };

    (nam_details, nams)
}

/// Shuffle the top-scoring NAMs. Input must be sorted by score.
/// This helps to ensure we pick a random location in case there are multiple
/// equally good ones.
fn shuffle_top_nams(nams: &mut [Nam], rng: &mut Rng) {
    if let Some(best) = nams.first() {
        let best_score = best.score;

        let pos = nams.iter().position(|nam| nam.score != best_score);
        let end = pos.unwrap_or(nams.len());
        if end > 1 {
            rng.shuffle(&mut nams[0..end]);
        }
    }
}
