use crate::fasta::RefSequence;
use crate::strobes::{RandstrobeIterator, RandstrobeParameters, DEFAULT_AUX_LEN};
use crate::syncmers::{SyncmerIterator, SyncmerParameters};
use log::debug;
use std::cmp::{max, min, Reverse};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::Instant;

/// Pre-defined index parameters that work well for a certain
/// "canonical" read length (and similar read lengths)
struct Profile {
    canonical_read_length: usize,
    r_threshold: usize,
    k: usize,
    s_offset: isize,
    l: isize,
    u: isize,
}

static PROFILES: [Profile; 7] = [
    Profile {
        canonical_read_length: 50,
        r_threshold: 70,
        k: 18,
        s_offset: -4,
        l: -2,
        u: 1,
    },
    Profile {
        canonical_read_length: 75,
        r_threshold: 90,
        k: 20,
        s_offset: -4,
        l: -3,
        u: 2,
    },
    Profile {
        canonical_read_length: 100,
        r_threshold: 110,
        k: 20,
        s_offset: -4,
        l: -2,
        u: 2,
    },
    Profile {
        canonical_read_length: 125,
        r_threshold: 135,
        k: 20,
        s_offset: -4,
        l: -1,
        u: 4,
    },
    Profile {
        canonical_read_length: 150,
        r_threshold: 175,
        k: 20,
        s_offset: -4,
        l: 1,
        u: 7,
    },
    Profile {
        canonical_read_length: 250,
        r_threshold: 375,
        k: 22,
        s_offset: -4,
        l: 2,
        u: 12,
    },
    Profile {
        canonical_read_length: 400,
        r_threshold: usize::MAX,
        k: 23,
        s_offset: -6,
        l: 2,
        u: 12,
    },
];

/* Settings that influence index creation */
#[derive(Debug,Clone)]
pub struct IndexParameters {
    pub canonical_read_length: usize,
    pub syncmer: SyncmerParameters,
    pub randstrobe: RandstrobeParameters,
}

impl IndexParameters {
    pub fn new(
        canonical_read_length: usize,
        k: usize,
        s: usize,
        l: isize,
        u: isize,
        q: u64,
        max_dist: u8,
        aux_len: u8,
    ) -> Self {
        let w_min = max(0, (k / (k - s + 1)) as isize + l) as usize;
        let w_max = ((k / (k - s + 1)) as isize + u) as usize;
        let main_hash_mask = !0u64 << (9 + aux_len);
        IndexParameters {
            canonical_read_length,
            syncmer: SyncmerParameters::new(k, s),
            randstrobe: RandstrobeParameters {
                w_min,
                w_max,
                q,
                max_dist,
                main_hash_mask,
            },
        }
    }
    /// Create an IndexParameters instance based on a given read length.
    /// k, s, l, u, c and max_seed_len can be used to override determined parameters
    pub fn from_read_length(
        read_length: usize,
        mut k: Option<usize>,
        mut s: Option<usize>,
        mut l: Option<isize>,
        mut u: Option<isize>,
        c: Option<u32>,
        max_seed_len: Option<usize>,
        aux_len: u8,
    ) -> IndexParameters {
        let default_c = 8;
        let mut canonical_read_length = 50;
        for profile in &PROFILES {
            if read_length <= profile.r_threshold {
                if k.is_none() {
                    k = Some(profile.k);
                }
                if s.is_none() {
                    s = Some((k.unwrap() as isize + profile.s_offset) as usize);
                }
                if l.is_none() {
                    l = Some(profile.l);
                }
                if u.is_none() {
                    u = Some(profile.u);
                }
                canonical_read_length = profile.canonical_read_length;
                break;
            }
        }

        let k = k.unwrap();
        let s = s.unwrap();
        let l = l.unwrap();
        let u = u.unwrap();

        let max_dist = match max_seed_len {
            Some(max_seed_len) => (max_seed_len - k) as u8,
            None => usize::clamp(canonical_read_length.saturating_sub(70), k, 255) as u8,
        };
        let q = 2u64.pow(c.unwrap_or(default_c)) - 1;

        IndexParameters::new(canonical_read_length, k, s, l, u, q, max_dist, aux_len)
    }

    pub fn default_from_read_length(read_length: usize) -> IndexParameters {
        Self::from_read_length(read_length, None, None, None, None, None, None, DEFAULT_AUX_LEN)
    }
}

impl SyncmerParameters {
    /// Pick a suitable number of bits for indexing randstrobe start indices
    fn pick_bits(&self, size: usize) -> u8 {
        let estimated_number_of_randstrobes = size / (self.k - self.s + 1) + 1;
        // Two randstrobes per bucket on average
        // TOOD checked_ilog2 or ilog2
        ((estimated_number_of_randstrobes as f64).log2() as u32).clamp(9, 32) as u8 - 1
    }
}

type RandstrobeHash = u64;
type BucketIndex = u64;

#[derive(Default)]
struct IndexCreationStatistics {
    tot_strobemer_count: u64,
    tot_occur_once: u64,
    tot_high_ab: u64,
    tot_mid_ab: u64,
    index_cutoff: usize,
    distinct_strobemers: u64,
}

#[derive(PartialEq, Eq, PartialOrd, Ord, Default, Clone)]
pub struct RefRandstrobe {
    /// packed representation of the hash and the strobe offset
    hash_offset: u64,
    position: u32,
    ref_index: u32,
}

pub const REF_RANDSTROBE_HASH_MASK: u64 = 0xFFFFFFFFFFFFFF00;
pub const REF_RANDSTROBE_MAX_NUMBER_OF_REFERENCES: usize = u32::MAX as usize;

impl RefRandstrobe {
    fn new(hash: RandstrobeHash, ref_index: u32, position: u32, offset: u8) -> Self {
        let hash_offset = hash & REF_RANDSTROBE_HASH_MASK | (offset as u64);
        RefRandstrobe {
            hash_offset,
            position,
            ref_index,
        }
    }

    pub fn hash(&self) -> RandstrobeHash {
        self.hash_offset & REF_RANDSTROBE_HASH_MASK
    }

    pub fn position(&self) -> usize {
        self.position as usize
    }

    pub fn reference_index(&self) -> usize {
        self.ref_index as usize
    }

    pub fn strobe2_offset(&self) -> usize {
        (self.hash_offset & 0xff) as usize
    }
}

fn count_randstrobes(seq: &[u8], parameters: &IndexParameters) -> usize {
    let syncmer_iterator = SyncmerIterator::new(
        seq,
        parameters.syncmer.k,
        parameters.syncmer.s,
        parameters.syncmer.t,
    );
    let n_syncmers = syncmer_iterator.count();

    // The last w_min syncmers do not result in a randstrobe
    n_syncmers.saturating_sub(parameters.randstrobe.w_min)
}

fn count_all_randstrobes(
    references: &[RefSequence],
    parameters: &IndexParameters,
    n_threads: usize,
) -> Vec<usize> {
    let counts = vec![0; references.len()];
    let mutex = Mutex::new(counts);
    let ref_index = AtomicUsize::new(0);
    thread::scope(|s| {
        for _ in 0..n_threads {
            s.spawn(|| loop {
                let j = ref_index.fetch_add(1, Ordering::SeqCst);
                if j >= references.len() {
                    break;
                }
                let count = count_randstrobes(&references[j].sequence, parameters);
                mutex.lock().unwrap()[j] = count;
            });
        }
    });

    mutex.into_inner().unwrap()
}

// TODO UnpopulatedStrobemerIndex
pub struct StrobemerIndex<'a> {
    references: &'a [RefSequence],
    pub parameters: IndexParameters,
    stats: IndexCreationStatistics,

    /// no. of bits of the hash to use when indexing a randstrobe bucket
    bits: u8,

    /// Regular (non-rescue) NAM finding ignores randstrobes that occur more often than
    /// this (see StrobemerIndex::is_filtered())
    pub filter_cutoff: usize,

    /// Filter partial seeds that occur more often than this
    pub partial_filter_cutoff: usize,

    pub rescue_cutoff: usize,

    /// The randstrobes vector contains all randstrobes sorted by hash.
    /// The randstrobe_start_indices vector points to entries in the
    /// randstrobes vector. `randstrobe_start_indices[x]` is the index of the
    /// first entry in randstrobes whose top *bits* bits of its hash value are
    /// greater than or equal to x.
    ///
    /// randstrobe_start_indices has one extra guard entry at the end that
    /// is always randstrobes.len().
    pub randstrobes: Vec<RefRandstrobe>,
    randstrobe_start_indices: Vec<BucketIndex>,
}

impl<'a> StrobemerIndex<'a> {
    pub fn new(
        references: &'a [RefSequence],
        parameters: IndexParameters,
        bits: Option<u8>,
    ) -> Self {
        let total_reference_length = references.iter().map(|r| r.sequence.len()).sum();
        let bits = bits.unwrap_or_else(|| parameters.syncmer.pick_bits(total_reference_length));
        let randstrobes = vec![];
        let randstrobe_start_indices = vec![];
        let stats = IndexCreationStatistics::default();
        let filter_cutoff = 0;
        let rescue_cutoff = 0;
        StrobemerIndex {
            references,
            parameters,
            stats,
            bits,
            filter_cutoff,
            partial_filter_cutoff: filter_cutoff,
            rescue_cutoff,
            randstrobes,
            randstrobe_start_indices,
        }
    }

    pub fn populate(&mut self, filter_fraction: f64, n_threads: usize) {
        let timer = Instant::now();
        let randstrobe_counts = count_all_randstrobes(self.references, &self.parameters, n_threads);
        debug!("  Counting hashes: {:.2} s", timer.elapsed().as_secs_f64());
        // stats.elapsed_counting_hashes = count_hash.duration();

        let total_randstrobes: usize = randstrobe_counts.iter().sum();
        // stats.tot_strobemer_count = total_randstrobes;

        debug!("  Total number of randstrobes: {}", total_randstrobes);
        let total_length: usize = self
            .references
            .iter()
            .map(|refseq| refseq.sequence.len())
            .sum();
        let memory_bytes: usize = total_length
            + size_of::<RefRandstrobe>() * total_randstrobes
            + size_of::<BucketIndex>() * (1usize << self.bits);
        debug!(
            "  Estimated total memory usage: {:.1} GB",
            memory_bytes as f64 / 1E9
        );

        if total_randstrobes > BucketIndex::MAX as usize {
            panic!("Too many randstrobes");
        }
        let timer = Instant::now();
        debug!("  Generating randstrobes ...");
        self.randstrobes = self.make_randstrobes_parallel(&randstrobe_counts, n_threads);
        debug!("  Generating seeds: {:.2} s", timer.elapsed().as_secs_f64());
        // stats.elapsed_generating_seeds = randstrobes_timer.duration();

        let timer = Instant::now();
        debug!("  Sorting ...");
        // TODO
        // ensure comparison function is branchless
        // Comment from C++ code:
        // Compare both hash and position to ensure that the order of the
        // RefRandstrobes in the index is reproducible no matter which sorting
        // function is used. This branchless comparison is faster than the
        // equivalent one using std::tie.
        // __uint128_t lhs = (static_cast<__uint128_t>(m_hash_offset_flag) << 64) | ((static_cast<uint64_t>(m_position) << 32) | m_ref_index);
        // __uint128_t rhs = (static_cast<__uint128_t>(other.m_hash_offset_flag) << 64) | ((static_cast<uint64_t>(other.m_position) << 32) | m_ref_index);
        self.randstrobes.sort_unstable_by_key(|r| (r.hash_offset, r.position, r.ref_index));
        debug!("    Took {:.2} s", timer.elapsed().as_secs_f64());
        // stats.elapsed_sorting_seeds = sorting_timer.duration();

        let timer = Instant::now();
        debug!("  Generating hash table index ...");

        let mut tot_high_ab = 0;
        let mut tot_mid_ab = 0;
        let mut strobemer_counts = Vec::<usize>::new();

        self.stats.tot_occur_once = 0;
        self.randstrobe_start_indices
            .reserve((1usize << self.bits) + 1);

        let mut unique_mers = u64::from(!self.randstrobes.is_empty());

        let mut prev_hash: RandstrobeHash = if self.randstrobes.is_empty() {
            0
        } else {
            self.randstrobes[0].hash()
        };
        let mut count = 1;

        if !self.randstrobes.is_empty() {
            self.randstrobe_start_indices.push(0);
        }
        for position in 1..self.randstrobes.len() {
            let cur_hash = self.randstrobes[position].hash();
            if cur_hash == prev_hash {
                count += 1;
                continue;
            }
            unique_mers += 1;

            if count == 1 {
                self.stats.tot_occur_once += 1;
            } else {
                if count > 100 {
                    tot_high_ab += 1;
                } else {
                    tot_mid_ab += 1;
                }
                strobemer_counts.push(count);
            }
            count = 1;
            let cur_hash_n = cur_hash >> (64 - self.bits);
            while self.randstrobe_start_indices.len() <= cur_hash_n as usize {
                self.randstrobe_start_indices
                    .push(position as BucketIndex);
            }
            prev_hash = cur_hash;
        }
        // wrap up last entry
        if count == 1 {
            self.stats.tot_occur_once += 1;
        } else {
            if count > 100 {
                tot_high_ab += 1;
            } else {
                tot_mid_ab += 1;
            }
            strobemer_counts.push(count);
        }
        while self.randstrobe_start_indices.len() < ((1usize << self.bits) + 1) {
            self.randstrobe_start_indices
                .push(self.randstrobes.len() as BucketIndex);
        }
        self.stats.tot_high_ab = tot_high_ab;
        self.stats.tot_mid_ab = tot_mid_ab;

        strobemer_counts.sort_unstable_by_key(|k| Reverse(*k));

        let index_cutoff = (unique_mers as f64 * filter_fraction) as usize;
        self.stats.index_cutoff = index_cutoff;
        self.filter_cutoff = usize::clamp(
            if index_cutoff < strobemer_counts.len() {
                strobemer_counts[index_cutoff]
            } else {
                *strobemer_counts.last().unwrap_or(&30)
            },
            30, // cutoff is around 30-50 on hg38. No reason to have a lower cutoff than this if aligning to a smaller genome or contigs.
            100, // limit upper cutoff for normal NAM finding - use rescue mode instead
        );
        self.rescue_cutoff = min(self.filter_cutoff * 2, 1000);
        //stats.elapsed_hash_index = hash_index_timer.duration();
        debug!("    Took {:.2} s", timer.elapsed().as_secs_f64());
        self.stats.distinct_strobemers = unique_mers;
    }

/*
    fn make_randstrobes(&self, randstrobe_counts: &[usize]) -> Vec<RefRandstrobe> {
        let mut randstrobes = vec![RefRandstrobe::default(); randstrobe_counts.iter().sum()];

        // Fill randstrobes vector
        let mut offset = 0;
        for ref_index in 0..self.references.len() {
            self.assign_randstrobes(ref_index, &mut randstrobes[offset..offset+randstrobe_counts[ref_index]]);
            offset += randstrobe_counts[ref_index];
        }
        randstrobes
    }
*/
    
    fn make_randstrobes_parallel(&mut self, randstrobe_counts: &[usize], n_threads: usize) -> Vec<RefRandstrobe> {
        let mut randstrobes = vec![RefRandstrobe::default(); randstrobe_counts.iter().sum()];
        let mut slices = vec![];
        {
            let mut slice = &mut randstrobes[..];
            for ref_index in 0..self.references.len()-1 {
                let (left, right) = slice.split_at_mut(randstrobe_counts[ref_index]);
                slices.push(Arc::new(Mutex::new(left)));
                slice = right;
            }
            slices.push(Arc::new(Mutex::new(slice)));
        }
        let ref_index = AtomicUsize::new(0);
        thread::scope(|s| {
            for _ in 0..n_threads {
                s.spawn(|| loop {
                    let j = ref_index.fetch_add(1, Ordering::SeqCst);
                    if j >= self.references.len() {
                        break;
                    }
                    self.assign_randstrobes(j, *slices[j].lock().unwrap());
                });
            }
        });

        randstrobes
    }

    /// Compute randstrobes of one reference contig and assign them to the provided slice
    fn assign_randstrobes(&self, ref_index: usize, randstrobes: &mut [RefRandstrobe]) {
        let seq = &self.references[ref_index].sequence;
        if seq.len() < self.parameters.randstrobe.w_max {
            return;
        }
        let mut syncmer_iter = SyncmerIterator::new(
            seq,
            self.parameters.syncmer.k,
            self.parameters.syncmer.s,
            self.parameters.syncmer.t,
        );
        let randstrobe_iter =
            RandstrobeIterator::new(&mut syncmer_iter, &self.parameters.randstrobe);

        let mut n= 0;
        for (i, randstrobe) in randstrobe_iter.enumerate() {
            n += 1;
            let offset = randstrobe.strobe2_pos - randstrobe.strobe1_pos;
            randstrobes[i] =
                RefRandstrobe::new(
                    randstrobe.hash,
                    ref_index as u32,
                    randstrobe.strobe1_pos as u32,
                    offset as u8,
                );

        }
        debug_assert_eq!(n, randstrobes.len());
    }

    pub fn get_full(&self, hash: RandstrobeHash) -> Option<usize> {
        self.get_masked(hash, REF_RANDSTROBE_HASH_MASK)
    }

    /// Find the first entry that matches the main hash
    pub fn get_partial(&self, hash: RandstrobeHash) -> Option<usize> {
        self.get_masked(hash, self.parameters.randstrobe.main_hash_mask)
    }

    /// Find index of first entry in randstrobe table that has the given
    /// hash value masked by the `hash_mask`
    pub fn get_masked(&self, hash: RandstrobeHash, hash_mask: RandstrobeHash) -> Option<usize> {
        let masked_hash = hash & hash_mask;
        const MAX_LINEAR_SEARCH: usize = 4;
        let top_n = (hash >> (64 - self.bits)) as usize;
        let position_start = self.randstrobe_start_indices[top_n];
        let position_end = self.randstrobe_start_indices[top_n + 1];
        let bucket = &self.randstrobes[position_start as usize..position_end as usize];
        if bucket.is_empty() {
            return None;
        } else if bucket.len() < MAX_LINEAR_SEARCH {
            for (pos, randstrobe) in bucket.iter().enumerate() {
                if randstrobe.hash() & hash_mask == masked_hash {
                    return Some(position_start as usize + pos);
                }
                if randstrobe.hash() & hash_mask > masked_hash {
                    return None;
                }
            }
            return None;
        }

        let pos = bucket.partition_point(|h| h.hash() & hash_mask < masked_hash);
        if pos < bucket.len() && bucket[pos].hash() & hash_mask == masked_hash {
            Some(position_start as usize + pos)
        } else {
            None
        }
    }

    pub fn k(&self) -> usize {
        self.parameters.syncmer.k
    }

    pub fn get_hash_partial(&self, position: usize) -> RandstrobeHash {
        self.randstrobes[position].hash() & self.parameters.randstrobe.main_hash_mask
    }

    pub fn strobe_extent_partial(&self, position: usize) -> (usize, usize) {
        let p = self.randstrobes[position].position;

        (p as usize, p as usize + self.k())
    }

    pub fn get_count_full(&self, position: usize) -> usize {
        self.get_count(position, REF_RANDSTROBE_HASH_MASK)
    }

    pub fn get_count_partial(&self, position: usize) -> usize {
        self.get_count(position, self.parameters.randstrobe.main_hash_mask)
    }

    pub fn get_count(&self, position: usize, hash_mask: u64) -> usize {
        const MAX_LINEAR_SEARCH: usize = 8;
        let key = self.randstrobes[position].hash();
        let masked_key = key & hash_mask;
        let top_n = (key >> (64 - self.bits)) as usize;
        let position_end = self.randstrobe_start_indices[top_n + 1] as usize;

        if position_end - position < MAX_LINEAR_SEARCH {
            let mut count = 1;
            for position_start in position + 1..position_end {
                if self.randstrobes[position_start].hash() & hash_mask == masked_key {
                    count += 1;
                } else {
                    break;
                }
            }
            count
        } else {
            let bucket = &self.randstrobes[position..position_end];
            bucket.partition_point(|h| h.hash() & hash_mask <= masked_key)
        }
    }

    /// Return whether the randstrobe at the given position occurs more often than cutoff
    pub fn is_too_frequent(&self, position: usize, cutoff: usize) -> bool {
        if position + self.filter_cutoff < self.randstrobes.len() {
            self.randstrobes[position].hash() == self.randstrobes[position + cutoff].hash()
        } else {
            false
        }
    }

    pub fn is_too_frequent_partial(&self, position: usize, cutoff: usize) -> bool {
        if position + cutoff < self.randstrobes.len() {
            self.randstrobes[position].hash() & self.parameters.randstrobe.main_hash_mask == self.randstrobes[position + cutoff].hash() & self.parameters.randstrobe.main_hash_mask
        } else {
            false
        }
    }
}

#[cfg(test)]
mod tests {
    use std::fs::File;
    use std::io::BufReader;
    use crate::fasta::read_fasta;
    use crate::revcomp::reverse_complement;
    use crate::syncmers::Syncmer;
    use super::*;

    #[test]
    fn test_ref_randstrobe() {
        let hash: u64 = 0x1234567890ABCDEFu64 & REF_RANDSTROBE_HASH_MASK;
        let ref_index: u32 = (REF_RANDSTROBE_MAX_NUMBER_OF_REFERENCES - 1) as u32;
        let offset = 255;
        let position = !0;
        let rr = RefRandstrobe::new(hash, ref_index, position, offset);

        assert_eq!(rr.hash(), hash);
        assert_eq!(rr.position(), position as usize);
        assert_eq!(rr.reference_index(), ref_index as usize);
        assert_eq!(rr.strobe2_offset(), offset as usize);
    }

    #[test]
    fn test_ref_randstrobe2() {
        let hash: u64 = 0x1234567890ABCDEFu64 & REF_RANDSTROBE_HASH_MASK;
        let ref_index: u32 = (REF_RANDSTROBE_MAX_NUMBER_OF_REFERENCES - 1) as u32;
        let offset = 255;
        let position = !0;
        let rr = RefRandstrobe::new(hash, ref_index, position, offset);

        assert_eq!(rr.hash(), hash);
        assert_eq!(rr.position(), position as usize);
        assert_eq!(rr.reference_index(), ref_index as usize);
        assert_eq!(rr.strobe2_offset(), offset as usize);
    }

    #[test]
    fn test_index_parameters() {
        let canonical_read_length = 250;
        let k = 22;
        let s = 18;
        let l = 2;
        let u = 12;
        let w_min = 6;
        let w_max = 16;
        let max_dist = 180;
        let q = 255;
        let main_hash_mask = 0xfffffffffc000000;
        let aux_len = 17;
        let sp = SyncmerParameters::new(k, s);
        let rp = RandstrobeParameters {w_min, w_max, q, max_dist, main_hash_mask};
        let ip = IndexParameters::new(
            canonical_read_length,
            k,
            s,
            l,
            u,
            q,
            max_dist,
            aux_len,
        );
        assert_eq!(ip.canonical_read_length, canonical_read_length);
        assert_eq!(ip.randstrobe, rp);
        assert_eq!(ip.syncmer, sp);

        let ip = IndexParameters::default_from_read_length(canonical_read_length + 1);
        assert_eq!(ip.canonical_read_length, canonical_read_length);
        assert_eq!(ip.randstrobe, rp);
        assert_eq!(ip.syncmer, sp);
    }
    
    fn syncmers_of(seq: &[u8], parameters: &SyncmerParameters) -> Vec<Syncmer> {
        SyncmerIterator::new(seq, parameters.k, parameters.s, parameters.t).collect()
    }

    #[test]
    fn test_canonical_syncmers() {
        let parameters = SyncmerParameters::new(20, 16);
        let f = File::open("tests/phix.fasta").unwrap();
        let mut reader = BufReader::new(f);
        let records = read_fasta(&mut reader).unwrap();
        let seqs = vec!["AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".as_bytes(), &records[0].sequence];
        for seq in seqs {
            let seq_revcomp = reverse_complement(seq);
            let syncmers_forward = syncmers_of(seq, &parameters);
            let mut syncmers_reverse = syncmers_of(&seq_revcomp, &parameters);
            syncmers_reverse.reverse();
            for syncmer_rev in &mut syncmers_reverse {
                syncmer_rev.position = seq.len() - parameters.k - syncmer_rev.position;
            }
            assert_eq!(syncmers_forward, syncmers_reverse);
        }
    }
    
    #[test]
    fn test_pick_bits() {
        let parameters = SyncmerParameters::new(20, 16);
        assert_eq!(parameters.pick_bits(0), 8);
    }
}
