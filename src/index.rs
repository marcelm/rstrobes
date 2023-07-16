use std::cmp::max;
use rayon::prelude::*;
use crate::strobes::{RandstrobeParameters, SyncmerIterator, SyncmerParameters};
use crate::fasta::RefSequence;

/* Settings that influence index creation */
pub struct IndexParameters {
    canonical_read_length: usize,
    pub syncmer: SyncmerParameters,
    pub randstrobe: RandstrobeParameters,
}

impl IndexParameters {
    pub fn new(canonical_read_length: usize, k: usize, s: usize, l: usize, u: usize, q: u64, max_dist: u8) -> Self {
        let w_min = max(1, k / (k - s + 1) + l);
        let w_max = k / (k - s + 1) + u;
        IndexParameters {
            canonical_read_length,
            syncmer: SyncmerParameters::new(k, s),
            randstrobe: RandstrobeParameters  { w_min, w_max, q, max_dist }
        }
    }
}

impl SyncmerParameters {
    /// Pick a suitable number of bits for indexing randstrobe start indices
    pub fn pick_bits(&self, size: usize) -> u8 {
        let estimated_number_of_randstrobes = size / (self.k - self.s + 1) as usize;
        // Two randstrobes per bucket on average
        // TOOD checked_ilog2 or ilog2
        ((estimated_number_of_randstrobes as f64).log2() as u32 - 1).clamp(8, 31) as u8
    }
}

fn count_randstrobes(seq: &[u8], parameters: &IndexParameters) -> usize {
    let syncmer_iterator = SyncmerIterator::new(seq, parameters.syncmer.k, parameters.syncmer.s, parameters.syncmer.t);
    let n_syncmers = syncmer_iterator.count();

    // The last w_min syncmers do not result in a randstrobe
    if n_syncmers < parameters.randstrobe.w_min {
        0
    } else {
        n_syncmers - parameters.randstrobe.w_min
    }
}

fn count_all_randstrobes(references: &Vec<RefSequence>, parameters: &IndexParameters) -> Vec<usize> {
    references.par_iter().map(|refseq| count_randstrobes(&refseq.sequence, parameters)).collect()
}

pub struct StrobemerIndex<'a> {
    references: &'a Vec<RefSequence>,
    parameters: IndexParameters,
    bits: u8,
}

impl<'a> StrobemerIndex<'a> {
    pub fn new(references: &'a Vec<RefSequence>, parameters: IndexParameters, bits: Option<u8>) -> Self {
        let total_reference_length = references.iter().map(|r| r.sequence.len()).sum();
        let bits = bits.unwrap_or_else(|| parameters.syncmer.pick_bits(total_reference_length));

        StrobemerIndex { references, parameters, bits }
    }

    pub fn populate(&self, filter_fraction: f32) {

    }
}
/*
void StrobemerIndex::populate(float f, size_t n_threads) {
    Timer count_hash;
    auto randstrobe_counts = count_all_randstrobes(references, parameters, n_threads);
    stats.elapsed_counting_hashes = count_hash.duration();

    uint64_t total_randstrobes = 0;
    for (auto& count : randstrobe_counts) {
        total_randstrobes += count;
    }
    stats.tot_strobemer_count = total_randstrobes;

    logger.debug() << "  Total number of randstrobes: " << total_randstrobes << '\n';
    uint64_t memory_bytes = references.total_length() + sizeof(RefRandstrobe) * total_randstrobes + sizeof(bucket_index_t) * (1u << bits);
    logger.debug() << "  Estimated total memory usage: " << memory_bytes / 1E9 << " GB\n";

    if (total_randstrobes > std::numeric_limits<bucket_index_t>::max()) {
        throw std::range_error("Too many randstrobes");
    }
    Timer randstrobes_timer;
    logger.debug() << "  Generating randstrobes ...\n";
    randstrobes.assign(total_randstrobes, RefRandstrobe{0, 0, 0});
    assign_all_randstrobes(randstrobe_counts, n_threads);
    stats.elapsed_generating_seeds = randstrobes_timer.duration();

    Timer sorting_timer;
    logger.debug() << "  Sorting ...\n";
    // sort by hash values
    pdqsort_branchless(randstrobes.begin(), randstrobes.end());
    stats.elapsed_sorting_seeds = sorting_timer.duration();

    Timer hash_index_timer;
    logger.debug() << "  Indexing ...\n";

    uint64_t tot_high_ab = 0;
    uint64_t tot_mid_ab = 0;
    std::vector<uint64_t> strobemer_counts;

    stats.tot_occur_once = 0;
    randstrobe_start_indices.reserve((1u << bits) + 1);

    uint64_t unique_mers = randstrobes.empty() ? 0 : 1;
    randstrobe_hash_t prev_hash = randstrobes.empty() ? 0 : randstrobes[0].hash;
    unsigned int count = 0;

    for (bucket_index_t position = 0; position < randstrobes.size(); ++position) {
        const randstrobe_hash_t cur_hash = randstrobes[position].hash;
        if (cur_hash == prev_hash) {
            ++count;
            continue;
        }

        ++unique_mers;

        if (count == 1) {
            ++stats.tot_occur_once;
        } else {
            if (count > 100) {
                ++tot_high_ab;
            } else {
                ++tot_mid_ab;
            }
            strobemer_counts.push_back(count);
        }
        count = 1;
        const unsigned int cur_hash_N = cur_hash >> (64 - bits);
        while (randstrobe_start_indices.size() <= cur_hash_N) {
            randstrobe_start_indices.push_back(position);
        }
        prev_hash = cur_hash;
    }
    // wrap up last entry
    if (count == 1) {
        ++stats.tot_occur_once;
    } else {
        if (count > 100) {
            tot_high_ab++;
        } else {
            tot_mid_ab++;
        }
        strobemer_counts.push_back(count);
    }
    while (randstrobe_start_indices.size() < ((1u << bits) + 1)) {
        randstrobe_start_indices.push_back(randstrobes.size());
    }
    stats.tot_high_ab = tot_high_ab;
    stats.tot_mid_ab = tot_mid_ab;

    std::sort(strobemer_counts.begin(), strobemer_counts.end(), std::greater<int>());

    uint64_t index_cutoff = unique_mers * f;
    stats.index_cutoff = index_cutoff;
    if (!strobemer_counts.empty()){
        filter_cutoff = index_cutoff < strobemer_counts.size() ?  strobemer_counts[index_cutoff] : strobemer_counts.back();
        filter_cutoff = std::max(30U, filter_cutoff); // cutoff is around 30-50 on hg38. No reason to have a lower cutoff than this if aligning to a smaller genome or contigs.
        filter_cutoff = std::min(100U, filter_cutoff); // limit upper cutoff for normal NAM finding - use rescue mode instead
    } else {
        filter_cutoff = 30;
    }
    stats.filter_cutoff = filter_cutoff;
    stats.elapsed_hash_index = hash_index_timer.duration();
    stats.distinct_strobemers = unique_mers;
}

void StrobemerIndex::assign_all_randstrobes(const std::vector<uint64_t>& randstrobe_counts, size_t n_threads) {
    // Compute offsets
    std::vector<size_t> offsets;
    size_t offset = 0;
    for (size_t ref_index = 0; ref_index < references.size(); ++ref_index) {
        offsets.push_back(offset);
        offset += randstrobe_counts[ref_index];
    }

    std::vector<std::thread> workers;
    std::atomic_size_t ref_index{0};
    for (size_t i = 0; i < n_threads; ++i) {
        workers.push_back(
            std::thread(
                [&]() {
                    while (true) {
                        size_t j = ref_index.fetch_add(1);
                        if (j >= references.size()) {
                            break;
                        }
                        assign_randstrobes(j, offsets[j]);
                    }
                })
        );
    }
    for (auto& worker : workers) {
        worker.join();
    }
}

/*
 * Compute randstrobes of one reference and assign them to the randstrobes
 * vector starting from the given offset
 */
void StrobemerIndex::assign_randstrobes(size_t ref_index, size_t offset) {
    auto seq = references.sequences[ref_index];
    if (seq.length() < parameters.randstrobe.w_max) {
        return;
    }
    RandstrobeGenerator randstrobe_iter{seq, parameters.syncmer, parameters.randstrobe};
    std::vector<Randstrobe> chunk;
    // TODO
    // Chunking makes this function faster, but the speedup is achieved even
    // with a chunk size of 1.
    const size_t chunk_size = 4;
    chunk.reserve(chunk_size);
    bool end = false;
    while (!end) {
        // fill chunk
        Randstrobe randstrobe;
        while (chunk.size() < chunk_size) {
            randstrobe = randstrobe_iter.next();
            if (randstrobe == randstrobe_iter.end()) {
                end = true;
                break;
            }
            chunk.push_back(randstrobe);
        }
        for (auto randstrobe : chunk) {
            RefRandstrobe::packed_t packed = ref_index << 8;
            packed = packed + (randstrobe.strobe2_pos - randstrobe.strobe1_pos);
            randstrobes[offset++] = RefRandstrobe{randstrobe.hash, randstrobe.strobe1_pos, packed};
        }
        chunk.clear();
    }
}

void StrobemerIndex::print_diagnostics(const std::string& logfile_name, int k) const {
    // Prins to csv file the statistics on the number of seeds of a particular length and what fraction of them them are unique in the index:
    // format:
    // seed_length, count, percentage_unique

    size_t max_size = 100000;
    std::vector<int> log_count(max_size, 0);  // stores count and each index represents the length
    std::vector<int> log_unique(max_size, 0);  // stores count unique and each index represents the length
    std::vector<int> log_repetitive(max_size, 0);  // stores count unique and each index represents the length


    std::vector<randstrobe_hash_t> log_count_squared(max_size,0);
    randstrobe_hash_t tot_seed_count = 0;
    randstrobe_hash_t tot_seed_count_sq = 0;

    std::vector<randstrobe_hash_t> log_count_1000_limit(max_size, 0);  // stores count and each index represents the length
    randstrobe_hash_t tot_seed_count_1000_limit = 0;

    size_t seed_length = 0;

    for (size_t it = 0; it < randstrobes.size(); it++) {
        seed_length = strobe2_offset(it) + k;
        auto count = get_count(it);

        if (seed_length < max_size){
            log_count[seed_length] ++;
            log_count_squared[seed_length] += count;
            tot_seed_count ++;
            tot_seed_count_sq += count;
            if (count <= 1000){
                log_count_1000_limit[seed_length] ++;
                tot_seed_count_1000_limit ++;
            }
        } else {
            // TODO This function should not log anything
            // logger.info() << "Detected seed size over " << max_size << " bp (can happen, e.g., over centromere): " << seed_length << std::endl;
        }

        if (count == 1 && seed_length < max_size) {
            log_unique[seed_length]++;
        }
        if (count >= 10 && seed_length < max_size) {
            log_repetitive[seed_length]++;
        }
    }

    // printing
    std::ofstream log_file;
    log_file.open(logfile_name);

    for (size_t i = 0; i < log_count.size(); ++i) {
        if (log_count[i] > 0) {
            double e_count = log_count_squared[i] / log_count[i];
            log_file << i << ',' << log_count[i] << ',' << e_count << std::endl;
        }
    }

    // Get median
    size_t n = 0;
    int median = 0;
    for (size_t i = 0; i < log_count.size(); ++i) {
        n += log_count[i];
        if (n >= tot_seed_count/2) {
            break;
        }
    }
    // Get median 1000 limit
    size_t n_lim = 0;
    for (size_t i = 0; i < log_count_1000_limit.size(); ++i) {
        n_lim += log_count_1000_limit[i];
        if (n_lim >= tot_seed_count_1000_limit/2) {
            break;
        }
    }

    log_file << "E_size for total seeding wih max seed size m below (m, tot_seeds, E_hits)" << std::endl;
    double e_hits = (double) tot_seed_count_sq/ (double) tot_seed_count;
    double fraction_masked = 1.0 - (double) tot_seed_count_1000_limit/ (double) tot_seed_count;
    log_file << median << ',' << tot_seed_count << ',' << e_hits << ',' << 100*fraction_masked << std::endl;
}
*/
