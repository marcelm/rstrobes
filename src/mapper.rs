use std::cmp::{min, Reverse};
use std::mem;
use fastrand::Rng;
use crate::aligner::{AlignmentInfo, hamming_align, hamming_distance};
use crate::cigar::{Cigar, CigarOperation};
use crate::fasta::RefSequence;
use crate::index::{IndexParameters, StrobemerIndex};
use crate::nam::{find_nams, find_nams_rescue, Nam, reverse_nam_if_needed};
use crate::revcomp::reverse_complement;
use crate::strobes::RandstrobeIterator;
use crate::syncmers::SyncmerIterator;
use crate::sam::{SamRecord, REVERSE, SECONDARY, UNMAP};
use crate::read::Read;
use crate::aligner::Aligner;
use crate::details::Details;
use crate::fastq::SequenceRecord;
use crate::insertsize::InsertSizeDistribution;


pub struct MappingParameters {
    r: usize,
    max_secondary: usize,
    dropoff_threshold: f32,
    rescue_level: usize,
    max_tries: usize,
    output_unmapped: bool,
}

impl Default for MappingParameters {
    fn default() -> Self {
        MappingParameters {
            r: 150,
            max_secondary: 0,
            dropoff_threshold: 0.5,
            rescue_level: 2,
            max_tries: 20,
            output_unmapped: true,
        }
    }
}

#[derive(Clone,Default)]
struct Alignment {
    reference_id: usize,
    ref_start: usize,
    cigar: Cigar,
    edit_distance: usize,
    soft_clip_left: usize,
    soft_clip_right: usize,
    score: u32,
    length: usize,
    is_revcomp: bool,
    is_unaligned: bool, // TODO get rid of this
    /// Whether a gapped alignment function was used to obtain this alignment
    /// (even if true, the alignment can still be without gaps)
    gapped: bool,
}

impl Alignment {
    fn global_edit_distance(&self) -> usize {
        self.edit_distance + self.soft_clip_left + self.soft_clip_right
    }
}

#[derive(Debug)]
pub struct QueryRandstrobe {
    pub hash: u64,
    pub start: usize,
    pub end: usize,
    pub is_revcomp: bool,
}

    /// Generate randstrobes for a query sequence and its reverse complement.
pub fn randstrobes_query(seq: &[u8], parameters: &IndexParameters) -> Vec<QueryRandstrobe> {
    let mut randstrobes= Vec::<QueryRandstrobe>::new();
    if seq.len() < parameters.randstrobe.w_max {
        return randstrobes;
    }

    // TODO
    // For the reverse complement, we could re-use the syncmers of the forward
    // sequence because canonical syncmers are invariant under reverse
    // complementing. Only the coordinates need to be adjusted.

    let seq_rc = reverse_complement(seq);
    for (s, is_reverse) in [(seq, false), (&seq_rc, true)] {
        // Generate randstrobes for the forward sequence
        let mut syncmer_iter = SyncmerIterator::new(s, parameters.syncmer.k, parameters.syncmer.s, parameters.syncmer.t);
        let randstrobe_iter = RandstrobeIterator::new(&mut syncmer_iter, &parameters.randstrobe);

        for randstrobe in randstrobe_iter {
            randstrobes.push(
                QueryRandstrobe {
                    hash: randstrobe.hash,
                    start: randstrobe.strobe1_pos,
                    end: randstrobe.strobe2_pos + parameters.syncmer.k,
                    is_revcomp: is_reverse
                }
            );
        }
    }

    randstrobes
}

/// Conversion of an Alignment into a SamRecord
#[derive(Default)]
pub struct SamOutput {
    cigar_eqx: bool,
    details: bool,
    rg_id: Option<String>,
}

impl SamOutput {
    pub fn new(details: bool, cigar_eqx: bool, rg_id: Option<String>) -> Self {
        SamOutput {
            cigar_eqx,
            details,
            rg_id,
        }
    }

    /// Convert the alignment into a SamRecord
    fn make_record(&self, alignment: &Alignment, references: &[RefSequence], record: &SequenceRecord, mut mapq: u8, is_primary: bool, details: Details) -> SamRecord {
        let mut flags = 0;

        if !alignment.is_unaligned && alignment.is_revcomp {
            flags |= REVERSE;
        }
        if !is_primary {
            mapq = 255;
            flags |= SECONDARY;
        }

        let query_sequence = if alignment.is_revcomp {
            reverse_complement(&record.sequence)
        } else {
            record.sequence.clone()
        };
        let query_qualities = if alignment.is_revcomp {
            let mut rev = record.qualities.clone();
            rev.reverse();
            rev
        } else {
            record.qualities.clone()
        };
        let mut cigar = Cigar::new();
        cigar.push(CigarOperation::Softclip, alignment.soft_clip_left);
        cigar.extend(&alignment.cigar);
        cigar.push(CigarOperation::Softclip, alignment.soft_clip_right);

        let reference_name = Some(references[alignment.reference_id].name.clone());
        let details = if self.details { Some(details) } else { None };

        let cigar = if self.cigar_eqx { Some(cigar) } else { Some(cigar.with_m()) };
        SamRecord {
            query_name: record.name.clone(),
            flags,
            reference_name,
            pos: Some(alignment.ref_start as u32),
            mapq,
            cigar,
            query_sequence: Some(query_sequence),
            query_qualities: Some(query_qualities),
            edit_distance: Some(alignment.edit_distance as u32),
            alignment_score: Some(alignment.score),
            details,
            rg_id: self.rg_id.clone(),
            ..SamRecord::default()
        }
    }

    fn make_unmapped_record(&self, record: &SequenceRecord, details: Details) -> SamRecord {
        let details = if self.details { Some(details) } else { None };
        SamRecord {
            query_name: record.name.clone(),
            flags: UNMAP,
            query_sequence: Some(record.sequence.clone()),
            query_qualities: Some(record.qualities.clone()),
            details,
            rg_id: self.rg_id.clone(),
            ..SamRecord::default()
        }
    }
}

pub fn map_single_end_read(
    record: &SequenceRecord,
    index: &StrobemerIndex,
    references: &[RefSequence],
    mapping_parameters: &MappingParameters,
    sam_output: &SamOutput,
    aligner: &Aligner,
) -> Vec<SamRecord> {
    let mut details = Details::default();

    let (was_rescued, mut nams) = get_nams(&record.sequence, &index, &mapping_parameters);
    details.nam_rescue = was_rescued;
    details.nams = nams.len();

    // Timer extend_timer;
    // align_SE(
    //     aligner, sam, nams, record, index.parameters.syncmer.k,
    //     references, details, map_param.dropoff_threshold, map_param.maxTries,
    //     map_param.max_secondary
    // );

    if nams.is_empty() {
        return vec![sam_output.make_unmapped_record(record, details)];
    }
    let mut sam_records = Vec::new();
    let mut alignments = Vec::new();
    let nam_max = nams[0];
    let mut best_edit_distance = usize::MAX;
    let mut best_score = 0;
    let mut second_best_score = 0;
    let mut best_alignment = None;

    let k = index.parameters.syncmer.k;
    let read = Read::new(&record.sequence);

    for (tries, nam) in nams.iter_mut().enumerate() {
        let score_dropoff = nam.n_hits as f32 / nam_max.n_hits as f32;

        // TODO iterate over slice of nams instead of tracking tries
        if tries >= mapping_parameters.max_tries || (tries > 1 && best_edit_distance == 0) || score_dropoff < mapping_parameters.dropoff_threshold {
            break;
        }
        let consistent_nam = reverse_nam_if_needed(nam, &read, references, k);
        details.nam_inconsistent += (!consistent_nam) as usize;
        let alignment = extend_seed(aligner, nam, references, &read, consistent_nam);
        if alignment.is_none() {
            continue;
        }
        let alignment = alignment.unwrap();
        details.tried_alignment += 1;
        details.gapped += alignment.gapped as usize;

        if alignment.score > best_score {
            second_best_score = best_score;
            best_score = alignment.score;
            best_alignment = Some(alignment.clone());
            if mapping_parameters.max_secondary == 0 {
                best_edit_distance = alignment.global_edit_distance();
            }
        } else if alignment.score > second_best_score {
            second_best_score = alignment.score;
        }
        if mapping_parameters.max_secondary > 0 {
            alignments.push(alignment);
        }
    }
    let mapq = ((60 * (best_score - second_best_score) + best_score - 1) / best_score) as u8;

    if best_alignment.is_none() {
        return vec![sam_output.make_unmapped_record(record, details)];
    }
    let best_alignment = best_alignment.unwrap();
    if mapping_parameters.max_secondary == 0 {
        sam_records.push(
            sam_output.make_record(&best_alignment, references, record, mapq, true, details)
        );
    } else {
        // Highest score first
        alignments.sort_by_key(|k| Reverse(k.score));

        let max_out = min(alignments.len(), mapping_parameters.max_secondary + 1);
        for (i, alignment) in alignments.iter().enumerate().take(max_out) {
            let is_primary = i == 0;
            if alignment.score - best_score > 2 * aligner.scores.mismatch as u32 + aligner.scores.gap_open as u32 {
                break;
            }
            // TODO .clone()
            sam_records.push(sam_output.make_record(alignment, references, record, mapq, is_primary, details.clone()));
        }
    }
    // statistics.tot_extend += extend_timer.duration();
    // statistics += details;
    sam_records
}

// TODO rename
fn get_nams(sequence: &[u8], index: &StrobemerIndex, mapping_parameters: &MappingParameters) -> (bool, Vec<Nam>) {
    //Timer strobe_timer;
    let query_randstrobes = randstrobes_query(&sequence, &index.parameters);
    //statistics.tot_construct_strobemers += strobe_timer.duration();

    // Timer nam_timer;
    let (nonrepetitive_fraction, mut nams) = find_nams(&query_randstrobes, index, index.filter_cutoff);
    // statistics.tot_find_nams += nam_timer.duration();

    let mut nam_rescue = false;
    if mapping_parameters.rescue_level > 1 {
        // Timer rescue_timer;
        if nams.is_empty() || nonrepetitive_fraction < 0.7 {
            nam_rescue = true;
            nams = find_nams_rescue(&query_randstrobes, index, index.rescue_cutoff);
        }
        // statistics.tot_time_rescue += rescue_timer.duration();
    }
    // Timer nam_sort_timer;

    nams.sort_by_key(|&k| -(k.score as i32));
    // TODO shuffle_top_nams(nams, random_engine);
    // statistics.tot_sort_nams += nam_sort_timer.duration();

    (nam_rescue, nams)
}

/// Extend a NAM so that it covers the entire read and return the resulting
/// alignment.
fn extend_seed(
    aligner: &Aligner,
    nam: &Nam,
    references: &[RefSequence],
    read: &Read,
    consistent_nam: bool,
) -> Option<Alignment> {
    let query = if nam.is_revcomp { read.rc() } else { read.seq() };
    let refseq = &references[nam.ref_id].sequence;

    let projected_ref_start = nam.ref_start.saturating_sub(nam.query_start);
    let projected_ref_end = min(nam.ref_end + query.len() - nam.query_end, refseq.len());

    // TODO ugly
    let mut info = AlignmentInfo {
        cigar: Default::default(),
        edit_distance: 0,
        ref_start: 0,
        ref_end: 0,
        query_start: 0,
        query_end: 0,
        score: 0,
    };
    let mut result_ref_start = 0;
    let mut gapped = true;
    if projected_ref_start + query.len() == projected_ref_end && consistent_nam {
        let ref_segm_ham = &refseq[projected_ref_start..projected_ref_end];
        if let Some(hamming_dist) = hamming_distance(query, ref_segm_ham) {
            if (hamming_dist as f32 / query.len() as f32) < 0.05 {
                // ungapped worked fine, no need to do gapped alignment
                info = hamming_align(query, ref_segm_ham, aligner.scores.match_, aligner.scores.mismatch, aligner.scores.end_bonus).expect(
                    "hamming_dist was successful, this should be as well"
                );
                result_ref_start = projected_ref_start + info.ref_start;
                gapped = false;
            }
        }
    }
    if gapped {
        let ref_start = projected_ref_start.saturating_sub(50);
        let ref_end = min(nam.ref_end + 50, refseq.len());
        let segment = &refseq[ref_start..ref_end];
        info = aligner.align(query, segment)?;
        result_ref_start = ref_start + info.ref_start;
    }
    Some(Alignment {
        cigar: info.cigar.clone(),
        edit_distance: info.edit_distance,
        soft_clip_left: info.query_start,
        soft_clip_right: query.len() - info.query_end,
        score: info.score,
        ref_start: result_ref_start,
        length: info.ref_span(),
        is_revcomp: nam.is_revcomp,
        is_unaligned: false,
        reference_id: nam.ref_id,
        gapped,
    })
}

// TODO alignment statistics
pub fn map_paired_end_read(
    r1: &SequenceRecord,
    r2: &SequenceRecord,
    index: &StrobemerIndex,
    references: &[RefSequence],
    mapping_parameters: &MappingParameters,
    index_parameters: &IndexParameters,
    insert_size_distribution: &mut InsertSizeDistribution,
    aligner: &Aligner,
    rng: &mut Rng,
) -> Vec<SamRecord> {
    let mut details = [Details::default(), Details::default()];
    let mut nams_pair = vec![vec![]];

    for is_revcomp in [0, 1] {
        let record = if is_revcomp == 0 { r1 } else { r2 };
        let (was_rescued, mut nams) = get_nams(&record.sequence, &index, &mapping_parameters);
        details[is_revcomp].nam_rescue = was_rescued;
        details[is_revcomp].nams = nams.len();
        nams_pair.push(nams);
    }

    // Timer extend_timer;
    let read1 = Read::new(&r1.sequence);
    let read2 = Read::new(&r2.sequence);
    let mut alignment_pairs = align_paired(
        aligner, &nams_pair[0], &nams_pair[1], &read1, &read2,
        index_parameters.syncmer.k, references, &mut details,
        mapping_parameters.dropoff_threshold, insert_size_distribution,
        mapping_parameters.max_tries
    );

    let mut sam_records = Vec::new();

    // -1 marks the typical case that both reads map uniquely and form a
    // proper pair. Then the mapping quality is computed based on the NAMs.
    if alignment_pairs.len() == 1 && alignment_pairs[0].score == -1 {
        let alignment1 = &alignment_pairs[0].alignment1;
        let alignment2 = &alignment_pairs[0].alignment2;
        let is_proper = is_proper_pair(alignment1, alignment2, insert_size_distribution.mu, insert_size_distribution.sigma);
        if is_proper
            && insert_size_distribution.sample_size < 400
            && alignment1.edit_distance + alignment2.edit_distance < 3
        {
            insert_size_distribution.update(alignment1.ref_start.abs_diff(alignment2.ref_start));
        }

        let mapq1 = proper_pair_mapq(&nams_pair[0]);
        let mapq2 = proper_pair_mapq(&nams_pair[1]);

        details[0].best_alignments = 1;
        details[1].best_alignments = 1;
        let is_primary = true;
        sam_records.extend(alignment_pairs[0].as_sam_records(&r1, &r2, read1.rc(), read2.rc(), mapq1, mapq2, is_proper, is_primary, details));
    } else {
        alignment_pairs.sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap());
        deduplicate_scored_pairs(&mut alignment_pairs);

        // If there are multiple top-scoring alignments (all with the same score),
        // pick one randomly and move it to the front.
        let i = count_best_alignment_pairs(&alignment_pairs);
        details[0].best_alignments = i;
        details[1].best_alignments = i;
        if i > 1 {
            let random_index = rng.usize(..i);
            if random_index != 0 {
                mem::swap(&mut alignment_pairs[0], &mut alignment_pairs[random_index]);
            }
        }

        let secondary_dropoff = 2 * aligner.scores.mismatch + aligner.scores.gap_open;
        output_aligned_pairs(
            alignment_pairs,
            sam,
            map_param.max_secondary,
            secondary_dropoff,
            &r1,
            &r2,
            &read1,
            &read2,
            insert_size_distribution.mu,
            insert_size_distribution.sigma,
            details
        );
    }
    // TODO
    // statistics.tot_extend += extend_timer.duration();
    // statistics += details[0];
    // statistics += details[1];

    sam_records
}

fn align_paired(
    aligner: &Aligner,
    nams1: &[Nam],
    nams2: &[Nam],
    read1: &Read,
    read2: &Read,
    k: usize,
    references: &[RefSequence],
    details: &mut [Details; 2],
    dropoff: f32,
    insert_size_distribution: &InsertSizeDistribution,
    max_tries: usize,
) -> Vec<ScoredAlignmentPair> {
    let mu = insert_size_distribution.mu;
    let sigma = insert_size_distribution.sigma;

    if nams1.is_empty() && nams2.is_empty() {
         // None of the reads have any NAMs
        return vec![];
    }

    if !nams1.is_empty() && nams2.is_empty() {
        // Only read 1 has NAMS: attempt to rescue read 2
        return rescue_read(
            read2,
            read1,
            aligner,
            references,
            nams1,
            max_tries,
            dropoff,
            details,
            k,
            mu,
            sigma
        );
    }

    if nams1.is_empty() && !nams2.is_empty() {
        // Only read 2 has NAMS: attempt to rescue read 1
        let mut swapped_details = [details[0].clone(), details[1].clone()];
        let pairs = rescue_read(
            read1,
            read2,
            aligner,
            references,
            nams2,
            max_tries,
            dropoff,
            swapped_details,
            k,
            mu,
            sigma
        );
        details[0] += swapped_details[1];
        details[1] += swapped_details[0];
        for mut pair in &pairs {
            mem::swap(&mut pair.alignment1, &mut pair.alignment2);
        }

        return pairs;
    }

    // Both reads have NAMs
    assert!(!nams1.is_empty() && !nams2.is_empty());

    // Deal with the typical case that both reads map uniquely and form a proper pair
    if top_dropoff(&nams1) < dropoff && top_dropoff(&nams2) < dropoff && is_proper_nam_pair(&nams1[0], &nams2[0], mu, sigma) {
        let mut n_max1 = nams1[0];
        let mut n_max2 = nams2[0];

        let consistent_nam1 = reverse_nam_if_needed(&mut n_max1, read1, references, k);
        details[0].nam_inconsistent += (!consistent_nam1).into();
        let consistent_nam2 = reverse_nam_if_needed(&mut n_max2, read2, references, k);
        details[1].nam_inconsistent += (!consistent_nam2).into();

        let alignment1 = extend_seed(aligner, &mut n_max1, references, read1, consistent_nam1);
        let alignment2 = extend_seed(aligner, &mut n_max2, references, read2, consistent_nam2);
        if let (Some(alignment1), Some(alignment2)) = (alignment1, alignment2) {
            details[0].tried_alignment += 1;
            details[0].gapped += alignment1.gapped as usize;
            details[1].tried_alignment += 1;
            details[1].gapped += alignment2.gapped as usize;

            return vec![ScoredAlignmentPair{score: -1, alignment1, alignment2}];
        }
        // TODO what if one of the alignments is None?
    }

    // Do a full search for highest-scoring pair
    // Get top hit counts for all locations.
    // The joint hit count is the sum of hits of the two mates.
    // Then align as long as score dropoff or cnt < 20

    let nam_pairs = get_best_scoring_nam_pairs(&nams1, &nams2, mu, sigma);

    // Cache for already computed alignments. Maps NAM ids to alignments.
    robin_hood::unordered_map<int,Alignment> is_aligned1;
    robin_hood::unordered_map<int,Alignment> is_aligned2;

    // These keep track of the alignments that would be best if we treated
    // the paired-end read as two single-end reads.
    Alignment a1_indv_max, a2_indv_max;
    {
        auto n1_max = nams1[0];
        bool consistent_nam1 = reverse_nam_if_needed(n1_max, read1, references, k);
        details[0].nam_inconsistent += !consistent_nam1;
        a1_indv_max = extend_seed(aligner, n1_max, references, read1, consistent_nam1);
        is_aligned1[n1_max.nam_id] = a1_indv_max;
        details[0].tried_alignment++;
        details[0].gapped += a1_indv_max.gapped;

        auto n2_max = nams2[0];
        bool consistent_nam2 = reverse_nam_if_needed(n2_max, read2, references, k);
        details[1].nam_inconsistent += !consistent_nam2;
        a2_indv_max = extend_seed(aligner, n2_max, references, read2, consistent_nam2);
        is_aligned2[n2_max.nam_id] = a2_indv_max;
        details[1].tried_alignment++;
        details[1].gapped += a2_indv_max.gapped;
    }

    // Turn pairs of high-scoring NAMs into pairs of alignments
    std::vector<ScoredAlignmentPair> high_scores;
    auto max_score = nam_pairs[0].n_hits;
    for (auto &[score_, n1, n2] : nam_pairs) {
        float score_dropoff = (float) score_ / max_score;

        if (high_scores.size() >= max_tries || score_dropoff < dropoff) {
            break;
        }

        // Get alignments for the two NAMs, either by computing the alignment,
        // retrieving it from the cache or by attempting a rescue (if the NAM
        // actually is a dummy, that is, only the partner is available)
        Alignment a1;
        // ref_start == -1 is a marker for a dummy NAM
        if (n1.ref_start >= 0) {
            if (is_aligned1.find(n1.nam_id) != is_aligned1.end() ){
                a1 = is_aligned1[n1.nam_id];
            } else {
                bool consistent_nam = reverse_nam_if_needed(n1, read1, references, k);
                details[0].nam_inconsistent += !consistent_nam;
                a1 = extend_seed(aligner, n1, references, read1, consistent_nam);
                is_aligned1[n1.nam_id] = a1;
                details[0].tried_alignment++;
                details[0].gapped += a1.gapped;
            }
        } else {
            details[1].nam_inconsistent += !reverse_nam_if_needed(n2, read2, references, k);
            a1 = rescue_align(aligner, n2, references, read1, mu, sigma, k);
            details[0].mate_rescue += !a1.is_unaligned;
            details[0].tried_alignment++;
        }
        if (a1.score > a1_indv_max.score) {
            a1_indv_max = a1;
        }

        Alignment a2;
        // ref_start == -1 is a marker for a dummy NAM
        if (n2.ref_start >= 0) {
            if (is_aligned2.find(n2.nam_id) != is_aligned2.end() ){
                a2 = is_aligned2[n2.nam_id];
            } else {
                bool consistent_nam = reverse_nam_if_needed(n2, read2, references, k);
                details[1].nam_inconsistent += !consistent_nam;
                a2 = extend_seed(aligner, n2, references, read2, consistent_nam);
                is_aligned2[n2.nam_id] = a2;
                details[1].tried_alignment++;
                details[1].gapped += a2.gapped;
            }
        } else {
            details[0].nam_inconsistent += !reverse_nam_if_needed(n1, read1, references, k);
            a2 = rescue_align(aligner, n1, references, read2, mu, sigma, k);
            details[1].mate_rescue += !a2.is_unaligned;
            details[1].tried_alignment++;
        }
        if (a2.score > a2_indv_max.score){
            a2_indv_max = a2;
        }

        bool r1_r2 = a2.is_rc && (a1.ref_start <= a2.ref_start) && ((a2.ref_start - a1.ref_start) < mu + 10*sigma); // r1 ---> <---- r2
        bool r2_r1 = a1.is_rc && (a2.ref_start <= a1.ref_start) && ((a1.ref_start - a2.ref_start) < mu + 10*sigma); // r2 ---> <---- r1

        double combined_score;
        if (r1_r2 || r2_r1) {
            // Treat a1/a2 as a pair
            float x = std::abs(a1.ref_start - a2.ref_start);
            combined_score = (double)a1.score + (double)a2.score + std::max(-20.0f + 0.001f, log(normal_pdf(x, mu, sigma)));
            //* (1 - s2 / s1) * min_matches * log(s1);
        } else {
            // Treat a1/a2 as two single-end reads
            // 20 corresponds to a value of log(normal_pdf(x, mu, sigma)) of more than 5 stddevs away (for most reasonable values of stddev)
            combined_score = (double)a1.score + (double)a2.score - 20;
        }

        ScoredAlignmentPair aln_pair{combined_score, a1, a2};
        high_scores.push_back(aln_pair);
    }

    // Finally, add highest scores of both mates as individually mapped
    double combined_score = (double)a1_indv_max.score + (double)a2_indv_max.score - 20; // 20 corresponds to  a value of log( normal_pdf(x, mu, sigma ) ) of more than 5 stddevs away (for most reasonable values of stddev)
    ScoredAlignmentPair aln_tuple{combined_score, a1_indv_max, a2_indv_max};
    high_scores.push_back(aln_tuple);

    return high_scores;
}

fn is_proper_pair(alignment1: &Alignment, alignment2: &Alignment, mu: f32, sigma: f32) -> bool {
    let dist = alignment2.ref_start as isize - alignment1.ref_start as isize;
    let same_reference = alignment1.reference_id == alignment2.reference_id;
    let both_aligned = same_reference && !alignment1.is_unaligned && !alignment2.is_unaligned;
    let r1_r2 = !alignment1.is_revcomp && alignment2.is_revcomp && dist >= 0; // r1 ---> <---- r2
    let r2_r1 = !alignment2.is_revcomp && alignment1.is_revcomp && dist <= 0; // r2 ---> <---- r1
    let rel_orientation_good = r1_r2 || r2_r1;
    let insert_good = dist.unsigned_abs() <= (mu + sigma * 6.0) as usize;

    both_aligned && insert_good && rel_orientation_good
}

fn is_proper_nam_pair(nam1: &Nam, nam2: &Nam, mu: f32, sigma: f32) -> bool {
    if nam1.ref_id != nam2.ref_id || nam1.is_revcomp == nam2.is_revcomp {
        return false;
    }
    let r1_ref_start = nam1.projected_ref_start();
    let r2_ref_start = nam2.projected_ref_start();

    // r1 ---> <---- r2
    let r1_r2 = nam2.is_revcomp && (r1_ref_start <= r2_ref_start) && ((r2_ref_start - r1_ref_start) as f32) < mu + 10.0 * sigma;

     // r2 ---> <---- r1
    let r2_r1 = nam1.is_revcomp && (r2_ref_start <= r1_ref_start) && ((r1_ref_start - r2_ref_start) as f32) < mu + 10.0 * sigma;

    return r1_r2 || r2_r1;
}

fn proper_pair_mapq(nams: &[Nam]) -> u8 {
    if nams.len() <= 1 {
        return 60;
    }
    let s1 = nams[0].score;
    let s2 = nams[1].score;
    // from minimap2: MAPQ = 40(1−s2/s1) ·min{1,|M|/10} · log s1
    let min_matches = min(nams[0].n_hits, 10) as f32 / 10.0;
    let uncapped_mapq = 40.0 * (1 - s2 / s1) as f32 * min_matches * (s1 as f32).ln();

    uncapped_mapq.min(60.0) as u8
}

pub struct ScoredAlignmentPair {
    score: f64,
    alignment1: Alignment,
    alignment2: Alignment,
}

fn count_best_alignment_pairs(pairs: &[ScoredAlignmentPair]) -> usize {
    if pairs.is_empty() {
        0
    } else {
        pairs.iter().take_while(|x| x.score == pairs[0].score).count()
    }
}

/// compute dropoff of the first (top) NAM
fn top_dropoff(nams: &[Nam]) -> f32 {
    let n_max = &nams[0];
    if n_max.n_hits <= 2 {
        1.0
    } else if nams.len() > 1 {
        nams[1].n_hits as f32 / n_max.n_hits as f32
    } else {
        0.0
    }
}

#[cfg(test)]
mod tests {
    use crate::mapper::{Alignment, count_best_alignment_pairs, ScoredAlignmentPair};

    #[test]
    fn test_count_best_alignment_pairs() {
        let mut pairs = vec![];
        fn add_alignment(pairs: &mut Vec<ScoredAlignmentPair>, score: f64) {
            pairs.push(ScoredAlignmentPair { score, alignment1: Alignment::default(), alignment2: Alignment::default() });
        }

        assert_eq!(count_best_alignment_pairs(&pairs), 0);
        add_alignment(&mut pairs, 10.0);
        assert_eq!(count_best_alignment_pairs(&pairs), 1);

        add_alignment(&mut pairs, 10.0);
        assert_eq!(count_best_alignment_pairs(&pairs), 2);

        add_alignment(&mut pairs, 5.0);
        assert_eq!(count_best_alignment_pairs(&pairs), 2);

        pairs[1].score = 5.0;
        assert_eq!(count_best_alignment_pairs(&pairs), 1);
    }
}
