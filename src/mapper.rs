use std::cmp::{min, Reverse};
use std::collections::{HashMap, HashSet};
use std::mem;
use std::ops::Index;
use fastrand::Rng;
use memchr::memmem;
use crate::aligner::{AlignmentInfo, hamming_align, hamming_distance};
use crate::cigar::{Cigar, CigarOperation};
use crate::fasta::RefSequence;
use crate::index::{IndexParameters, StrobemerIndex};
use crate::nam::{find_nams, find_nams_rescue, Nam, reverse_nam_if_needed};
use crate::revcomp::reverse_complement;
use crate::strobes::RandstrobeIterator;
use crate::syncmers::SyncmerIterator;
use crate::sam::{SamRecord, REVERSE, SECONDARY, UNMAP, PAIRED, READ1, READ2, PROPER_PAIR, MUNMAP, MREVERSE};
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
    fn make_record(
        &self, alignment: &Alignment, references: &[RefSequence], record: &SequenceRecord, mut mapq: u8, is_primary: bool, details: Details
    ) -> SamRecord {
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

    // alignment: &Alignment, references: &[RefSequence], record: &SequenceRecord, mut mapq: u8, is_primary: bool, details: Details) -> SamRecord {
    fn make_paired_records(
        &self,
        alignments: [&Alignment; 2],
        references: &[RefSequence],
        records: [&SequenceRecord; 2],
        mapq: [u8; 2],
        details: &[Details; 2],
        is_primary: bool,
        is_proper: bool,
    ) -> [SamRecord; 2] {
        let mut sam_records = [
            self.make_record(&alignments[0], references, &records[0], mapq[0], is_primary, details[0].clone()),
            self.make_record(&alignments[1], references, &records[1], mapq[1], is_primary, details[1].clone()),
        ];
        sam_records[0].flags |= PAIRED | READ1;
        sam_records[1].flags |= PAIRED | READ2;

        let mut template_len1 = 0;
        let both_aligned = !alignments[0].is_unaligned && !alignments[1].is_unaligned;
        if both_aligned && alignments[0].reference_id == alignments[1].reference_id {
            let dist = alignments[1].ref_start - alignments[0].ref_start;
            if alignments[1].ref_start > alignments[0].ref_start {
                template_len1 = (alignments[1].ref_start - alignments[0].ref_start + alignments[1].length) as isize;
            }
            else {
                template_len1 = -((alignments[0].ref_start - alignments[1].ref_start) as isize) - alignments[0].length as isize;
            }
        }
        if is_proper {
            sam_records[0].flags |= PROPER_PAIR;
            sam_records[1].flags |= PROPER_PAIR;
        }

        let reference_name1;
        let mut pos1 = Some(alignments[0].ref_start);
        let edit_distance1 = alignments[0].edit_distance;
        if alignments[0].is_unaligned {
            sam_records[0].flags |= UNMAP;
            sam_records[1].flags |= MUNMAP;
            pos1 = None;
            reference_name1 = "*";
        } else {
            if alignments[0].is_revcomp {
                sam_records[0].flags |= REVERSE;
                sam_records[1].flags |= MREVERSE;
            }
            reference_name1 = &references[alignments[0].reference_id].name;
        }

        let reference_name2;
        let mut pos2 = Some(alignments[1].ref_start);
        if alignments[1].is_unaligned {
            sam_records[1].flags |= UNMAP;
            sam_records[0].flags |= MUNMAP;
            pos2 = None;
            reference_name2 = "*";
        } else {
            if alignments[1].is_revcomp {
                sam_records[0].flags |= MREVERSE;
                sam_records[1].flags |= REVERSE;
            }
            reference_name2 = &references[alignments[1].reference_id].name;
        }

        // Reference name as used in the RNEXT field;
        // set to "=" if identical to reference_name
        let mut mate_reference_name1 = reference_name1;
        let mut mate_reference_name2 = reference_name2;
        if
            (!alignments[0].is_unaligned && !alignments[1].is_unaligned && alignments[0].reference_id == alignments[1].reference_id)
            || (alignments[0].is_unaligned != alignments[1].is_unaligned)
        {
            mate_reference_name1 = "=";
            mate_reference_name2 = "=";
        }

        if alignments[0].is_unaligned != alignments[1].is_unaligned {
            if alignments[0].is_unaligned {
                pos1 = pos2;
            } else {
                pos2 = pos1;
            }
        }

        /* TODO
        if alignments[0].is_unaligned {
            add_unmapped_mate(record1, f1, reference_name2, pos2);
        } else {
            add_record(record1.name, record1.comment, f1, reference_name1, alignment1.ref_start, mapq1, alignment1.cigar, mate_reference_name2, pos2, template_len1, record1.seq, read1_rc, record1.qual, edit_distance1, alignment1.score, details[0]);
        }
        if alignments[1].is_unaligned {
            add_unmapped_mate(record2, f2, reference_name1, pos1);
        } else {
            add_record(record2.name, record2.comment, f2, reference_name2, alignment2.ref_start, mapq2, alignment2.cigar, mate_reference_name1, pos1, -template_len1, record2.seq, read2_rc, record2.qual, edit_distance2, alignment2.score, details[1]);
        }
        */
        sam_records
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
    sam_output: &SamOutput,
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
        aligner, &mut nams_pair[0], &mut nams_pair[1], &read1, &read2,
        index_parameters.syncmer.k, references, &mut details,
        mapping_parameters.dropoff_threshold, insert_size_distribution,
        mapping_parameters.max_tries
    );

    let mut sam_records = Vec::new();

    // score set to None marks the
    match alignment_pairs {
        // Typical case: both reads map uniquely and form a proper pair.
        // Then the mapping quality is computed based on the NAMs.
        AlignedPairs::Proper((alignment1, alignment2)) => {
            let is_proper = is_proper_pair(&alignment1, &alignment2, insert_size_distribution.mu, insert_size_distribution.sigma);
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

            sam_records.extend(
                sam_output.make_paired_records([&alignment1, &alignment2], references, [r1, r2], [mapq1, mapq2], &details, is_primary, is_proper)
            );
        },
        AlignedPairs::WithScores(mut alignment_pairs) => {
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
            sam_records.extend(aligned_pairs_to_sam(
                &alignment_pairs,
                mapping_parameters.max_secondary,
                secondary_dropoff as f64,
                &r1,
                &r2,
                &read1,
                &read2,
                insert_size_distribution.mu,
                insert_size_distribution.sigma,
                &details
            ));
        }
    }
    // TODO
    // statistics.tot_extend += extend_timer.duration();
    // statistics += details[0];
    // statistics += details[1];

    sam_records
}

enum AlignedPairs {
    Proper((Alignment, Alignment)),
    WithScores(Vec<ScoredAlignmentPair>),
}

fn align_paired(
    aligner: &Aligner,
    nams1: &mut [Nam],
    nams2: &mut [Nam],
    read1: &Read,
    read2: &Read,
    k: usize,
    references: &[RefSequence],
    details: &mut [Details; 2],
    dropoff: f32,
    insert_size_distribution: &InsertSizeDistribution,
    max_tries: usize,
) -> AlignedPairs {
    let mu = insert_size_distribution.mu;
    let sigma = insert_size_distribution.sigma;

    if nams1.is_empty() && nams2.is_empty() {
         // None of the reads have any NAMs
        return AlignedPairs::WithScores(vec![]);
    }

    if !nams1.is_empty() && nams2.is_empty() {
        // Only read 1 has NAMS: attempt to rescue read 2
        return AlignedPairs::WithScores(rescue_read(
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
        ));
    }

    if nams1.is_empty() && !nams2.is_empty() {
        // Only read 2 has NAMS: attempt to rescue read 1
        let mut swapped_details = [details[0].clone(), details[1].clone()];
        let mut pairs = rescue_read(
            read1,
            read2,
            aligner,
            references,
            nams2,
            max_tries,
            dropoff,
            &mut swapped_details,
            k,
            mu,
            sigma
        );
        details[0] += swapped_details[1].clone();
        details[1] += swapped_details[0].clone();
        for pair in &mut pairs {
            mem::swap(&mut pair.alignment1, &mut pair.alignment2);
        }

        return AlignedPairs::WithScores(pairs);
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

            return AlignedPairs::Proper((alignment1, alignment2));
        }
        // TODO what if one of the alignments is None?
    }

    // Do a full search for highest-scoring pair
    // Get top hit counts for all locations.
    // The joint hit count is the sum of hits of the two mates.
    // Then align as long as score dropoff or cnt < 20

    let nam_pairs = get_best_scoring_nam_pairs(&nams1, &nams2, mu, sigma);

    // Cache for already computed alignments. Maps NAM ids to alignments.
    let mut is_aligned1 = HashMap::new();
    let mut is_aligned2 = HashMap::new();

    // These keep track of the alignments that would be best if we treated
    // the paired-end read as two single-end reads.
    let a1_indv_max;
    let a2_indv_max;
    {
        let mut n1_max = nams1[0];
        let consistent_nam1 = reverse_nam_if_needed(&mut n1_max, read1, references, k);
        details[0].nam_inconsistent += !(consistent_nam1 as usize);
        // TODO unwrap
        a1_indv_max = extend_seed(aligner, &n1_max, references, read1, consistent_nam1).unwrap();
        details[0].tried_alignment += 1;
        details[0].gapped += a1_indv_max.gapped as usize;
        is_aligned1.insert(n1_max.nam_id, a1_indv_max);

        let mut n2_max = nams2[0];
        let consistent_nam2 = reverse_nam_if_needed(&mut n2_max, read2, references, k);
        details[1].nam_inconsistent += !(consistent_nam2 as usize);
        a2_indv_max = extend_seed(aligner, &mut n2_max, references, read2, consistent_nam2).unwrap();
        details[1].tried_alignment += 1;
        details[1].gapped += a2_indv_max.gapped as usize;
        is_aligned2.insert(n2_max.nam_id, a2_indv_max);
    }

    // Turn pairs of high-scoring NAMs into pairs of alignments
    let mut high_scores = vec![];
    let max_score = nam_pairs[0].n_hits;
    for nam_pair in nam_pairs {
        let score_ = nam_pair.n_hits;
        let mut n1 = nam_pair.nam1;
        let mut n2 = nam_pair.nam2;
        let score_dropoff = score_ as f32 / max_score as f32;

        if high_scores.len() >= max_tries || score_dropoff < dropoff {
            break;
        }

        // Get alignments for the two NAMs, either by computing the alignment,
        // retrieving it from the cache or by attempting a rescue (if the NAM
        // actually is a dummy, that is, only the partner is available)
        let mut a1;
        // TODO get rid of dummy NAMs
        // ref_start == -1 is a marker for a dummy NAM
        if n1.ref_start >= 0 {
            // TODO combine .contains_key + .get, avoid unwrap
            if is_aligned1.contains_key(&n1.nam_id) {
                a1 = is_aligned1.get(&n1.nam_id).unwrap();
            } else {
                let consistent_nam = reverse_nam_if_needed(&mut n1, read1, references, k);
                details[0].nam_inconsistent += !(consistent_nam as usize);
                if let Some(a1) = extend_seed(aligner, &n1, references, read1, consistent_nam) {
                    details[0].tried_alignment += 1;
                    details[0].gapped += a1.gapped as usize;
                    is_aligned1.insert(n1.nam_id, a1);
                }
            }
        } else {
            details[1].nam_inconsistent += !(reverse_nam_if_needed(&mut n2, read2, references, k) as usize);
            a1 = rescue_align(aligner, n2, references, read1, mu, sigma, k);
            details[0].mate_rescue += !(a1.is_unaligned as usize);
            details[0].tried_alignment += 1;
        }
        if a1.score > a1_indv_max.score {
            a1_indv_max = a1;
        }

        let a2;
        // ref_start == -1 is a marker for a dummy NAM
        if n2.ref_start >= 0 {
            if is_aligned2.find(n2.nam_id) != is_aligned2.end() {
                a2 = is_aligned2[n2.nam_id];
            } else {
                let consistent_nam = reverse_nam_if_needed(&mut n2, read2, references, k);
                details[1].nam_inconsistent += !(consistent_nam as usize);
                a2 = extend_seed(aligner, n2, references, read2, consistent_nam);
                is_aligned2[n2.nam_id] = a2;
                details[1].tried_alignment++;
                details[1].gapped += a2.gapped;
            }
        } else {
            details[0].nam_inconsistent += !(reverse_nam_if_needed(&mut n1, read1, references, k) as usize);
            a2 = rescue_align(aligner, n1, references, read2, mu, sigma, k);
            details[1].mate_rescue += !(a2.is_unaligned as usize);
            details[1].tried_alignment += 1;
        }
        if (a2.score > a2_indv_max.score){
            a2_indv_max = a2;
        }

        let r1_r2 = a2.is_revcomp && (a1.ref_start <= a2.ref_start) && ((a2.ref_start - a1.ref_start) < mu + 10*sigma); // r1 ---> <---- r2
        let r2_r1 = a1.is_revcomp && (a2.ref_start <= a1.ref_start) && ((a1.ref_start - a2.ref_start) < mu + 10*sigma); // r2 ---> <---- r1

        let combined_score;
        if r1_r2 || r2_r1 {
            // Treat a1/a2 as a pair
            let x = a1.ref_start.abs_diff(a2.ref_start);
            combined_score = a1.score as f64 + a2.score as f64 + (-20.0f64 + 0.001).max(normal_pdf(x as f32, mu, sigma).ln() as f64);
            //* (1 - s2 / s1) * min_matches * log(s1);
        } else {
            // Treat a1/a2 as two single-end reads
            // 20 corresponds to a value of log(normal_pdf(x, mu, sigma)) of more than 5 stddevs away (for most reasonable values of stddev)
            combined_score = a1.score as f64 + a2.score as f64 - 20.0;
        }

        let aln_pair = ScoredAlignmentPair { score: combined_score, alignment1: a1, alignment2: a2 };
        high_scores.push(aln_pair);
    }

    // Finally, add highest scores of both mates as individually mapped
    // 20 corresponds to  a value of log( normal_pdf(x, mu, sigma ) ) of more than 5 stddevs away (for most reasonable values of stddev)
    let combined_score = a1_indv_max.score as f64 + a2_indv_max.score as f64 - 20;
    high_scores.push_back(ScoredAlignmentPair {score: combined_score, alignment1: a1_indv_max, alignment2: a2_indv_max});

    return high_scores;
}

/*
 * Align a pair of reads for which only one has NAMs. For the other, rescue
 * is attempted by aligning it locally.
 */
fn rescue_read(
    read2: &Read,  // read to be rescued
    read1: &Read,  // read that has NAMs
    aligner: &Aligner,
    references: &[RefSequence],
    nams1: &mut [Nam],
    max_tries: usize,
    dropoff: f32,
    details: &mut [Details; 2],
    k: usize,
    mu: f32,
    sigma: f32
) -> Vec<ScoredAlignmentPair> {
    let n_max1 = nams1[0];
    let mut tries = 0;

    let alignments1 = vec![];
    let alignments2 = vec![];
    for nam in nams1 {
        let score_dropoff1 = nam.n_hits as f32 / n_max1.n_hits as f32;
        // only consider top hits (as minimap2 does) and break if below dropoff cutoff.
        if tries >= max_tries || score_dropoff1 < dropoff {
            break;
        }

        let consistent_nam = reverse_nam_if_needed(nam, read1, references, k);
        details[0].nam_inconsistent += !consistent_nam;
        if let Some(alignment) = extend_seed(aligner, nam, references, read1, consistent_nam) {
            details[0].gapped += alignment.gapped;
            alignments1.emplace_back(alignment);
            details[0].tried_alignment+ +;

            // Force SW alignment to rescue mate
            Alignment
            a2 = rescue_align(aligner, nam, references, read2, mu, sigma, k);
            details[1].mate_rescue += !a2.is_unaligned;
            alignments2.emplace_back(a2);
        }
        tries += 1;
    }
    std::sort(alignments1.begin(), alignments1.end(), by_score<Alignment>);
    std::sort(alignments2.begin(), alignments2.end(), by_score<Alignment>);

    // Calculate best combined score here
    let high_scores = get_best_scoring_pairs(alignments1, alignments2, mu, sigma );

    return high_scores;
}

/// Align a read to the reference given the mapping location of its mate.
/// Return true if rescue by alignment was actually attempted
fn rescue_align(
    aligner: &Aligner,
    mate_nam: &Nam,
    references: &[RefSequence],
    read: &Read,
    mu: f32,
    sigma: f32,
    k: usize,
) -> Alignment {
    let alignment;
    let read_len = read.len();

    let (r_tmp, a, b) =
        if mate_nam.is_revcomp {
            (
                read.seq(),
                mate_nam.projected_ref_start().saturating_sub(mu+5*sigma),
                mate_nam.projected_ref_start() + read_len/2  // at most half read overlap
            )
        } else {
            (
                read.rc(), // mate is rc since fr orientation
                mate_nam.ref_end + (read_len - mate_nam.query_end) - read_len/2,  // at most half read overlap
                mate_nam.ref_end + (read_len - mate_nam.query_end) + (mu+5*sigma)
            )
        };

    let ref_len = references[mate_nam.ref_id].sequence.len();
    let ref_start = std::max(0, std::min(a, ref_len));
    let ref_end = std::min(ref_len, std::max(0, b));

    if ref_end < ref_start + k {
        alignment.cigar = Cigar::new();
        alignment.edit_distance = read_len;
        alignment.score = 0;
        alignment.ref_start =  0;
        alignment.is_revcomp = mate_nam.is_revcomp;
        alignment.reference_id = mate_nam.ref_id;
        alignment.is_unaligned = true;
//        std::cerr << "RESCUE: Caught Bug3! ref start: " << ref_start << " ref end: " << ref_end << " ref len:  " << ref_len << std::endl;
        return Alignment {
            cigar: Cigar::new(),
            alignment.edit_distance = read_len;
            alignment.score = 0;
            alignment.ref_start =  0;
            alignment.is_revcomp = mate_nam.is_revcomp;
            alignment.reference_id = mate_nam.ref_id;
            alignment.is_unaligned = true;
        }
    }
    std::string ref_segm = references.sequences[mate_nam.ref_id].substr(ref_start, ref_end - ref_start);

    if (!has_shared_substring(r_tmp, ref_segm, k)) {
        alignment.cigar = Cigar();
        alignment.edit_distance = read_len;
        alignment.score = 0;
        alignment.ref_start =  0;
        alignment.is_rc = mate_nam.is_rc;
        alignment.ref_id = mate_nam.ref_id;
        alignment.is_unaligned = true;
        return alignment;
    }
    auto info = aligner.align(r_tmp, ref_segm);

    alignment.cigar = info.cigar;
    alignment.edit_distance = info.edit_distance;
    alignment.score = info.sw_score;
    alignment.ref_start = ref_start + info.ref_start;
    alignment.is_rc = !mate_nam.is_rc;
    alignment.ref_id = mate_nam.ref_id;
    alignment.is_unaligned = info.cigar.empty();
    alignment.length = info.ref_span();

    return alignment;
}

/// Determine (roughly) whether the read sequence has some l-mer (with l = k*2/3)
/// in common with the reference sequence
fn has_shared_substring(
    read_seq: &[u8],
    ref_seq: &[u8],
    k: usize,
) -> bool {
    let sub_size = 2 * k / 3;
    let step_size = k / 3;
    let submer;
    for i in (0..read_seq.len().saturating_sub(sub_size)).step_by(step_size) {
        submer = read_seq[i..i+sub_size];
        if memmem::find(ref_seq, submer).is_some() {
            return true;
        }
    }

    false
}

fn get_best_scoring_pairs(
    alignments1: &[Alignment], alignments2: &[Alignment], mu: f32, sigma: f32
) -> Vec<ScoredAlignmentPair> {
    let mut pairs = vec![];
    for a1 in alignments1 {
        for a2 in alignments2 {
            let dist = a1.ref_start.abs_diff(a2.ref_start);
            let mut score = (a1.score + a2.score) as f32;
            if (a1.is_revcomp ^ a2.is_revcomp) && (dist as f32) < mu + 4.0 * sigma {
                score += normal_pdf(dist as f32, mu, sigma).ln();
            }
            else { // individual score
                // 10 corresponds to a value of log(normal_pdf(dist, mu, sigma)) of more than 4 stddevs away
                score -= 10.0;
            }
            pairs.push(ScoredAlignmentPair {score: score as f64, alignment1: a1.clone(), alignment2: a2.clone()});
        }
    }

    pairs
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

/// Find high-scoring NAMs and NAM pairs. Proper pairs are preferred, but also
/// high-scoring NAMs that could not be paired up are returned (these get a
/// "dummy" NAM as partner in the returned vector).
fn get_best_scoring_nam_pairs(
    nams1: &[Nam],
    nams2: &[Nam],
    mu: f32,
    sigma: f32,
) -> Vec<NamPair> {
    let mut nam_pairs = vec![];
    if nams1.is_empty() && nams2.is_empty() {
        return nam_pairs;
    }

    // Find NAM pairs that appear to be proper pairs
    let mut added_n1 = HashSet::new();
    let mut added_n2 = HashSet::new();
    let mut best_joint_hits = 0;
    for nam1 in nams1 {
        for nam2 in nams2 {
            let joint_hits = nam1.n_hits + nam2.n_hits;
            if joint_hits < best_joint_hits / 2 {
                break;
            }
            if is_proper_nam_pair(nam1, nam2, mu, sigma) {
                nam_pairs.push(NamPair{n_hits: joint_hits, nam1: *nam1, nam2: *nam2});
                added_n1.insert(nam1.nam_id);
                added_n2.insert(nam2.nam_id);
                best_joint_hits = joint_hits.max(best_joint_hits);
            }
        }
    }

    // Find high-scoring R1 NAMs that are not part of a proper pair
    let mut dummy_nam;
    dummy_nam.ref_start = -1;
    if (!nams1.empty()) {
        let best_joint_hits1 = if best_joint_hits > 0 { best_joint_hits } else { nams1[0].n_hits };
        for nam1 : nams1 {
            if nam1.n_hits < best_joint_hits1 / 2 {
                break;
            }
            if (added_n1.find(nam1.nam_id) != added_n1.end()) {
                continue;
            }
//            int n1_penalty = std::abs(nam1.query_span() - nam1.ref_span());
            nam_pairs.push_back(NamPair{nam1.n_hits, nam1, dummy_nam});
        }
    }

    // Find high-scoring R2 NAMs that are not part of a proper pair
    if (!nams2.empty()) {
        int best_joint_hits2 = best_joint_hits > 0 ? best_joint_hits : nams2[0].n_hits;
        for (auto &nam2 : nams2) {
            if (nam2.n_hits < best_joint_hits2 / 2) {
                break;
            }
            if (added_n2.find(nam2.nam_id) != added_n2.end()){
                continue;
            }
//            int n2_penalty = std::abs(nam2.query_span() - nam2.ref_span());
            nam_pairs.push_back(NamPair{nam2.n_hits, dummy_nam, nam2});
        }
    }

    std::sort(
        nam_pairs.begin(),
        nam_pairs.end(),
        [](const NamPair& a, const NamPair& b) -> bool { return a.n_hits > b.n_hits; }
    ); // Sort by highest score first

    return nam_pairs;
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

struct NamPair {
    n_hits: usize,
    nam1: Nam,
    nam2: Nam,
}

#[derive(Clone)]
pub struct ScoredAlignmentPair {
    score: f64,
    alignment1: Alignment,
    alignment2: Alignment,
}

/// Remove consecutive identical alignment pairs and leave only the first.
fn deduplicate_scored_pairs(pairs: &mut Vec<ScoredAlignmentPair>) {
    // TODO use Vec::dedup(_by...)
    if pairs.len() < 2 {
        return;
    }
    let mut prev_ref_start1 = pairs[0].alignment1.ref_start;
    let mut prev_ref_start2 = pairs[0].alignment2.ref_start;
    let mut prev_ref_id1 = pairs[0].alignment1.reference_id;
    let mut prev_ref_id2 = pairs[0].alignment2.reference_id;
    let mut j = 1;
    for i in 1..pairs.len() {
        let ref_start1 = pairs[i].alignment1.ref_start;
        let ref_start2 = pairs[i].alignment2.ref_start;
        let ref_id1 = pairs[i].alignment1.reference_id;
        let ref_id2 = pairs[i].alignment2.reference_id;
        if  ref_start1 != prev_ref_start1 ||
            ref_start2 != prev_ref_start2 ||
            ref_id1 != prev_ref_id1 ||
            ref_id2 != prev_ref_id2
        {
            prev_ref_start1 = ref_start1;
            prev_ref_start2 = ref_start2;
            prev_ref_id1 = ref_id1;
            prev_ref_id2 = ref_id2;
            pairs[j] = pairs[i].clone();
            j += 1;
        }
    }
    pairs.truncate(j);
}

fn count_best_alignment_pairs(pairs: &[ScoredAlignmentPair]) -> usize {
    if pairs.is_empty() {
        0
    } else {
        pairs.iter().take_while(|x| x.score == pairs[0].score).count()
    }
}

fn aligned_pairs_to_sam(
    high_scores: &[ScoredAlignmentPair],
    max_secondary: usize,
    secondary_dropoff: f64,
    record1: &SequenceRecord,
    record2: &SequenceRecord,
    read1: &Read,
    read2: &Read,
    mu: f32,
    sigma: f32,
    details: &[Details; 2],
) -> Vec<SamRecord> {

    if high_scores.is_empty() {
        sam.add_unmapped_pair(record1, record2);
        return;
    }

    let (mapq1, mapq2) = joint_mapq_from_high_scores(&high_scores);
    let best_aln_pair = &high_scores[0];


    // append both alignments to string here
    if max_secondary == 0 {
        let alignment1 = best_aln_pair.alignment1;
        let alignment2 = best_aln_pair.alignment2;

        sam.add_pair(alignment1, alignment2, record1, record2, read1.rc, read2.rc, mapq1, mapq2, is_proper_pair(alignment1, alignment2, mu, sigma), true, details);
    } else {
        let max_out = std::min(high_scores.size(), max_secondary);
        let mut is_primary = true;
        let s_max = best_aln_pair.score;
        for i in 0..max_out {
            let aln_pair = high_scores[i];
            let alignment1 = aln_pair.alignment1;
            let alignment2 = aln_pair.alignment2;
            let s_score = aln_pair.score;
            if i > 0 {
                is_primary = false;
                mapq1 = 0;
                mapq2 = 0;
            }
            if s_max - s_score < secondary_dropoff {
                let is_proper = is_proper_pair(alignment1, alignment2, mu, sigma);
                sam.add_pair(alignment1, alignment2, record1, record2, read1.rc, read2.rc, mapq1, mapq2, is_proper, is_primary, details);
            } else {
                break;
            }
        }
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

const INV_SQRT_PI: f32 = 0.3989422804014327;

fn normal_pdf(x: f32, mu: f32, sigma: f32) -> f32 {
    let a = (x - mu) / sigma;

    INV_SQRT_PI / sigma * (-0.5 * a * a).exp()
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
