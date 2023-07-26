use std::cmp::{max, min, Reverse};
use crate::cigar::Cigar;
use crate::fasta::RefSequence;
use crate::index::{IndexParameters, StrobemerIndex};
use crate::nam::{find_nams, Nam, reverse_nam_if_needed};
use crate::revcomp::reverse_complement;
use crate::strobes::RandstrobeIterator;
use crate::syncmers::SyncmerIterator;
use crate::sam::SamRecord;
use crate::read::Read;

enum CigarMode {
    M, EQX
}

struct MappingParameters {
    r: usize, // { 150 },
    max_secondary: usize,
    dropoff_threshold: f32,
    rescue_level: usize,
    max_tries: usize,
    // rescue_cutoff,
    cigar_mode: CigarMode,
    output_unmapped: bool,
}

impl MappingParameters {
    pub fn new() -> Self {
        MappingParameters {
            r: 150,
            max_secondary: 0,
            dropoff_threshold: 0.5,
            rescue_level: 2,
            max_tries: 20,
            // rescue_cutoff: ...,
            cigar_mode: CigarMode::M,
            output_unmapped: true,
            // details: false,
        }
    }
}

struct Alignment {
    reference_id: usize,
    ref_start: usize,
    cigar: Cigar,
    edit_distance: usize,
    score: usize,
    mapq: u8,
    length: usize,
    is_revcomp: bool,
    is_unaligned: bool,
    /// Whether a gapped alignment function was used to obtain this alignment
    /// (even if true, the alignment can still be without gaps)
    gapped: bool,
}

impl Alignment {
    fn new() -> Self {
        Alignment {

        }
    }
}

pub fn map_single_end_read(
    seq: &Vec<u8>,
    name: &String,
    index: &StrobemerIndex,
    references: &Vec<RefSequence>,
    mapping_parameters: &MappingParameters
) -> Vec<SamRecord> {
    //Details details;
    //Timer strobe_timer;
    let query_randstrobes = randstrobes_query(&seq, &index.parameters);
    //statistics.tot_construct_strobemers += strobe_timer.duration();

    // Timer nam_timer;
    let (nonrepetitive_fraction, mut nams) = find_nams(&query_randstrobes, index);
    println!("nonrepetitive fraction: {}. nams.len(): {}", nonrepetitive_fraction, nams.len());
    for nam in &nams {
        println!("{:?}", nam);
    }

    // statistics.tot_find_nams += nam_timer.duration();
    /*
    TODO
    if (map_param.R > 1) {
        Timer rescue_timer;
        if (nams.empty() || nonrepetitive_fraction < 0.7) {
            details.nam_rescue = true;
            nams = find_nams_rescue(query_randstrobes, index, map_param.rescue_cutoff);
        }
        statistics.tot_time_rescue += rescue_timer.duration();
    }
    details.nams = nams.size();
    Timer nam_sort_timer;
    */
    nams.sort_by_key(|&k| k.score);
    // statistics.tot_sort_nams += nam_sort_timer.duration();

    // Timer extend_timer;
    // align_SE(
    //     aligner, sam, nams, record, index.parameters.syncmer.k,
    //     references, details, map_param.dropoff_threshold, map_param.maxTries,
    //     map_param.max_secondary
    // );

    if nams.is_empty() {
        return Vec::new();
    }
    let sam_records = Vec::new();

    let mut alignments = Vec::new();
    let mut tries = 0;
    let nam_max = nams[0];
    let mut best_edit_distance = u32::MAX;
    let mut best_score = 0;
    let mut second_best_score = 0;

    let mut best_alignment = Alignment::new();

    let k = index.parameters.syncmer.k;
    let read = Read::new(seq);

    best_alignment.is_unaligned = true;
    for nam in &mut nams {
        let score_dropoff = nam.n_hits as f32 / nam_max.n_hits as f32;

        // TODO iterate over slice of nams instead of tracking tries
        if tries >= mapping_parameters.max_tries || (tries > 1 && best_edit_distance == 0) || score_dropoff < mapping_parameters.dropoff_threshold {
            break;
        }
        let consistent_nam = reverse_nam_if_needed(nam, &read, &references, k);
        // details.nam_inconsistent += !consistent_nam;
        let alignment = get_alignment(&aligner, nam, &references, &read, consistent_nam);
        // details.tried_alignment += 1;
        // details.gapped += sam_aln.gapped;

        if alignment.score > best_score {
            second_best_score = best_score;
            best_score = alignment.score;
            best_alignment = alignment.clone();
            if mapping_parameters.max_secondary == 0 {
                best_edit_distance = best_alignment.global_edit_distance;
            }
        } else if alignment.score > second_best_score {
            second_best_score = alignment.score;
        }
        if mapping_parameters.max_secondary > 0 {
            alignments.push(alignment);
        }
        tries += 1;
    }
    let mapq = (60.0 * (best_score - second_best_score) as f32 / best_score as f32) as u8;

    if mapping_parameters.max_secondary == 0 {
        best_alignment.mapq = mapq;
        sam.add(best_alignment, seq, name, read.rc(), true); // TODO details
        return;
    }
    // Highest score first
    alignments.sort_by_key(|&k| Reverse(*k));

    let max_out = min(alignments.len(), mapping_parameters.max_secondary + 1);
    for i in 0..max_out {
        let is_primary = i == 0;
        let alignment = &alignments[i];
        if alignment.score - best_score > 2 * aligner.parameters.mismatch + aligner.parameters.gap_open {
            break;
        }
        if is_primary {
            alignment.mapq = mapq;
        } else {
            alignment.mapq = 255;
        }
        // sam_records.push(SamRecord { ... });
        sam.add(alignment, seq, name, read.rc(), is_primary); // TODO details
    }

    // statistics.tot_extend += extend_timer.duration();
    // statistics += details;
    sam_records
}

/*
 Extend a NAM so that it covers the entire read and return the resulting
 alignment.
*/
fn get_alignment(
    aligner: &Aligner,
    nam: &Nam,
    references: &Vec<RefSequence>,
    read: &Read,
    consistent_nam: bool,
) -> Alignment {
    let query = if nam.is_revcomp { read.rc() } else { read.seq() };
    let refseq = &references[nam.ref_id].sequence;

    let projected_ref_start = max(0, nam.ref_start - nam.query_start);
    let projected_ref_end = min(nam.ref_end + query.len() - nam.query_end, refseq.len());

    let info;
    let result_ref_start;
    let gapped = true;
    if projected_ref_end - projected_ref_start == query.len() && consistent_nam {
        let ref_segm_ham = &refseq[projected_ref_start..projected_ref_start+query.len()];
        let hamming_dist = hamming_distance(query, ref_segm_ham);

        if hamming_dist >= 0 && (((float) hamming_dist / query.size()) < 0.05) {
            // Hamming distance worked fine, no need to ksw align
            info = hamming_align(query, ref_segm_ham, aligner.parameters.match, aligner.parameters.mismatch, aligner.parameters.end_bonus);
            result_ref_start = projected_ref_start + info.ref_start;
            gapped = false;
        }
    }
    if (gapped) {
        const int diff = std::abs(nam.ref_span() - nam.query_span());
        const int ext_left = std::min(50, projected_ref_start);
        const int ref_start = projected_ref_start - ext_left;
        const int ext_right = std::min(std::size_t(50), ref.size() - nam.ref_end);
        const auto ref_segm_size = read.size() + diff + ext_left + ext_right;
        const auto ref_segm = ref.substr(ref_start, ref_segm_size);
        info = aligner.align(query, ref_segm);
        result_ref_start = ref_start + info.ref_start;
    }
    int softclipped = info.query_start + (query.size() - info.query_end);
    Alignment alignment;
    alignment.cigar = std::move(info.cigar);
    alignment.edit_distance = info.edit_distance;
    alignment.global_ed = info.edit_distance + softclipped;
    alignment.score = info.sw_score;
    alignment.ref_start = result_ref_start;
    alignment.length = info.ref_span();
    alignment.is_rc = nam.is_rc;
    alignment.is_unaligned = false;
    alignment.ref_id = nam.ref_id;
    alignment.gapped = gapped;

    alignment
}

#[derive(Debug)]
pub struct QueryRandstrobe {
    pub hash: u64,
    pub start: usize,
    pub end: usize,
    pub is_reverse: bool,
}

/// Generate randstrobes for a query sequence and its reverse complement.
pub fn randstrobes_query(seq: &Vec<u8>, parameters: &IndexParameters) -> Vec<QueryRandstrobe> {
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
                    is_reverse
                }
            );
        }
    }

    randstrobes
}
