use fastrand::Rng;
use crate::details::Details;
use crate::fasta::RefSequence;
use crate::fastq::SequenceRecord;
use crate::index::StrobemerIndex;
use crate::insertsize::InsertSizeDistribution;
use crate::mapper::{get_best_scoring_nam_pairs, NamPair};
use crate::nam::{get_nams, Nam};
use crate::paf::PafRecord;

/// Map a single-end read to the reference and return PAF records
///
/// This implements "mapping-only" mode in which no base-level alignments are computed
pub fn map_single_end_read(
    record: &SequenceRecord,
    index: &StrobemerIndex,
    references: &[RefSequence],
    rescue_level: usize,
    use_mcs: bool,
    rng: &mut Rng,
) -> (Vec<PafRecord>, Details) {
    let (nam_details, nams) = get_nams(&record.sequence, index, rescue_level, use_mcs, rng);

    if nams.is_empty() {
        (vec![], nam_details.into())
    } else {
        (vec![paf_record_from_nam(&nams[0], &record.name, references, record.sequence.len())], nam_details.into())
    }
}

/// Map a single-end read to the reference and estimate abundances
///
/// This implements abundance estimation mode (`--aemb`)
pub fn abundances_single_end_read(
    record: &SequenceRecord,
    index: &StrobemerIndex,
    abundances: &mut [f64],
    rescue_level: usize,
    use_mcs: bool,
    rng: &mut Rng,
) {
    let (_, nams) = get_nams(&record.sequence, index, rescue_level, use_mcs, rng);
    let n_best = nams.iter().take_while(|nam| nam.score == nams[0].score).count();
    let weight = record.sequence.len() as f64 / n_best as f64;
    for nam in &nams[0..n_best] {
        abundances[nam.ref_id] += weight;
    }
}

/// Convert Nam into PAF record
fn paf_record_from_nam(nam: &Nam, name: &str, references: &[RefSequence], query_length: usize) -> PafRecord {
    PafRecord {
        query_name: name.into(),
        query_length: query_length as u64,
        query_start: nam.query_start as u64,
        query_end: nam.query_end as u64,
        is_revcomp: nam.is_revcomp,
        target_name: references[nam.ref_id].name.clone(),
        target_length: references[nam.ref_id].sequence.len() as u64,
        target_start: nam.ref_start as u64,
        target_end: nam.ref_end as u64,
        n_matches: nam.n_matches as u64,
        alignment_length: (nam.ref_end - nam.ref_start) as u64,
        mapping_quality: None,
    }
}

/// Map a paired-end read pair to the reference and return PAF records
///
/// This implements "mapping-only" mode in which no base-level alignments are computed
pub fn map_paired_end_read(
    r1: &SequenceRecord,
    r2: &SequenceRecord,
    index: &StrobemerIndex,
    references: &[RefSequence],
    rescue_level: usize,
    insert_size_distribution: &mut InsertSizeDistribution,
    use_mcs: bool,
    rng: &mut Rng,
) -> (Vec<PafRecord>, Details) {
    let (mut nam_details1, nams1) = get_nams(&r1.sequence, index, rescue_level, use_mcs, rng);
    let (nam_details2, nams2) = get_nams(&r2.sequence, index, rescue_level, use_mcs, rng);

    let nam_pairs = get_best_scoring_nam_pairs(&nams1, &nams2, insert_size_distribution.mu, insert_size_distribution.sigma);
    let mapped_nam = get_best_paired_map_location(
        &nam_pairs,
        &nams1,
        &nams2,
        insert_size_distribution,
    );
    let mut records = vec![];

    match mapped_nam {
        MappedNams::Individual(nam1, nam2) => {
            if let Some(nam) = nam1 {
                records.push(paf_record_from_nam(&nam, &r1.name, references, r1.sequence.len()))
            }
            if let Some(nam) = nam2 {
                records.push(paf_record_from_nam(&nam, &r2.name, references, r2.sequence.len()))
            }
        }
        MappedNams::Pair(nam1, nam2) => {
            records.push(paf_record_from_nam(&nam1, &r1.name, references, r1.sequence.len()));
            records.push(paf_record_from_nam(&nam2, &r2.name, references, r2.sequence.len()));
        }
        MappedNams::Unmapped => {},
    }
    nam_details1 += nam_details2;
    (records, nam_details1.into())
}

/// Map a paired-end read pair to the reference and estimate abundances
///
/// This implements abundance estimation mode (`--aemb`)
pub fn abundances_paired_end_read(
    r1: &SequenceRecord,
    r2: &SequenceRecord,
    index: &StrobemerIndex,
    abundances: &mut [f64],
    rescue_level: usize,
    insert_size_distribution: &mut InsertSizeDistribution,
    use_mcs: bool,
    rng: &mut Rng,
) {
    let nams1 = get_nams(&r1.sequence, index, rescue_level, use_mcs, rng).1;
    let nams2 = get_nams(&r2.sequence, index, rescue_level, use_mcs, rng).1;

    let nam_pairs = get_best_scoring_nam_pairs(&nams1, &nams2, insert_size_distribution.mu, insert_size_distribution.sigma);
    let mapped_nam = get_best_paired_map_location(
        &nam_pairs,
        &nams1,
        &nams2,
        insert_size_distribution,
    );

    match mapped_nam {
        MappedNams::Pair(nam1, nam2) => {
            let joint_score = nam1.score + nam2.score;
            let n_best = nam_pairs.iter()
                .take_while(|nam_pair| 
                    nam_pair.nam1.as_ref().map_or(0, |nam| nam.score) + 
                    nam_pair.nam2.as_ref().map_or(0, |nam| nam.score) == joint_score)
                .count();
            let weight_r1 = r1.sequence.len() as f64 / n_best as f64;
            let weight_r2 = r2.sequence.len() as f64 / n_best as f64;
            for nam_pair in &nam_pairs[..n_best] {
                if let Some(nam) = &nam_pair.nam1 {
                    abundances[nam.ref_id] += weight_r1;
                } 
                if let Some(nam) = &nam_pair.nam2 {
                    abundances[nam.ref_id] += weight_r2;
                }
            }
        }
        MappedNams::Individual(_, _) => {
            for (nams, read_len) in [(&nams1, r1.sequence.len()), (&nams2, r2.sequence.len())] {
                let n_best = nams.iter().take_while(|nam| nam.score == nams[0].score).count();
                let weight = read_len as f64 / n_best as f64;
                for nam in &nams[0..n_best] {
                    abundances[nam.ref_id] += weight;
                }
            }
        }
        MappedNams::Unmapped => {},
    }
}

enum MappedNams {
    Individual(Option<Nam>, Option<Nam>),
    Pair(Nam, Nam),
    Unmapped,
}

/// Given two lists of NAMs from R1 and R2, find the best location (preferably a proper pair).
/// This is used for mapping-only (PAF) mode and abundances output
fn get_best_paired_map_location(
    nam_pairs: &[NamPair],
    nams1: &[Nam],
    nams2: &[Nam],
    insert_size_distribution: &mut InsertSizeDistribution,
) -> MappedNams {
    if nam_pairs.is_empty() && nams1.is_empty() && nams2.is_empty() {
        return MappedNams::Unmapped;
    }

    // Find first NAM pair that is a proper pair.
    // The first one is also the one with the highest score
    // since nam_pairs is sorted descending by score
    let best_joint_pair = nam_pairs
        .iter()
        .find(|&nam_pair| nam_pair.nam1.is_some() && nam_pair.nam2.is_some());

    let joint_score = if let Some(nam_pair) = best_joint_pair {
        nam_pair.nam1.as_ref().map_or(0, |nam| nam.score) +
        nam_pair.nam2.as_ref().map_or(0, |nam| nam.score)
    } else { 0 };

    // Get individual best scores.
    // nams1 and nams2 are also sorted descending by score.
    let best_individual_nam1 = nams1.first();
    let best_individual_nam2 = nams2.first();

    let individual_score =
        best_individual_nam1.map_or(0, |nam| nam.score) +
        best_individual_nam2.map_or(0, |nam| nam.score);

    // Divisor 2 is penalty for being mapped individually
    if joint_score > individual_score / 2 {
        let best_joint_pair= best_joint_pair.unwrap();
        let best = (best_joint_pair.nam1.clone(), best_joint_pair.nam2.clone());
        if insert_size_distribution.sample_size < 400 {
            insert_size_distribution.update(best.0.as_ref().unwrap().ref_start.abs_diff(best.1.as_ref().unwrap().ref_start));
        }

        // TODO unwrap should not be needed
        MappedNams::Pair(best.0.unwrap(), best.1.unwrap())
    } else {
        MappedNams::Individual(best_individual_nam1.cloned(), best_individual_nam2.cloned())
    }
}
