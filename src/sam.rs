use std::fmt::{Display, Formatter};
use crate::cigar::Cigar;

struct SamRecord {
    query_name: String,
    flags: u32,
    reference_name: Option<String>,
    pos: Option<u32>,
    mapq: u8,
    cigar: Cigar,
    mate_reference_name: Option<String>,
    mate_pos: Option<u32>,
    template_len: Option<u32>,
    query_sequence: Option<Vec<u8>>,
    qual: Option<Vec<u8>>,
    edit_distance: u32,
    aln_score: u32,
    // Details details;
}

impl Display for SamRecord {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let mate_pos = match self.mate_pos {
            Some(pos) => pos + 1,
            None => 0,
        };
        let query_sequence = match &self.query_sequence {
            Some(seq) if self.flags & SECONDARY != 0 => std::str::from_utf8(seq).unwrap(),
            _ => "*",
        };
        let qual = match &self.qual {
            Some(qual) if self.flags & SECONDARY != 0 => std::str::from_utf8(qual).unwrap(),
            _ => "*",
        };
        let pos = match self.pos {
            Some(pos) => pos + 1,
            None => 0,
        };
        write!(f, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.query_name,
            self.flags,
            self.reference_name.unwrap_or("*".to_string()),
            pos,
            self.mapq,
            self.cigar,
            self.mate_reference_name.unwrap_or("*".to_string()),
            mate_pos,
            self.template_len.unwrap_or(0),
            query_sequence,
            qual,
        )
    }
}
