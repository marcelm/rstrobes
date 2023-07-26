use std::fmt::{Display, Formatter};
use crate::cigar::Cigar;

const PAIRED: u16 = 1;
const PROPER_PAIR: u16 = 2;
const UNMAP: u16 = 4;
const MUNMAP: u16 = 8;
const REVERSE: u16 = 0x10;
const MREVERSE: u16 = 0x20;
const READ1: u16 = 0x40;
const READ2: u16 = 0x80;
const SECONDARY: u16 = 0x100;
const QCFAIL: u16 = 0x200;
const DUP: u16 = 0x400;
const SUPPLEMENTARY: u16 = 0x800;

pub struct SamRecord {
    query_name: String,
    flags: u16,
    reference_name: Option<String>,
    /// 0-based position, converted to 1-based for output
    pos: Option<u32>,
    mapq: u8,
    cigar: Option<Cigar>,
    mate_reference_name: Option<String>,
    mate_pos: Option<u32>,
    template_len: Option<i32>,
    query_sequence: Option<Vec<u8>>,
    qual: Option<Vec<u8>>,
    edit_distance: u32,
    aln_score: u32,
    // details: Details,
}

impl Display for SamRecord {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        // SAM columns
        //
        // 1 QNAME query template name
        // 2 FLAG  bitwise flag
        // 3 RNAME reference sequence name
        // 4 POS   1-based leftmost mapping position (set to 0 for unmapped read w/o pos)
        // 5 MAPQ  mapping quality (255 if unavailable)
        // 6 CIGAR
        // 7 RNEXT reference name of mate/next read
        // 8 PNEXT position of mate/next read
        // 9 TLEN template length
        // 10 SEQ
        // 11 QUAL
        // 12 optional fields

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
        let cigar = match &self.cigar {
            Some(cigar) => cigar.to_string(),
            None => "*".to_string(),
        };
        write!(f, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.query_name,
            self.flags,
            self.reference_name.unwrap_or("*".to_string()),
            pos,
            self.mapq,
            cigar,
            self.mate_reference_name.unwrap_or("*".to_string()),
            mate_pos,
            self.template_len.unwrap_or(0),
            query_sequence,
            qual,
        )
    }
}
