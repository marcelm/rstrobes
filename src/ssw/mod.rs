// partially based on ssw_cpp.h/.cpp by Wan-Ping Lee and Mengyao Zhao

use std::marker::PhantomData;
use crate::cigar::Cigar;

mod raw;

const ENCODED_OPS: [u8; 128] = [
	0,         0,         0,         0,
	0,         0,         0,         0,
	0,         0,         0,         0,
	0,         0,         0,         0,
	0,         0,         0,         0,
	0,         0,         0,         0,
	0,         0,         0,         0,
	0,         0,         0,         0,
	0 /*   */, 0 /* ! */, 0 /* " */, 0 /* # */,
	0 /* $ */, 0 /* % */, 0 /* & */, 0 /* ' */,
	0 /* ( */, 0 /* ) */, 0 /* * */, 0 /* + */,
	0 /* , */, 0 /* - */, 0 /* . */, 0 /* / */,
	0 /* 0 */, 0 /* 1 */, 0 /* 2 */, 0 /* 3 */,
	0 /* 4 */, 0 /* 5 */, 0 /* 6 */, 0 /* 7 */,
	0 /* 8 */, 0 /* 9 */, 0 /* : */, 0 /* ; */,
	0 /* < */, 7 /* = */, 0 /* > */, 0 /* ? */,
	0 /* @ */, 0 /* A */, 0 /* B */, 0 /* C */,
	2 /* D */, 0 /* E */, 0 /* F */, 0 /* G */,
	5 /* H */, 1 /* I */, 0 /* J */, 0 /* K */,
	0 /* L */, 0 /* M */, 3 /* N */, 0 /* O */,
	6 /* P */, 0 /* Q */, 0 /* R */, 4 /* S */,
	0 /* T */, 0 /* U */, 0 /* V */, 0 /* W */,
	8 /* X */, 0 /* Y */, 0 /* Z */, 0 /* [ */,
	0 /* \ */, 0 /* ] */, 0 /* ^ */, 0 /* _ */,
	0 /* ` */, 0 /* a */, 0 /* b */, 0 /* c */,
	0 /* d */, 0 /* e */, 0 /* f */, 0 /* g */,
	0 /* h */, 0 /* i */, 0 /* j */, 0 /* k */,
	0 /* l */, 0 /* m */, 0 /* n */, 0 /* o */,
	0 /* p */, 0 /* q */, 0 /* r */, 0 /* s */,
	0 /* t */, 0 /* u */, 0 /* v */, 0 /* w */,
	0 /* x */, 0 /* y */, 0 /* z */, 0 /* { */,
	0 /* | */, 0 /* } */, 0 /* ~ */, 0 /*  */
];

const MAPSTR: &[u8] = b"MIDNSHP=X";

/*const BAM_CIGAR_SHIFT: u32 = 4;

fn to_cigar_int(length: u32, op_letter: c_char) -> u32 {
    (length << BAM_CIGAR_SHIFT) | (ENCODED_OPS[op_letter])
}

fn cigar_int_to_op(cigar_int: u32) -> u8 {
	if (cigar_int & 0xf) > 8 {
        b'M'
    } else {
        MAPSTR[(cigar_int & 0xf) as u8]
    }
}

fn cigar_int_to_len(cigar_int: u32) -> u32 {
	cigar_int >> BAM_CIGAR_SHIFT
}*/


  // =========
  // @function Align the query againt the reference.
  //           [NOTICE] The reference won't replace the reference
  //                      set by SetReferenceSequence.
  // @param    query     The query sequence.
  // @param    ref       The reference sequence.
  //                     [NOTICE] It is not necessary null terminated.
  // @param    ref_len   The length of the reference sequence.
  // @param    filter    The filter for the alignment.
  // @param    alignment The container contains the result.
  // @param    maskLen   The distance between the optimal and suboptimal alignment ending position will >= maskLen. We suggest to
  //                     use readLen/2, if you don't have special concerns. Note: maskLen has to be >= 15, otherwise this function
  //                     will NOT return the suboptimal alignment information.
  // @return   If the alignment path is accurate (or has missing part). 0: accurate; 1: banded_sw is totally failed; 2: banded_sw returned path has missing part
  // =========
  // uint16_t Align(const char* query, const char* ref, const int& ref_len,
  //            const Filter& filter, Alignment* alignment, const int32_t maskLen) const;

/*fn doit() {
    ssw_aligner(StripedSmithWaterman::Aligner(parameters.match, parameters.mismatch, parameters.gap_open, parameters.gap_extend));
    flag = ssw_aligner.Align(query.c_str(), ref.c_str(), ref.size(), filter, &alignment_ssw, maskLen);
    if flag != 0 {
        aln.edit_distance = 100000;
        aln.ref_start = 0;
        aln.sw_score = -100000;
        return aln;
    }
}*/

fn translate(query: &[u8]) -> Vec<i8> {
    static TABLE: [i8; 128] = [
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
    ];

    query.iter().map(|c| TABLE[*c as usize]).collect()
}

#[derive(Debug)]
struct Profile<'a> {
    profile: *mut raw::s_profile,
    query: &'a [i8],
    score_matrix: &'a [i8]
}

#[derive(Debug)]
pub struct SswAlignment {
    pub sw_score: u16,           // The best alignment score
    pub sw_score_next_best: u16, // The next best alignment score
    pub ref_begin: i32,          // Reference begin position of the best alignment
    pub ref_end: i32,            // Reference end position of the best alignment
    pub query_begin: i32,        // Query begin position of the best alignment
    pub query_end: i32,          // Query end position of the best alignment
    pub ref_end_next_best: i32,  // Reference end position of the next best alignment
    pub mismatches: i32,         // Number of mismatches of the alignment
    //pub cigar_string: &'a str,   // Cigar string of the best alignment
    pub cigar: Cigar,
    //pub cigar: &'a [u32],        // Cigar stored in the BAM format
                             //   high 28 bits: length
                             //   low 4 bits: M/I/D/S/X (0/1/2/4/8);
}

struct Alignment<'a> {
    raw: *mut raw::s_align,
    _align: PhantomData<&'a raw::s_align>,
}

impl<'a> Alignment<'a> {
    fn is_valid(&self) -> bool {
        unsafe { (*self.raw).flag == 0 }
    }
}

impl<'a> Drop for Alignment<'a> {
    fn drop(&mut self) {
        unsafe { raw::align_destroy(self.raw); }
    }
}

impl<'a> From<Alignment<'a>> for SswAlignment {
    fn from(alignment: Alignment) -> SswAlignment {
        let raw = unsafe { alignment.raw.as_ref().expect("hm") };
        let cigar_slice = unsafe { std::slice::from_raw_parts(raw.cigar, raw.cigar_length as usize) };
        let cigar = Cigar::try_from(cigar_slice).expect("Invalid CIGAR");


            /*al->cigar.clear();
  al->cigar_string.clear();

  if (s_al.cigarLen > 0) {
    std::ostringstream cigar_string;
    if (al->query_begin > 0) {
      uint32_t cigar = to_cigar_int(al->query_begin, 'S');
      al->cigar.push_back(cigar);
      cigar_string << al->query_begin << 'S';
    }

    for (int i = 0; i < s_al.cigarLen; ++i) {
      al->cigar.push_back(s_al.cigar[i]);
      cigar_string << cigar_int_to_len(s_al.cigar[i]) << cigar_int_to_op(s_al.cigar[i]);
    }

    int end = query_len - al->query_end - 1;
    if (end > 0) {
      uint32_t cigar = to_cigar_int(end, 'S');
      al->cigar.push_back(cigar);
      cigar_string << end << 'S';
    }

  } // end if
}
*/
/*  ConvertAlignment(*s_al, query_len, alignment);
  alignment->mismatches = CalculateNumberMismatch(&*alignment, translated_ref, translated_query, query_len);
*/

        SswAlignment {
            sw_score: raw.score1,
            sw_score_next_best: raw.score2,
            ref_begin: raw.ref_begin1,
            ref_end: raw.ref_end1,
            query_begin: raw.read_begin1,
            query_end: raw.read_end1,
            ref_end_next_best: raw.ref_end2,
            cigar,
            mismatches: 0,  // TODO
        }
    }
}

impl<'a> Profile<'a> {
    /// query must be translated
    fn new(translated_query: &'a [i8], score_matrix: &'a [i8]) -> Self {
        // TODO should return an error if query.is_empty()
        let score_size = 2;
        let profile = unsafe {
            raw::ssw_init(translated_query.as_ptr(), translated_query.len() as i32, score_matrix.as_ptr(), score_matrix.len() as i32, score_size)
        };
        dbg!("Profile::new");
        Profile { profile, query: translated_query, score_matrix }
    }

    fn align(&self, translated_refseq: &[i8], gap_open_penalty: u8, gap_extend_penalty: u8, flag: u8, score_filter: u16, distance_filter: i32, mask_len: i32) -> Alignment {
        dbg!("align");
        let alignment = unsafe {
            raw::ssw_align(self.profile, translated_refseq.as_ptr(), translated_refseq.len() as i32, gap_open_penalty, gap_extend_penalty, flag, score_filter, distance_filter, mask_len)
        };
        unsafe { dbg!(&*alignment); }
        Alignment { raw: alignment, _align: PhantomData }
    }
}


impl<'a> Drop for Profile<'a> {
    fn drop(&mut self) {
        dbg!("Calling drop for Profile");
        unsafe { raw::init_destroy(self.profile); }
    }
}

pub struct Aligner {
    score_matrix: Vec<i8>,
    match_score: u8,
    mismatch_penalty: u8,
    gap_open_penalty: u8,
    gap_extend_penalty: u8,
//    ssw_profile: *Profile,
}

impl Aligner {
    pub fn new(match_score: u8, mismatch_penalty: u8, gap_open_penalty: u8, gap_extend_penalty: u8) -> Self {
        let mat = match_score as i8;
        let mis = -(mismatch_penalty as i8);
        let score_matrix = vec![
            mat, mis, mis, mis, mis,
            mis, mat, mis, mis, mis,
            mis, mis, mat, mis, mis,
            mis, mis, mis, mat, mis,
            mis, mis, mis, mis, mis,
        ];
        Aligner { score_matrix, match_score, mismatch_penalty, gap_open_penalty, gap_extend_penalty }
    }

    pub fn align(&self, query: &[u8], refseq: &[u8]) -> Option<SswAlignment> {
        if query.is_empty() {
            return None;
        }

        let query = translate(query);
        let refseq = translate(refseq);
        let profile = Profile::new(&query, &self.score_matrix);
        dbg!(&profile);
        let flag = 0x0f;
        let score_filter = 0;
        let distance_filter = i32::MAX;
        let mask_len = std::cmp::max(query.len() / 2, 15);

        let alignment = profile.align(&refseq, self.gap_open_penalty, self.gap_extend_penalty, flag, score_filter, distance_filter, mask_len as i32);
        if !alignment.is_valid() {
            return None;
        }

        Some(SswAlignment::from(alignment))
    }
}
