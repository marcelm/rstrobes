// partially based on ssw_cpp.h/.cpp by Wan-Ping Lee and Mengyao Zhao

mod raw;

use raw::{ssw_init, s_profile};

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

struct SswAlignment<'a> {
    sw_score: u16,           // The best alignment score
    sw_score_next_best: u16, // The next best alignment score
    ref_begin: i32,          // Reference begin position of the best alignment
    ref_end: i32,            // Reference end position of the best alignment
    query_begin: i32,        // Query begin position of the best alignment
    query_end: i32,          // Query end position of the best alignment
    ref_end_next_best: i32,  // Reference end position of the next best alignment
    mismatches: i32,         // Number of mismatches of the alignment
    cigar_string: &'a str,   // Cigar string of the best alignment
    cigar: &'a [u32],        // Cigar stored in the BAM format
                             //   high 28 bits: length
                             //   low 4 bits: M/I/D/S/X (0/1/2/4/8);
}

struct Filter {
  // NOTE: No matter the filter, those five fields of Alignment will be given anyway.
  //       sw_score; sw_score_next_best; ref_end; query_end; ref_end_next_best.
  // NOTE: Only need score of alignments, please set 'report_begin_position'
  //       and 'report_cigar' false.

  bool report_begin_position;    // Give ref_begin and query_begin.
                                 //   If it is not set, ref_begin and query_begin are -1.
  bool report_cigar;             // Give cigar_string and cigar.
                                 //   report_begin_position is automatically TRUE.

  // When *report_cigar* is true and alignment passes these two filters,
  //   cigar_string and cigar will be given.
  uint16_t score_filter;         // score >= score_filter
  uint16_t distance_filter;      // ((ref_end - ref_begin) < distance_filter) &&
                                 // ((query_end - read_begin) < distance_filter)

  Filter()
    : report_begin_position(true)
    , report_cigar(true)
    , score_filter(0)
    , distance_filter(32767)
  {};

  Filter(const bool& pos, const bool& cigar, const uint16_t& score, const uint16_t& dis)
    : report_begin_position(pos)
    , report_cigar(cigar)
    , score_filter(score)
    , distance_filter(dis)
    {};
}

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


fn doit() {
    ssw_aligner(StripedSmithWaterman::Aligner(parameters.match, parameters.mismatch, parameters.gap_open, parameters.gap_extend));
    flag = ssw_aligner.Align(query.c_str(), ref.c_str(), ref.size(), filter, &alignment_ssw, maskLen);
    if flag != 0 {
        aln.edit_distance = 100000;
        aln.ref_start = 0;
        aln.sw_score = -100000;
        return aln;
    }
}

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

    query.iter().map(|c| TABLE[c as usize]).collect()
}

struct Profile<'a> {
    profile: *mut s_profile,
    query: &'a [i8],
    score_matrix: &'a [i8]
}

impl Profile {
    /// query must be translated
    fn new(query: &[i8], score_matrix: &[i8]) -> Self {
        // TODO should return an error if query.is_empty()
        let score_size = 2;
        let profile = unsafe {
            raw::ssw_init(query.as_ptr(), query.len() as i32, score_matrix.as_ptr(), score_matrix.len() as i32, score_size)
        };
        Profile { profile, query, score_matrix }
    }
}

impl Drop for Profile {
    fn drop(&mut self) {
        unsafe { raw::init_destroy(self.profile); }
    }
}

pub struct SswAligner {
    score_matrix: *const i8,
    match_score: u8,
    mismatch_penalty: u8,
    gap_open_penalty: u8,
    gap_extend_penalty: u8,
    ssw_profile: *Profile,
}

impl SswAligner {
    pub fn new(match_score: u8, mismatch_penalty: u8, gap_open_penalty: u8, gap_extend_penalty: u8) -> Self {
        /*
    : score_matrix_(NULL)
    , score_matrix_size_(5)
    , translation_matrix_(NULL)
    , translated_reference_(NULL)
    , reference_length_(0)
         */
        let mat = match_score as i8;
        let mis = -mismatch_penalty as i8;
        let score_matrix: [[i8; 5]; 5] = [
            [mat, mis, mis, mis, mis],
            [mis, mat, mis, mis, mis],
            [mis, mis, mat, mis, mis],
            [mis, mis, mis, mat, mis],
            [mis, mis, mis, mis, mis],
        ];
        SswAligner { score_matrix: score_matrix as *const i8, match_score, mismatch_penalty, gap_open_penalty, gap_extend_penalty, ssw_profile }
    }

    pub fn align(query: &[u8], refseq: &[u8]) -> Option<SswAlignment> {
        if query.is_empty() {
            return None;
        }

        let query = translate(query);
        let refseq = translate(refseq);
        // 1. translate query
        // 2. translate ref

        //uint16_t Aligner::Align(const char* query, const char* ref, const int& ref_len,
        //            const Filter& filter, Alignment* alignment, const int32_t maskLen) const
/*        if !translation_matrix_ return false;

  int8_t* translated_query = new int8_t[query_len];
  TranslateBase(query, query_len, translated_query);

  // calculate the valid length
  int valid_ref_len = ref_len;
  int8_t* translated_ref = new int8_t[valid_ref_len];
  TranslateBase(ref, valid_ref_len, translated_ref);


  const int8_t score_size = 2;
  */
  s_profile* profile = ssw_init(translated_query, query_len, score_matrix_, score_matrix_size_, score_size);

  uint8_t flag = 0;
  SetFlag(filter, &flag);
  s_align* s_al = ssw_align(profile, translated_ref, valid_ref_len,
                                 static_cast<int>(gap_opening_penalty_),
				 static_cast<int>(gap_extending_penalty_),
				 flag, filter.score_filter, filter.distance_filter, maskLen);

  alignment->Clear();
  ConvertAlignment(*s_al, query_len, alignment);
  alignment->mismatches = CalculateNumberMismatch(&*alignment, translated_ref, translated_query, query_len);
  uint16_t align_flag = s_al->flag;

  // Free memory
  delete [] translated_query;
  delete [] translated_ref;
  align_destroy(s_al);
  init_destroy(profile);

  return align_flag;
}


    }
        translation_matrix_ = new int8_t[SizeOfArray(kBaseTranslation)];
        memcpy(translation_matrix_, kBaseTranslation, sizeof(int8_t) * SizeOfArray(kBaseTranslation));

        unsafe { ssw_init(query, query_length, mat, n, score_size); }
    score_matrix: *const i8,
    match_score: u8,
    mismatch_penalty: u8,
    gap_open_penalty: u8,
    gap_extend_penalty: u8,
    ssw_profile: *SswProfile,

    }

}
