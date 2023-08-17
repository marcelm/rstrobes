use std::ffi::c_char;

/*
   @typedef	structure of the alignment result
   @field	score1	the best alignment score
   @field	score2	sub-optimal alignment score
   @field	ref_begin1	0-based best alignment beginning position on reference;	ref_begin1 = -1 when the best alignment beginning
                    position is not available
   @field	ref_end1	0-based best alignment ending position on reference
   @field	read_begin1	0-based best alignment beginning position on read; read_begin1 = -1 when the best alignment beginning
                    position is not available
   @field	read_end1	0-based best alignment ending position on read
   @field	ref_end2	0-based sub-optimal alignment ending position on reference
   @field	cigar	best alignment cigar; stored the same as that in BAM format, high 28 bits: length, low 4 bits: M/I/D (0/1/2);
                cigar = 0 when the best alignment path is not available
   @field	cigarLen	length of the cigar string; cigarLen = 0 when the best alignment path is not available
   @field  flag  If the alignment path is accurate (or has missing part). 0: accurate; 1: banded_sw is totally failed; 2: banded_sw returned path has missing part
 */

#[repr(C)]
struct SswAlignment {
    score1: u16,
    score2: u16,
    ref_begin1: i32,
    ref_end1: i32,
    read_begin1: i32,
    read_end1: i32,
    ref_end2: i32,
    cigar: *const u32,
    cigar_length: i32,
    flag: u16,
}

pub struct SswProfile {
    _data: [u8; 0],
    _marker:
        core::marker::PhantomData<(*mut u8, core::marker::PhantomPinned)>,
}

#[link(name = "ssw")]
extern {
    fn ssw_init(read: *const i8, read_length: i32, mat: *const i8, n: i32, score_size: i8) -> *mut SswProfile;

    fn init_destroy(p: *mut SswProfile);

    fn ssw_align(
        prof: *const SswProfile,
        refseq: *const i8,
        ref_length: i32,
        weight_gap_open: u8,
        weight_gap_extend: u8,
        flag: u8,
        filters: u16,
        filterd: i32,
        mask_length: i32
    ) -> *const SswAlignment;

    fn align_destroy(a: *mut SswAlignment);
}
