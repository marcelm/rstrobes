use std::collections::VecDeque;
use crate::index::REF_RANDSTROBE_HASH_MASK;
use crate::syncmers::{Syncmer, SyncmerIterator};

pub const DEFAULT_AUX_LEN: u8 = 17;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RandstrobeParameters {
    pub w_min: usize,
    pub w_max: usize,
    pub q: u64,
    pub max_dist: u8,
    
    /// Mask for bits of the hash that represent the main hash
    pub main_hash_mask: u64,

    // TODO ensure aux_len <= 63
// TODO
// void verify() const {
//     if (max_dist > 255) {
//         throw BadParameter("maximum seed length (-m <max_dist>) is larger than 255");
//     }
//     if (w_min > w_max) {
//         throw BadParameter("w_min is greater than w_max (choose different -l/-u parameters)");
//     }
// }
}

#[derive(Debug, PartialEq, Eq)]
pub struct Randstrobe {
    pub hash: u64,
    pub strobe1_pos: usize,
    pub strobe2_pos: usize,
}

impl Randstrobe {
    pub fn from_strobes(strobe1: Syncmer, strobe2: Syncmer, main_hash_mask: u64) -> Self {
        Randstrobe {
            hash: Randstrobe::hash(strobe1.hash, strobe2.hash, main_hash_mask),
            strobe1_pos: strobe1.position,
            strobe2_pos: strobe2.position,
        }
    }

    pub fn hash(hash1: u64, hash2: u64, main_hash_mask: u64) -> u64 {
        ((hash1 & main_hash_mask) | (hash2 & !main_hash_mask)) & REF_RANDSTROBE_HASH_MASK
    }
}

pub struct RandstrobeIterator<'a> {
    parameters: &'a RandstrobeParameters,
    syncmers: VecDeque<Syncmer>,
    syncmer_iterator: &'a mut SyncmerIterator<'a>,
}

impl<'a> RandstrobeIterator<'a> {
    pub fn new(syncmer_iterator: &'a mut SyncmerIterator<'a>, parameters: &'a RandstrobeParameters) -> RandstrobeIterator<'a> {
        RandstrobeIterator {
            parameters,
            syncmers: VecDeque::<Syncmer>::new(),
            syncmer_iterator,
        }
    }
}

impl<'a> Iterator for RandstrobeIterator<'a> {
    type Item = Randstrobe;
    fn next(&mut self) -> Option<Self::Item> {
        while self.syncmers.len() <= self.parameters.w_max {
            if let Some(syncmer) = self.syncmer_iterator.next() {
                self.syncmers.push_back(syncmer);
            } else {
                break;
            }
        }
        if self.syncmers.len() <= self.parameters.w_min {
            return None;
        }
        let strobe1 = self.syncmers[0];
        let max_position = strobe1.position + self.parameters.max_dist as usize;
        let mut min_val = u64::MAX;
        let mut strobe2 = self.syncmers[0]; // Defaults if no nearby syncmer

        for i in self.parameters.w_min .. self.syncmers.len() {
            debug_assert!(i <= self.parameters.w_max);
            if self.syncmers[i].position > max_position {
                break;
            }
            let b = (strobe1.hash ^ self.syncmers[i].hash) & self.parameters.q;
            let ones = b.count_ones() as u64;
            if ones < min_val {
                min_val = ones;
                strobe2 = self.syncmers[i];
            }
        }
        self.syncmers.pop_front();

        Some(Randstrobe::from_strobes(strobe1, strobe2, self.parameters.main_hash_mask))
    }
}

#[cfg(test)]
mod test {
    use std::fs::File;
    use std::io::BufReader;
    use crate::fasta::{read_fasta, RefSequence};
    use crate::index::IndexParameters;
    use super::*;

    fn read_phix() -> RefSequence {
        let f = File::open("tests/phix.fasta").unwrap();
        let mut reader = BufReader::new(f);

        read_fasta(&mut reader).unwrap().first().unwrap().clone()
    }

    #[test]
    fn test_randstrobe_iterator() {
        let refseq = read_phix().sequence;
        let parameters = IndexParameters::default_from_read_length(300);
        let mut syncmer_iter = SyncmerIterator::new(&refseq, parameters.syncmer.k, parameters.syncmer.s, parameters.syncmer.t);
        let randstrobe_iter = RandstrobeIterator::new(&mut syncmer_iter, &parameters.randstrobe);

        for randstrobe in randstrobe_iter {
            assert!(randstrobe.hash > 0);
            assert!(randstrobe.strobe2_pos >= randstrobe.strobe2_pos);
            assert!(randstrobe.strobe1_pos < refseq.len());
            assert!(randstrobe.strobe2_pos < refseq.len());
        }
    }

    // This tests an assumption that we need to hold for count_randstrobes()
    #[test]
    fn test_syncmer_and_randstrobe_iterator_same_count() {
        let refseq = read_phix().sequence;
        let parameters = IndexParameters::default_from_read_length(100);
        let syncmer_iter = SyncmerIterator::new(&refseq, parameters.syncmer.k, parameters.syncmer.s, parameters.syncmer.t);
        let syncmer_count = syncmer_iter.count();

        let mut syncmer_iter = SyncmerIterator::new(&refseq, parameters.syncmer.k, parameters.syncmer.s, parameters.syncmer.t);
        let randstrobe_iter = RandstrobeIterator::new(&mut syncmer_iter, &parameters.randstrobe);
        let randstrobe_count = randstrobe_iter.count();

        assert_eq!(randstrobe_count + parameters.randstrobe.w_min, syncmer_count);
    }
}
