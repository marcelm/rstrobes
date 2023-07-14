use std::cmp::{max, min};
use std::collections::VecDeque;
use std::{env, io};
use std::collections::hash_map::RandomState;
use std::fs::File;
use std::hash::Hasher;
use std::io::{BufReader, BufWriter, BufRead, Error, Write};
use std::path::Path;
use clap::Parser;
use fxhash;
use rstrobes::fastq::FastqReader;
use rstrobes::strobes::{SyncmerIterator,RandstrobeIterator,RandstrobeParameters};
use strobes::strobes::SyncmerIterator;

#[derive(Parser, Debug)]
#[command(long_about = None)]
struct Args {

    /// Print syncmers instead of randstrobes
    //#[arg(long, default_value_t = false)]
    //syncmers: bool,

    /// Path to input FASTA
    ref_path: String,

    fastq_path: String,
}

#[derive(Debug)]
struct FastaRecord {
    name: String,
    sequence: Vec<u8>,
}

fn read_fasta<R: BufRead>(reader: &mut R) -> Result<Vec<FastaRecord>, Error> {
    let mut records = Vec::<FastaRecord>::new();
    let mut name = String::new();
    let mut sequence = Vec::new();
    let mut has_record = false;
    for line in reader.lines() {
        let line = line.unwrap();
        let line = line.as_bytes();
        if line.is_empty() {
            continue;
        }
        if line[0] == b'>' {
            if has_record {
                records.push(FastaRecord{name, sequence});
            }
            name = String::from_utf8(line[1..].to_vec()).unwrap();
            if let Some(i) = name.find(|c: char| c.is_ascii_whitespace()) {
                name = name[..i].to_string();
            }
            sequence = Vec::new();
            has_record = true;
        } else {
            sequence.extend(line);
        }
    }
    if has_record {
        records.push(FastaRecord{name, sequence});
    }

    Ok(records)
}

struct SyncmerParameters {
    k: u8,
    s: u8,
    t: u8,
}

impl SyncmerParameters {
    pub fn new(k: u8, s: u8) -> Self {
        SyncmerParameters { k, s, t: (k - s) / 2 + 1 }
    }

// TODO
// void verify() const {
//     if (k <= 7 || k > 32) {
//         throw BadParameter("k not in [8,32]");
//     }
//     if (s > k) {
//         throw BadParameter("s is larger than k");
//     }
//     if ((k - s) % 2 != 0) {
//         throw BadParameter("(k - s) must be an even number to create canonical syncmers. Please set s to e.g. k-2, k-4, k-6, ...");
//     }
// }
}

struct RandstrobeParameters {
    l: u32,
    u: u32,
    q: u64,
    max_dist: usize,
    w_min: u32,
    w_max: u32,

// TODO
// void verify() const {
//     if (max_dist > 255) {
//         throw BadParameter("maximum seed length (-m <max_dist>) is larger than 255");
//     }
// }
}

/* Settings that influence index creation */
struct IndexParameters {
    canonical_read_length: u32,
    syncmer: SyncmerParameters,
    randstrobe: RandstrobeParameters,
}

impl IndexParameters {
    pub fn new(canonical_read_length: u32, k: u8, s: u8, l: u32, u: u32, q: u64, max_dist: usize) -> Self {
        let w_min = max(1, (k / (k - s + 1) + l)) as u32;
        let w_max = (k / (k - s + 1) + u) as u32;
        IndexParameters {
            canonical_read_length,
            syncmer: SyncmerParameters::new(k, s),
            randstrobe: RandstrobeParameters { l, u, q, max_dist, w_min, w_max },
        }
    }
}

fn main() -> Result<(), Error> {
    let args = Args::parse();
    let path = Path::new(&args.ref_path);
    let f = File::open(path)?;
    let mut reader = BufReader::new(f);
    let records = read_fasta(&mut reader).unwrap();
    let k = 20usize;
    let s = 16usize;

    // IndexParameters(r=150, k=20, s=16, l=1, u=7, q=255, max_dist=80, t_syncmer=3, w_min=5, w_max=11)
    let parameters = IndexParameters::new(150, 20, 16, 1, 7, 255, 80);
    debug_assert_eq!(parameters.w_min, 5);
    debug_assert_eq!(parameters.w_max, 11);

    let index = StrobemerIndex::new(references, index_parameters, bits);

    index.populate(opt.f, opt.n_threads);

    Ok(())
}

fn dumpstrobes(parameters: RandstrobeParameters) {
    let mut writer = BufWriter::new(io::stdout());

    for record in &records {
        let name = &record.name;
        let mut syncmers = SyncmerIterator::new(&record.sequence, k, s, 3);
        for randstrobe in RandstrobeIterator::new(&mut syncmers, &parameters) {
            writeln!(writer, "{}\t{}\t{}", name, randstrobe.strobe1_pos, randstrobe.strobe2_pos + k)?;
        }
    }
}

#[test]
fn verify_cli() {
    use clap::CommandFactory;
    Args::command().debug_assert()
}
