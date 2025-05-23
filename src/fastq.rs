use std::collections::VecDeque;
use std::fmt;
use std::io::{BufRead, BufReader, Read, Result};
use std::str;
use crate::fasta::split_header;
use crate::io::xopen;

#[derive(Debug, Clone)]
pub struct SequenceRecord {
    pub name: String,
    pub comment: Option<String>,
    pub sequence: Vec<u8>,
    pub qualities: Option<Vec<u8>>,
}

impl SequenceRecord {
    pub fn new(name: String, comment: Option<String>, sequence: Vec<u8>, qualities: Option<Vec<u8>>) -> Self {
        SequenceRecord { name, comment, sequence, qualities }
    }

    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }
}

impl fmt::Display for SequenceRecord {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let s = str::from_utf8(&self.sequence).unwrap();
        let e = vec![];
        let q = str::from_utf8(self.qualities.as_ref().unwrap_or(&e)).unwrap();
        write!(
            f,
            "name={}, length={}, sequence={}, qualities={}",
            self.name,
            self.len(),
            s,
            q
        )
    }
}

#[derive(Debug)]
pub struct PeekableSequenceReader<R: Read + Send> {
    fastq_reader: FastqReader<R>,
    buffer: VecDeque<SequenceRecord>,
}

impl<R: Read + Send> PeekableSequenceReader<R> {
    pub fn new(reader: R) -> Self {
        let fastq_reader = FastqReader::new(reader);
        PeekableSequenceReader {
            fastq_reader,
            buffer: VecDeque::new(),
        }
    }

    // Retrieve up to n records
    pub fn peek(&mut self, n: usize) -> Result<Vec<SequenceRecord>> {
        while self.buffer.len() < n {
            if let Some(item) = self.fastq_reader.next() {
                self.buffer.push_back(item?)
            } else {
                break;
            }
        }

        Ok(self.buffer.make_contiguous().to_vec())
    }
}

impl<R: Read + Send> Iterator for PeekableSequenceReader<R> {
    type Item = Result<SequenceRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(item) = self.buffer.pop_front() {
            Some(Ok(item))
        } else {
            self.fastq_reader.next()
        }
    }
}


#[derive(Debug)]
pub struct FastqReader<R: Read> {
    reader: BufReader<R>,
    err: bool,
}

impl<R: Read> FastqReader<R> {
    pub fn new(reader: R) -> FastqReader<R> {
        FastqReader {
            reader: BufReader::new(reader),
            err: false,
        }
    }
}

impl<R: Read> Iterator for FastqReader<R> {
    type Item = Result<SequenceRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.err {
            return None;
        }
        let mut name = String::new();
        match self.reader.read_line(&mut name) {
            Ok(0) => {
                return None;
            }
            Ok(_) => {}
            Err(e) => {
                self.err = true;
                return Some(Err(e));
            }
        }
        let name = name[1..].trim_end();
        let (name, comment) = split_header(name);

        if name.is_empty() {
            //self.err = true;
            //return Some(Err("error"));
        }

        let mut sequence = Vec::new();
        self.reader.read_until(b'\n', &mut sequence).unwrap();
        sequence.pop();
        if !sequence.is_empty() && sequence[sequence.len() - 1] == b'\r' {
            sequence.pop();
        }

        let mut name2 = String::new();
        self.reader.read_line(&mut name2).unwrap();

        let mut qualities = Vec::new();
        self.reader.read_until(b'\n', &mut qualities).unwrap();
        qualities.pop();
        if !qualities.is_empty() && qualities[qualities.len() - 1] == b'\r' {
            qualities.pop();
        }
        assert_eq!(sequence.len(), qualities.len());
        Some(Ok(SequenceRecord {
            name,
            comment,
            sequence: sequence.iter().map(|&c| c.to_ascii_uppercase()).collect(),
            qualities: Some(qualities),
        }))
    }
}

type RecordPair = (SequenceRecord, Option<SequenceRecord>);

/// Iterate over paired-end *or* single-end reads
pub fn record_iterator<'a, R: Read + Send + 'a>(
    fastq_reader1: PeekableSequenceReader<R>, path_r2: Option<&str>
) -> Result<Box<dyn Iterator<Item=Result<RecordPair>> + Send + 'a>> 
{
    if let Some(r2_path) = path_r2 {
        let fastq_reader2 = FastqReader::new(xopen(r2_path)?);
        Ok(
            Box::new(
                fastq_reader1.zip(fastq_reader2).map(
                    |p|
                        match p {
                            (Ok(r1), Ok(r2)) => { Ok((r1, Some(r2))) }
                            (Err(e), _) => Err(e),
                            (_, Err(e)) => Err(e),
                        }
                )
            )
        )
    } else {
        Ok(
            Box::new(
                fastq_reader1.map(
                    |r| r.map(|sr| (sr, None))
                )
            )
        )
    }
}


struct InterleavedIterator<R: Read + Send> {
    fastq_reader: PeekableSequenceReader<R>,
    next_record: Option<SequenceRecord>,
}

impl<R: Read + Send> InterleavedIterator<R> {
    pub fn new(fastq_reader: PeekableSequenceReader<R>) -> Self {
        InterleavedIterator { fastq_reader, next_record: None }
    }
}

impl<R: Read + Send> Iterator for InterleavedIterator<R> { 
    type Item = Result<RecordPair>;
    
    fn next(&mut self) -> Option<<Self as Iterator>::Item> {
        let record1 =    
            if let Some(record) = self.next_record.take() {
                record
            } else if let Some(record) = self.fastq_reader.next() {
                match record {
                    Ok(record) => { record }
                    Err(e) => { return Some(Err(e)); }
                }           
            } else {
                return None;
            };
        
        match self.fastq_reader.next() {
            Some(Err(e)) => { Some(Err(e)) }
            Some(Ok(record2)) => {
                if record1.name != record2.name {
                    self.next_record = Some(record2);
                    
                    Some(Ok((record1, None)))
                } else {
                    Some(Ok((record1, Some(record2))))
                }
            }
            None => Some(Ok((record1, None))),
        }
    }
}

pub fn interleaved_record_iterator<'a, R: Read + Send + 'a>(
    fastq_reader: PeekableSequenceReader<R>,
) -> Box<dyn Iterator<Item=Result<RecordPair>> + Send + 'a>
{
    Box::new(InterleavedIterator::new(fastq_reader))
}

#[cfg(test)]
mod tests {
    use std::fs::File;
    use std::io::Result;
    use crate::fastq::{interleaved_record_iterator, PeekableSequenceReader, RecordPair};

    #[test]
    fn test_peekable_sequence_reader() {
        let f = File::open("tests/phix.1.fastq").unwrap();
        let mut reader = PeekableSequenceReader::new(f);

        assert_eq!(reader.next().unwrap().unwrap().name, "SRR1377138.1");
        assert_eq!(
            reader.peek(2).unwrap().iter().map(|record| record.name.clone()).collect::<Vec<String>>(),
            vec![String::from("SRR1377138.2"), String::from("SRR1377138.3/1")]
        );
        assert_eq!(reader.next().unwrap().unwrap().name, "SRR1377138.2");
    }

    #[test]
    fn test_interleaved_record_iterator() {
        let f = File::open("tests/interleaved.fq").unwrap();
        let reader = PeekableSequenceReader::new(f);
        let it = interleaved_record_iterator(reader);
        
        let record_pairs: Vec<RecordPair> = it.collect::<Result<Vec<RecordPair>>>().unwrap();
        
        assert_eq!(record_pairs.len(), 6);
        assert_eq!(record_pairs[0].0.name, "SRR4052021.2");
        assert_eq!(record_pairs[0].1.as_ref().unwrap().name, "SRR4052021.2");
        assert!(record_pairs[0].0.sequence.starts_with(b"GTCGCCCA"));
        assert!(record_pairs[0].1.as_ref().unwrap().sequence.starts_with(b"ATGTATTA"));
        
        assert_eq!(record_pairs[1].0.name, "SRR4052021.3");
        assert_eq!(record_pairs[1].1.as_ref().unwrap().name, "SRR4052021.3");
        
        assert_eq!(record_pairs[2].0.name, "SRR4052021.13852607");
        assert!(record_pairs[2].1.is_none());
    }
}