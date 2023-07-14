use std::fmt;
use std::io::{BufRead, BufReader, BufWriter, Read, Result, Write};
use std::str;

#[derive(Debug, Clone)]
pub struct SequenceRecord {
    pub name: String,
    pub sequence: Vec<u8>,
    pub qualities: Vec<u8>,
}

impl SequenceRecord {
    pub fn new(name: String, sequence: Vec<u8>, qualities: Vec<u8>) -> Self {
        SequenceRecord { name, sequence, qualities }
    }

    pub fn len(&self) -> usize {
        self.sequence.len()
    }
}

impl fmt::Display for SequenceRecord {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let s = str::from_utf8(&self.sequence).unwrap();
        let q = str::from_utf8(&self.qualities).unwrap();
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
        if qualities.len() > 0 && qualities[qualities.len() - 1] == b'\r' {
            qualities.pop();
        }
        assert_eq!(sequence.len(), qualities.len());
        Some(Ok(SequenceRecord {
            name: name.to_string(),
            sequence,
            qualities,
        }))
    }
}

pub struct FastqWriter<W: Write> {
    writer: BufWriter<W>,
}

impl<W: Write> FastqWriter<W> {
    pub fn new(writer: W) -> Self {
        FastqWriter {
            writer: BufWriter::new(writer),
        }
    }

    pub fn write_record(&mut self, record: &SequenceRecord) {
        self.writer.write(b"@").unwrap();
        self.writer.write(&record.name.as_bytes()).unwrap();
        self.writer.write(b"\n").unwrap();
        self.writer.write(&record.sequence).unwrap();
        self.writer.write(b"\n+\n").unwrap();
        self.writer.write(&record.qualities).unwrap();
        self.writer.write(b"\n").unwrap();
        //write!(self.writer, "@{}\n{:?}\n+\n{:?}\n", record.name, record.sequence, record.qualities);
    }
}
