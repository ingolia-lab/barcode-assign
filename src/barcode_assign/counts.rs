use std::collections::HashMap;
use std::fmt::Display;
use std::io::BufRead;
use std::io::BufReader;
use std::io::Read;
use std::io::Write;
use std::iter::FromIterator;
use std::iter::Sum;
use std::path::Path;

use bio::io::fasta;
use bio::io::fastq;
use failure;

/// Tabulation of barcode counts in a sample
#[derive(Debug,Clone)]
pub struct SampleCounts(HashMap<Vec<u8>, usize>);

impl SampleCounts {
    /// Returns a `HashMap` of barcodes and their counts
    pub fn count_map(self) -> HashMap<Vec<u8>, usize> { self.0 }

    /// Reads a barcode count table from a file. 
    ///
    /// # Arguments
    ///
    /// * `filename` is the name of the file that will be read.
    pub fn from_file<P: AsRef<Path> + Display>(filename: P) -> Result<Self, failure::Error> {
        Self::read(std::fs::File::open(filename.as_ref())?)
            .map_err(|e| format_err!("Reading file {}: {}", filename, e))
    }

    /// Reads and parses a barcode count table.
    ///
    /// A barcode count table is a tab-delimited table of barcode<tab>count.
    ///
    /// # Arguments
    ///
    /// * `input` is an input source that will be read.
    pub fn read<R: Read>(input: R) -> Result<Self, failure::Error> {
        let mut counts = HashMap::new();

        for (line_no, line_res) in BufReader::new(input).lines().enumerate() {
            let line = line_res?;
            let mut field_iter = line.split("\t");

            let barcode = field_iter
                .next()
                .ok_or_else(|| format_err!("Missing barcode line {}", line_no))?;
            let count = field_iter
                .next()
                .ok_or_else(|| format_err!("Missing count line {} barcode {}", line_no, barcode))?
                .parse::<usize>()
                .map_err(|e| {
                    format_err!(
                        "Malformed count line {} barcode {}: {}",
                        line_no,
                        barcode,
                        e
                    )
                })?;
            if field_iter.next().is_some() {
                bail!(
                    "Extra fields after count line {} barcode {}\n{:?}",
                    line_no,
                    barcode,
                    line
                );
            }

            if counts.insert(barcode.as_bytes().to_vec(), count).is_some() {
                return Err(format_err!("Duplicate entry line {} barcode {}", line_no, barcode));
            }
        }
        Ok(SampleCounts(counts))
    }

    /// Writes a barcode count table.
    ///
    /// # Arguments
    ///
    /// * The table will be written to the output desntiation `barcode_out`
    pub fn write<W: Write>(&self, barcode_out: W) -> Result<(), failure::Error> {
        let mut bcout = std::io::BufWriter::new(barcode_out);
        
        for (barcode, count) in self.0.iter() {
            write!(bcout, "{}\t{}\n", String::from_utf8_lossy(barcode), count)?;
        }
        
        Ok(())        
    }

    /// Writes a barcode frequency table.
    ///
    /// The frequency table lists the number of barcodes seen 1 time,
    /// 2 times, etc. The format of the table is
    /// times-seen<tab>number-of-barcodes.
    ///
    /// # Arguments
    ///
    /// * The table will be written to the output desntiation `barcode_out`
    pub fn write_freq_table<W: Write>(&self, freq_out: W) -> Result<(), failure::Error> {
        let mut fout = std::io::BufWriter::new(freq_out);
        
        let mut freq_counts = HashMap::new();
        
        for freq in self.0.values() {
            let freq_count = freq_counts.entry(freq).or_insert(0);
            *freq_count += 1;
        }
        
        let mut freqs: Vec<usize> = freq_counts.keys().map(|&&k| k).collect();
        freqs.sort();
        
        for freq in freqs {
            write!(fout, "{}\t{}\n", freq, freq_counts.get(&freq).unwrap_or(&0))?;
        }
        
        Ok(())
    }

    /// Returns the count for a barcode in the sample.
    ///
    /// If the barcode does not occur, then 0 is returned.
    ///
    /// # Arguments
    ///
    /// `barcode` specifies the barcode to look up.
    pub fn barcode_count(&self, barcode: &[u8]) -> usize {
        self.0.get(barcode).copied().unwrap_or(0)
    }

    /// Returns a vector of counts for a barcode across a set of samples.
    ///
    /// The vector of counts is the same length as the iterator, with
    /// `0` for samples where the barcode is absent.
    ///
    /// # Arguments
    ///
    /// * `counts_iter` is an iterator over sample count data
    /// * `barcode` specifies the barcode to count.
    pub fn barcode_count_vec<'a, I: Iterator<Item = &'a SampleCounts>>(
        counts_iter: I,
        barcode: &'a [u8],
    ) -> Vec<usize> {
        counts_iter
            .map(|counts| counts.barcode_count(barcode))
            .collect()
    }
}

impl <'a> Sum<&'a SampleCounts> for SampleCounts {
    fn sum<I>(iter: I) -> Self
        where I: Iterator<Item = &'a SampleCounts>
    {
        let mut total_counts = HashMap::new();
        for SampleCounts(count_map) in iter {
            for (barcode, count) in count_map.iter() {
                let barcode_count = total_counts.entry(barcode.to_vec()).or_insert(0);
                *barcode_count += *count;
            }
        }
        SampleCounts(total_counts)
    }
}

impl IntoIterator for SampleCounts {
    type Item = (Vec<u8>, usize);
    type IntoIter = ::std::collections::hash_map::IntoIter<Vec<u8>, usize>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl <'a> FromIterator<&'a str> for SampleCounts {
    fn from_iter<I>(iter: I) -> SampleCounts
        where I: IntoIterator<Item = &'a str>
    {
        let mut barcode_counts = HashMap::new();
        
        for bc in iter {
            if barcode_counts.contains_key(bc.as_bytes()) {
                *(barcode_counts.get_mut(bc.as_bytes()).unwrap()) += 1;
            } else {
                barcode_counts.insert(bc.as_bytes().to_vec(), 1);
            }
        }

        SampleCounts(barcode_counts)
    }
}

impl <'a> FromIterator<&'a [u8]> for SampleCounts {
    fn from_iter<I>(iter: I) -> SampleCounts
        where I: IntoIterator<Item = &'a [u8]>
    {
        let mut barcode_counts = HashMap::new();
        
        for bc in iter {
            if barcode_counts.contains_key(bc) {
                *(barcode_counts.get_mut(bc).unwrap()) += 1;
            } else {
                barcode_counts.insert(bc.to_vec(), 1);
            }
        }
        
        SampleCounts(barcode_counts)
    }
}

impl <'a> FromIterator<&'a fastq::Record> for SampleCounts {
    fn from_iter<I>(iter: I) -> SampleCounts
        where I: IntoIterator<Item = &'a fastq::Record>
    {
        iter.into_iter().map(fastq::Record::seq).collect()
    }
}

impl FromIterator<fastq::Record> for SampleCounts {
    fn from_iter<I>(iter: I) -> SampleCounts
        where I: IntoIterator<Item = fastq::Record>
    {
        let mut barcode_counts = HashMap::new();
        
        for rec in iter {
            if barcode_counts.contains_key(rec.seq()) {
                *(barcode_counts.get_mut(rec.seq()).unwrap()) += 1;
            } else {
                barcode_counts.insert(rec.seq().to_vec(), 1);
            }
        }

        SampleCounts(barcode_counts)
    }
}

impl <'a> FromIterator<&'a fasta::Record> for SampleCounts {
    fn from_iter<I>(iter: I) -> SampleCounts
        where I: IntoIterator<Item = &'a fasta::Record>
    {
        iter.into_iter().map(fasta::Record::seq).collect()
    }
}

impl FromIterator<fasta::Record> for SampleCounts {
    fn from_iter<I>(iter: I) -> SampleCounts
        where I: IntoIterator<Item = fasta::Record>
    {
        let mut barcode_counts = HashMap::new();
        
        for rec in iter {
            if barcode_counts.contains_key(rec.seq()) {
                *(barcode_counts.get_mut(rec.seq()).unwrap()) += 1;
            } else {
                barcode_counts.insert(rec.seq().to_vec(), 1);
            }
        }

        SampleCounts(barcode_counts)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn from_iterator() {
        let bc1 = b"ACGTACGT".to_vec();
        let bc2 = b"CATGCATG".to_vec();
        let bc3 = b"TACGTACG".to_vec();

        let ct1: SampleCounts = vec![bc1.as_slice(), bc2.as_slice(), bc1.as_slice()].into_iter().collect();
        assert_eq!(ct1.barcode_count(&bc1), 2);
        assert_eq!(ct1.barcode_count(&bc2), 1);
        assert_eq!(ct1.barcode_count(&bc3), 0);

        let ct2: SampleCounts = vec![bc2.as_slice(), bc1.as_slice(), bc3.as_slice(), bc3.as_slice()].into_iter().collect();
        assert_eq!(ct2.barcode_count(&bc1), 1);
        assert_eq!(ct2.barcode_count(&bc2), 1);
        assert_eq!(ct2.barcode_count(&bc3), 2);

        let bc1s = "ACGTACGT".to_string();
        let bc2s = "CATGCATG".to_string();

        let ct1s: SampleCounts = vec![bc1s.as_str(), bc2s.as_str(), bc1s.as_str()].into_iter().collect();
        assert_eq!(ct1s.barcode_count(&bc1), 2);
        assert_eq!(ct1s.barcode_count(&bc2), 1);
        assert_eq!(ct1s.barcode_count(&bc3), 0);
    }

    #[test]
    fn from_fasta() {
        let fastafile = b">one\nACGTA\n>two\nCGTAC\n>three\nACGTA\n>four\nGTACG\n>five\nACGTA\n".to_vec();
        let fain = fasta::Reader::new(fastafile.as_slice());
        let counts_res: Result<SampleCounts, std::io::Error>
            = fain.records().collect();
        let counts = counts_res.unwrap();
        let mut cvec: Vec<(Vec<u8>, usize)> = counts.into_iter().collect();
        cvec.sort();
        assert_eq!(cvec, vec![(b"ACGTA".to_vec(), 3), (b"CGTAC".to_vec(), 1), (b"GTACG".to_vec(), 1)]);
    }
    
    #[test]
    fn from_fastq() {
        let fastafile = b"@one\nACGTA\n+\n~~~~~\n@two\nCGTAC\n+\n~~~~~\n@three\nGTACG\n+\n~~~~~\n@four\nCGTAC\n+\n~~~~~\n".to_vec();
        let fqin = fastq::Reader::new(fastafile.as_slice());
        let counts_res: Result<SampleCounts, std::io::Error>
            = fqin.records().collect();
        let counts = counts_res.unwrap();
        let mut cvec: Vec<(Vec<u8>, usize)> = counts.into_iter().collect();
        cvec.sort();
        assert_eq!(cvec, vec![(b"ACGTA".to_vec(), 1), (b"CGTAC".to_vec(), 2), (b"GTACG".to_vec(), 1)]);
    }
    
    #[test]
    fn read_table() {
        let table = "TACGGA\t3\nCAGTA\t2\nAATTA\t6\n";
        let cts = SampleCounts::read(table.as_bytes()).unwrap();
        assert_eq!(cts.barcode_count(b"TACGGA"), 3);
        assert_eq!(cts.barcode_count(b"CAGTA"), 2);
        assert_eq!(cts.barcode_count(b"AATTA"), 6);
        assert_eq!(cts.barcode_count(b"ACTGA"), 0);
    }

    #[test]
    fn read_written() {
        let bc1 = b"ACGTACGT".to_vec();
        let bc2 = b"CATGCATG".to_vec();
        let bc3 = b"TACGTACG".to_vec();
        let bc4 = b"GACTAGTC".to_vec();
        
        let cts: SampleCounts = vec![bc1.as_slice(), bc2.as_slice(), bc1.as_slice(), bc3.as_slice(), bc2.as_slice(), bc1.as_slice()].into_iter().collect();
        assert_eq!(cts.barcode_count(&bc1), 3);
        assert_eq!(cts.barcode_count(&bc2), 2);
        assert_eq!(cts.barcode_count(&bc3), 1);
        assert_eq!(cts.barcode_count(&bc4), 0);

        let mut table_out: Vec<u8> = Vec::new();
        cts.write(&mut table_out).unwrap();

        assert_eq!(cts.into_iter().count(), 3);
        assert_eq!(BufReader::new(table_out.as_slice()).lines().count(), 3);

        let ctsnew = SampleCounts::read(table_out.as_slice()).unwrap();
        assert_eq!(ctsnew.barcode_count(&bc1), 3);
        assert_eq!(ctsnew.barcode_count(&bc2), 2);
        assert_eq!(ctsnew.barcode_count(&bc3), 1);
        assert_eq!(ctsnew.barcode_count(&bc4), 0);
        assert_eq!(ctsnew.into_iter().count(), 3);
    }

    #[test]
    fn count_vec() {
        let table1 = "ACGTA\t3\nCGTAC\t7\nGTACG\t4\n";
        let table2 = "ACGTA\t4\nTACGT\t5\n";
        let table3 = "ACGTA\t5\nCGTAC\t2\nTGCAT\t8\n";

        let cts1 = SampleCounts::read(table1.as_bytes()).unwrap();
        let cts2 = SampleCounts::read(table2.as_bytes()).unwrap();
        let cts3 = SampleCounts::read(table3.as_bytes()).unwrap();

        let ctses = vec![cts1, cts2, cts3];

        assert_eq!(SampleCounts::barcode_count_vec(ctses.iter(), b"ACGTA"),
                   vec![3, 4, 5]);
        assert_eq!(SampleCounts::barcode_count_vec(ctses.iter(), b"CGTAC"),
                   vec![7, 0, 2]);
        assert_eq!(SampleCounts::barcode_count_vec(ctses.iter(), b"GTACG"),
                   vec![4, 0, 0]);
        assert_eq!(SampleCounts::barcode_count_vec(ctses.iter(), b"TACGT"),
                   vec![0, 5, 0]);
        assert_eq!(SampleCounts::barcode_count_vec(ctses.iter(), b"TGCAT"),
                   vec![0, 0, 8]);
        assert_eq!(SampleCounts::barcode_count_vec(ctses.iter(), b"GACTG"),
                   vec![0, 0, 0]);
    }

    #[test]
    fn sum() {
        let table1 = "ACGTA\t3\nCGTAC\t7\nGTACG\t4\n";
        let table2 = "ACGTA\t4\nTACGT\t5\n";
        let table3 = "ACGTA\t5\nCGTAC\t2\nTGCAT\t8\n";

        let cts1 = SampleCounts::read(table1.as_bytes()).unwrap();
        let cts2 = SampleCounts::read(table2.as_bytes()).unwrap();
        let cts3 = SampleCounts::read(table3.as_bytes()).unwrap();
        
        let ctssum: SampleCounts = vec![cts1, cts2, cts3].iter().sum();

        assert_eq!(ctssum.barcode_count(b"ACGTA"), 12);
        assert_eq!(ctssum.barcode_count(b"CGTAC"), 9);
        assert_eq!(ctssum.barcode_count(b"GTACG"), 4);
        assert_eq!(ctssum.barcode_count(b"TACGT"), 5);
        assert_eq!(ctssum.barcode_count(b"TGCAT"), 8);
        assert_eq!(ctssum.barcode_count(b"GACTG"), 0);
        assert_eq!(ctssum.into_iter().count(), 5);
    }

    #[test]
    fn into_iterator() {
        let table1 = "ACGTA\t3\nCGTAC\t7\nGTACG\t4\n";
        let cts1 = SampleCounts::read(table1.as_bytes()).unwrap();
        let mut cvec: Vec<(Vec<u8>, usize)> = cts1.into_iter().collect();
        cvec.sort();
        assert_eq!(cvec, vec![(b"ACGTA".to_vec(), 3), (b"CGTAC".to_vec(), 7), (b"GTACG".to_vec(), 4)]);
    }
}
