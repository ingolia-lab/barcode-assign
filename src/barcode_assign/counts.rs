use std::collections::HashMap;
use std::fmt::Display;
use std::io::BufRead;
use std::io::BufReader;
use std::io::Read;
use std::io::Write;
use std::iter::FromIterator;
use std::path::Path;

use bio::io::fasta;
use bio::io::fastq;
use failure;

pub struct SampleCounts(HashMap<Vec<u8>, usize>);

impl SampleCounts {
    pub fn count_map(self) -> HashMap<Vec<u8>, usize> { self.0 }

    pub fn from_file<P: AsRef<Path> + Display>(filename: P) -> Result<Self, failure::Error> {
        Self::read(std::fs::File::open(filename.as_ref())?)
            .map_err(|e| format_err!("Reading file {}: {}", filename, e))
    }

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

            counts.insert(barcode.as_bytes().to_vec(), count);
        }
        Ok(SampleCounts(counts))
    }

    pub fn write<W: Write>(&self, barcode_out: W) -> Result<(), failure::Error> {
        let mut bcout = std::io::BufWriter::new(barcode_out);
        
        for (barcode, count) in self.0.iter() {
            write!(bcout, "{}\t{}\n", String::from_utf8_lossy(barcode), count)?;
        }
        
        Ok(())        
    }

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

    pub fn total_counts<'a, I: Iterator<Item = &'a SampleCounts>>(counts_iter: I) -> Self {
        let mut total_counts = HashMap::new();
        for SampleCounts(count_map) in counts_iter {
            for (barcode, count) in count_map.iter() {
                let barcode_count = total_counts.entry(barcode.to_vec()).or_insert(0);
                *barcode_count += *count;
            }
        }
        SampleCounts(total_counts)
    }

    pub fn barcode_count_vec<'a, I: Iterator<Item = &'a SampleCounts>>(
        counts_iter: I,
        barcode: &'a [u8],
    ) -> Vec<usize> {
        counts_iter
            .map(|SampleCounts(count_table)| count_table.get(barcode).map_or(0, |ct| *ct))
            .collect()
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
