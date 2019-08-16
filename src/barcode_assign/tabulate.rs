use std::collections::HashMap;
use std::fmt::Display;
use std::io::BufRead;
use std::io::BufReader;
use std::io::Read;
use std::io::Write;
use std::path::Path;

use failure;

pub struct CLI {
    pub inputs: Vec<String>,
    pub output: String,
    pub mintotal: Option<usize>,
    pub minsamples: Option<usize>,
    pub mininsample: Option<usize>,
    pub omitfile: Option<String>
}

impl CLI {
    pub fn run(&self) -> Result<(), failure::Error> {
        let mut counts = Vec::new();
        
        for input in self.inputs.iter() {
            let input_counts = SampleCounts::from_file(&input)?;
            counts.push((input.to_string(), input_counts));
        }
        
        let total_counts = SampleCounts::total_counts(counts.iter().map(|(_input, counts)| counts));
        let mut barcodes: Vec<(String, usize)> = total_counts.into_iter().collect();
        barcodes.sort_by_key(|(_barcode, counts)| -(*counts as isize));

        let mut out = std::fs::File::create(&self.output)?;
        write!(out, "barcode")?;
        for (input, _input_counts) in counts.iter() {
            write!(out, "\t{}", input)?;
        }
        write!(out, "\n")?;

        let mut omit_out: Box<dyn Write> = match &self.omitfile {
            Some(f) => Box::new(std::fs::File::create(f)?),
            None => Box::new(std::io::sink()),
        };
        
        for (barcode, _counts) in barcodes.iter() {
            let count_vec = SampleCounts::barcode_count_vec(counts.iter().map(|(_input, counts)| counts), barcode);
            
            if self.is_omitted(&count_vec) {
                write!(omit_out, "{}\n", barcode)?;
            } else {
                write!(out, "{}", barcode)?;
                for ct in count_vec.iter() {
                    write!(out, "\t{}", ct)?;
                }
                write!(out, "\n")?;
            }
        }

        Ok(())
    }

    pub fn is_omitted(&self, count_vec: &Vec<usize>) -> bool {
        if let Some(mintotal) = self.mintotal {
            let total: usize = count_vec.iter().copied().sum();
            if total < mintotal {
                return true;
            }
        }

        if let Some(minsamples) = self.minsamples {
            let nsamples: usize = count_vec.iter().map(|ct| if *ct > 0 { 1 } else { 0 }).sum();
            if nsamples < minsamples {
                return true;
            }
        }

        if let Some(mininsample) = self.mininsample {
            if count_vec.iter().max().copied().unwrap_or(0) < mininsample {
                return true;
            }
        }
        
        false
    }
    
}

pub struct SampleCounts (HashMap<String, usize>);

impl SampleCounts {
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
                .map_err(|e| format_err!("Malformed count line {} barcode {}: {}", line_no, barcode, e))?;
            if field_iter.next().is_some() {
                bail!("Extra fields after count line {} barcode {}\n{:?}", line_no, barcode, line);
            }
            
            counts.insert(barcode.to_string(), count);
        }
        Ok( SampleCounts(counts) )
    }

    pub fn total_counts<'a, I: Iterator<Item = &'a SampleCounts>>(counts_iter: I) -> Self {
        let mut total_counts = HashMap::new();
        for SampleCounts (count_map) in counts_iter {
            for (barcode, count) in count_map.iter() {
                let barcode_count = total_counts.entry(barcode.to_string()).or_insert(0);
                *barcode_count += *count;
            }
        }
        SampleCounts(total_counts)
    }

    pub fn barcode_count_vec<'a, I: Iterator<Item = &'a SampleCounts>>(counts_iter: I, barcode: &'a str) -> Vec<usize>
    {
        counts_iter.map(|SampleCounts(count_table)| count_table.get(barcode).map_or(0, |ct| *ct)).collect()
    }
    
}

impl IntoIterator for SampleCounts
{
    type Item = (String, usize);
    type IntoIter = ::std::collections::hash_map::IntoIter<String, usize>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}
