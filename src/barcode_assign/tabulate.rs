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
            let input_counts = Self::read_counts(&input)?;
            counts.push((input.to_string(), input_counts));
        }
        
        let total_counts = Self::total_read_counts(counts.iter().flat_map(|(_input, counts)| counts.iter()));
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
            let count_vec = Self::barcode_count_vec(&counts, barcode);
            
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
            let total: usize = count_vec.iter().map(|ct| *ct).sum();
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
            if count_vec.iter().max().unwrap_or(&0) < &mininsample {
                return true;
            }
        }
        
        false
    }
    
    pub fn barcode_count_vec(counts: &Vec<(String, HashMap<String, usize>)>,
                             barcode: &str)
                             -> Vec<usize>
    {
        let mut cvec = Vec::new();
        for (_input, input_counts) in counts.iter() {
            cvec.push(input_counts.get(barcode).map_or(0, |ct| *ct));
        }
        cvec
    }
    
    pub fn total_read_counts<'a, I: Iterator<Item = (&'a String, &'a usize)>>(read_counts: I) -> HashMap<String, usize> {
        let mut total_counts = HashMap::new();
        for (barcode, count) in read_counts {
            let barcode_count = total_counts.entry(barcode.to_string()).or_insert(0);
            *barcode_count += *count;
        }
        total_counts
    }
    
    pub fn read_counts<P: AsRef<Path> + Display>(input: P) -> Result<HashMap<String,usize>, failure::Error> {
        let mut counts = HashMap::new();
        
        for (line_no, line_res) in BufReader::new(std::fs::File::open(input.as_ref())?).lines().enumerate() {
            let line = line_res?;
            let mut field_iter = line.split("\t");

            let barcode = field_iter
                .next()
                .ok_or_else(|| format_err!("Missing barcode in {} line {}", input, line_no))?;
            let count = field_iter
                .next()
                .ok_or_else(|| format_err!("Missing count in {} line {} barcode {}", input, line_no, barcode))?
                .parse::<usize>()
                .map_err(|e| format_err!("Malformed count in {} line {} barcode {}: {}", input, line_no, barcode, e))?;
            if field_iter.next().is_some() {
                bail!("Extra fields after count in {} line {} barcode {}\n{:?}", input, line_no, barcode, line);
            }
            
            counts.insert(barcode.to_string(), count);
        }
        Ok(counts)
    }
}
