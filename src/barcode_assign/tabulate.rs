use std::collections::HashMap;
use std::fmt::Display;
use std::io::BufRead;
use std::io::BufReader;
use std::io::Read;
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
        
        unimplemented!()
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
