use std::collections::HashMap;
use std::fmt::Display;
use std::io::BufRead;
use std::io::BufReader;
use std::io::Read;
use std::path::Path;

use failure;

pub struct SampleCounts(HashMap<String, usize>);

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

            counts.insert(barcode.to_string(), count);
        }
        Ok(SampleCounts(counts))
    }

    pub fn total_counts<'a, I: Iterator<Item = &'a SampleCounts>>(counts_iter: I) -> Self {
        let mut total_counts = HashMap::new();
        for SampleCounts(count_map) in counts_iter {
            for (barcode, count) in count_map.iter() {
                let barcode_count = total_counts.entry(barcode.to_string()).or_insert(0);
                *barcode_count += *count;
            }
        }
        SampleCounts(total_counts)
    }

    pub fn barcode_count_vec<'a, I: Iterator<Item = &'a SampleCounts>>(
        counts_iter: I,
        barcode: &'a str,
    ) -> Vec<usize> {
        counts_iter
            .map(|SampleCounts(count_table)| count_table.get(barcode).map_or(0, |ct| *ct))
            .collect()
    }
}

impl IntoIterator for SampleCounts {
    type Item = (String, usize);
    type IntoIter = ::std::collections::hash_map::IntoIter<String, usize>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}
