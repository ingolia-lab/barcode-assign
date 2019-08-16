use std::io::Write;

use failure;

use counts::*;

pub struct CLI {
    pub inputs: Vec<String>,
    pub output: String,
    pub mintotal: Option<usize>,
    pub minsamples: Option<usize>,
    pub mininsample: Option<usize>,
    pub omitfile: Option<String>,
}

impl CLI {
    pub fn run(&self) -> Result<(), failure::Error> {
        let mut counts = Vec::new();

        for input in self.inputs.iter() {
            let input_counts = SampleCounts::from_file(&input)?;
            counts.push((input.to_string(), input_counts));
        }

        let total_counts = SampleCounts::total_counts(counts.iter().map(|(_input, counts)| counts));
        let mut barcodes: Vec<(Vec<u8>, usize)> = total_counts.into_iter().collect();
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
            let count_vec = SampleCounts::barcode_count_vec(
                counts.iter().map(|(_input, counts)| counts),
                barcode,
            );

            if self.is_omitted(&count_vec) {
                write!(omit_out, "{}\n", String::from_utf8_lossy(barcode))?;
            } else {
                write!(out, "{}", String::from_utf8_lossy(barcode))?;
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
            let nsamples: usize = count_vec.iter().filter(|&&ct| ct > 0).count();
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
