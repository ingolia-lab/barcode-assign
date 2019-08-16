use std::collections::HashMap;
use std::io::BufRead;
use std::io::BufReader;
use std::io::Read;
use std::path::{Path,PathBuf};

use neighborhood::Neighborhood;

pub struct CLI {
    pub input: String,
    pub output_base: String,
}

impl CLI {
    pub fn output_filename(&self, name: &str) -> PathBuf {
        let base_ref: &Path = self.output_base.as_ref();
        let mut namebase = base_ref
            .file_name()
            .map_or(std::ffi::OsString::new(), std::ffi::OsStr::to_os_string);
        namebase.push(name);
        base_ref.with_file_name(namebase)
    }

    pub fn run(&self) -> Result<(), failure::Error> {
        let reader: Box<dyn Read> = if self.input == "-" {
            Box::new(std::io::stdin())
        } else {
            Box::new(std::fs::File::open(&self.input)?)
        };
        let mut barcode_reader = BufReader::new(reader);

        let barcode_counts = CLI::count_barcodes(&mut barcode_reader)?;

        let mut nbhds = Neighborhood::gather_neighborhoods(barcode_counts);
        for nbhd in nbhds.iter_mut() {
            nbhd.sort_by_counts();
        }

        let mut barcode_to_nbhd_out = std::fs::File::create(self.output_filename("-barcode-to-nbhd.txt"))?;
        let mut nbhd_count_out = std::fs::File::create(self.output_filename("-nbhd-count.txt"))?;
        let mut nbhds_out = std::fs::File::create(self.output_filename("-nbhds.txt"))?;
        
        for nbhd in nbhds.iter() {
            nbhd.write_total_counts(&mut nbhd_count_out)?;
            nbhd.write_barcode_counts(&mut barcode_to_nbhd_out)?;
            nbhd.write_nbhd_counts(&mut nbhds_out)?;
        }
        
        Ok(())
    }

    pub fn count_barcodes<R: BufRead>(barcode_reader: &mut R) -> Result<HashMap<Vec<u8>, usize>, failure::Error> {
        let mut barcode_counts = HashMap::new();
        
        for line_res in barcode_reader.lines() {
            let line = line_res?;
            
            let barcode = line.into_bytes();
            let barcode_count = barcode_counts.entry(barcode).or_insert(0);
            *barcode_count += 1;
        }
        
        Ok(barcode_counts)
    }
}

