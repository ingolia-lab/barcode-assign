use std::fs::File;
use std::io::{self, Read, Write};
use std::path::{Path,PathBuf};

use bio::io::fastq;

use neighborhood::Neighborhood;
use counts::SampleCounts;

#[derive(Debug)]
pub struct Config {
    pub barcode_fastq: String,
    pub out_barcodes: String,
    pub freq_filename: Option<String>,
    pub neighborhood: Option<String>,
}

pub fn bc_count(config: Config) -> Result<(), failure::Error> {
    let reader: Box<Read> = if config.barcode_fastq == "-" {
        Box::new(io::stdin())
    } else {
        Box::new(File::open(&config.barcode_fastq)?)
    };
    let barcode_reader = fastq::Reader::new(reader);

    let barcode_counts_res: Result<SampleCounts, std::io::Error>
        = barcode_reader.records().collect();
    let barcode_counts = barcode_counts_res?;
    
    let writer: Box<Write> = if config.out_barcodes == "-" {
        Box::new(io::stdout())
    } else {
        Box::new(File::create(&config.out_barcodes)?)
    };
    barcode_counts.write(writer)?;
    
    if let Some(freq_filename) = config.freq_filename {
        barcode_counts.write_freq_table(File::create(freq_filename)?)?;
    }

    if let Some(nbhd_filename) = config.neighborhood {
        neighborhood_counts(barcode_counts, &nbhd_filename)?;
    }
    
    Ok(())
}

fn neighborhood_counts(barcode_counts: SampleCounts, nbhd_filename: &str) -> Result<(), failure::Error>
{
    let mut barcode_to_nbhd_out = std::fs::File::create(output_filename(nbhd_filename, "-barcode-to-nbhd.txt"))?;
    let mut nbhd_count_out = std::fs::File::create(output_filename(nbhd_filename, "-nbhd-count.txt"))?;
    let mut nbhds_out = std::fs::File::create(output_filename(nbhd_filename, "-nbhds.txt"))?;

    let mut nbhds = Neighborhood::gather_neighborhoods(barcode_counts.count_map());
    for nbhd in nbhds.iter_mut() {
        nbhd.sort_by_counts();
    }
    
    for nbhd in nbhds.iter() {
        nbhd.write_total_counts(&mut nbhd_count_out)?;
        nbhd.write_barcode_counts(&mut barcode_to_nbhd_out)?;
        nbhd.write_nbhd_counts(&mut nbhds_out)?;
    }

    Ok(())
}

fn output_filename(output_base: &str, name: &str) -> PathBuf {
    let base_ref: &Path = output_base.as_ref();
    let mut namebase = base_ref
        .file_name()
        .map_or(std::ffi::OsString::new(), std::ffi::OsStr::to_os_string);
    namebase.push(name);
    base_ref.with_file_name(namebase)
}

