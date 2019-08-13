use std::collections::HashMap;
use std::fs::File;
use std::io::{self, Read, Write};

use bio::io::fastq;

#[derive(Debug)]
pub struct Config {
    pub barcode_fastq: String,
    pub out_barcodes: String,
    pub freq_filename: Option<String>,
}

pub fn bc_count(config: Config) -> Result<(), failure::Error> {
    let reader: Box<Read> = if config.barcode_fastq == "-" {
        Box::new(io::stdin())
    } else {
        Box::new(File::open(&config.barcode_fastq)?)
    };
    let barcode_reader = fastq::Reader::new(reader);

    let barcode_counts = count_barcodes(barcode_reader)?;

    let writer: Box<Write> = if config.out_barcodes == "-" {
        Box::new(io::stdout())
    } else {
        Box::new(File::create(&config.out_barcodes)?)
    };
    write_barcode_table(writer, &barcode_counts)?;

    if let Some(freq_filename) = config.freq_filename {
        let freq_writer = File::create(freq_filename)?;
        write_freq_table(freq_writer, &barcode_counts)?;
    }

    Ok(())
}

pub fn count_barcodes<R: Read>(barcode_reader: fastq::Reader<R>) -> Result<HashMap<Vec<u8>, usize>, failure::Error> {
    let mut barcode_counts = HashMap::new();
    
    for (_recno, rec_res) in barcode_reader.records().enumerate() {
        let rec = rec_res?;
        
        let barcode = rec.seq().to_vec();
        let barcode_count = barcode_counts.entry(barcode).or_insert(0);
        *barcode_count += 1;
    }

    Ok(barcode_counts)
}

fn write_barcode_table<W>(
    barcode_out: W,
    barcode_counts: &HashMap<Vec<u8>, usize>,
) -> Result<(), failure::Error>
where
    W: std::io::Write,
{
    let mut bcout = std::io::BufWriter::new(barcode_out);

    for (barcode, count) in barcode_counts.iter() {
        write!(bcout, "{}\t{}\n",
               String::from_utf8_lossy(barcode), count)?;
    }

    Ok(())
}

fn write_freq_table<W>(
    freq_out: W,
    barcode_counts: &HashMap<Vec<u8>, usize>,
) -> Result<(), failure::Error>
where
    W: std::io::Write,
{
    let mut fout = std::io::BufWriter::new(freq_out);

    let mut freq_counts = HashMap::new();

    for freq in barcode_counts.values() {
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
