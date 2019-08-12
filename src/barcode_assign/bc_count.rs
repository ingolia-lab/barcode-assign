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

    let mut barcode_counts = HashMap::new();

    let records = barcode_reader.records();

    for result in records {
        let barcode_record = result?;

        let barcode = String::from_utf8(barcode_record.seq().to_vec())?;
        let barcode_count = barcode_counts.entry(barcode.to_string()).or_insert(0);
        *barcode_count += 1;
    }

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

fn write_barcode_table<W>(
    barcode_out: W,
    barcode_counts: &HashMap<String, usize>,
) -> Result<(), failure::Error>
where
    W: std::io::Write,
{
    let mut bcout = std::io::BufWriter::new(barcode_out);

    for (barcode, count) in barcode_counts.iter() {
        bcout.write(barcode.as_bytes())?;
        bcout.write("\t".as_bytes())?;
        bcout.write(count.to_string().as_bytes())?;
        bcout.write("\n".as_bytes())?;
    }

    Ok(())
}

fn write_freq_table<W>(
    freq_out: W,
    barcode_counts: &HashMap<String, usize>,
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
