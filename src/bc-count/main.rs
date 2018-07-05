extern crate bio;
extern crate clap;

#[macro_use]
extern crate error_chain;

use std::fs::File;
use std::io::{self,Read,Write};
use std::collections::HashMap;

use bio::io::fastq;
use clap::{Arg, App};

#[derive(Debug)]
struct Config {
    barcode_fastq: String,
    out_barcodes: String,
    freq_filename: Option<String>,
}

mod errors {
    error_chain!{
        foreign_links {
            IO(::std::io::Error);
            FromUtf8(::std::string::FromUtf8Error);
            Utf8(::std::str::Utf8Error);
        }
    }
}

use errors::*;

fn main() {
    let matches = App::new("bc-count")
        .version("1.0")
        .author("Nick Ingolia <ingolia@berkeley.edu>")
        .about("Count barcode sequences")
        .arg(Arg::with_name("fastq")
             .short("f")
             .long("fastq")
             .value_name("BARCODE-FQ")
             .help("FastQ file of barcode sequences")
             .takes_value(true)
             .required(true))
        .arg(Arg::with_name("output")
             .short("o")
             .long("output")
             .value_name("OUTPUT-TXT")
             .help("Tab-delimited text file of barcode counts")
             .takes_value(true)
             .required(true))
        .get_matches();

    let config = Config {
        barcode_fastq: matches.value_of("fastq").unwrap().to_string(),
        out_barcodes: matches.value_of("output").unwrap().to_string(),
        freq_filename: None,
    };

    match run(config) {
        Ok(_) => (),
        Err(e) => panic!(e),
    }
}

fn run(config: Config) -> Result<()> {
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
        let mut freq_writer = File::create(freq_filename)?;
        write_freq_table(freq_writer, &barcode_counts)?;
    }
    
    Ok(())
}

fn write_barcode_table<W>(barcode_out: W, barcode_counts: &HashMap<String, usize>) -> Result<()>
    where W: std::io::Write
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

fn write_freq_table<W>(freq_out: W, barcode_counts: &HashMap<String, usize>) -> Result<()>
    where W: std::io::Write
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

