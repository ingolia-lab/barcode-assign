extern crate clap;
#[macro_use]
extern crate failure;
extern crate rust_htslib;

use std::str::FromStr;
use std::path::{Path,PathBuf};

use clap::{Arg, App};

use rust_htslib::bam;
use rust_htslib::bam::Read;

fn main() {
    let matches = App::new("bc-bam-chunk")
        .version("1.0")
        .author("Nick Ingolia <ingolia@berkeley.edu>")
        .about("Split BAM file into chunks of a specified size")
        .arg(Arg::with_name("bam_input")
             .short("i")
             .long("input")
             .value_name("INPUT.BAM")
             .help("Large BAM file")
             .takes_value(true)
             .required(true))
        .arg(Arg::with_name("outbase")
             .short("o")
             .long("outbase")
             .value_name("OUTBASE.BAM")
             .help("Output filename base (sequences in OUTBASE_001.bam, OUTBASE_002.bam, &c.)")
             .takes_value(true)
             .required(true))
        .arg(Arg::with_name("nreads")
             .short("n")
             .long("nreads")
             .value_name("#READS")
             .help("Number of reads per chunk")
             .takes_value(true)
             .required(true))
        .get_matches();

    let input = matches.value_of("bam_input").unwrap();
    let outbase = matches.value_of("outbase").unwrap();
    let nreads_str = matches.value_of("nreads").unwrap();
    
    match bc_bam_chunk(&input, &outbase, &nreads_str) {
        Ok(_) => (),
        Err(e) => panic!(e),
    }
}

fn bc_bam_chunk(input_filename: &str, outbase: &str, nreads_str: &str) -> Result<(), failure::Error> {
    let nreads = usize::from_str(nreads_str)?;

    let mut input = bam::Reader::from_path(input_filename)?;
    let header = bam::Header::from_template(input.header());

    let mut chunk_no = 0;
    let mut read_no = 0;
    
    let mut output = output_writer(outbase, &header, chunk_no)?;
    for rr in input.records() {
        let r = rr?;

        if read_no >= nreads {
            chunk_no += 1;
            read_no = 0;
            output = output_writer(outbase, &header, chunk_no)?;
        }

        output.write(&r)?;
        read_no += 1;
    }
    
    Ok(())
}

fn output_writer(outbase: &str, input_header: &bam::header::Header, chunk_no: usize) -> Result<bam::Writer, failure::Error> {
    let mut chunk_header = input_header.clone();
    chunk_header.push_record(bam::header::HeaderRecord::new(b"PG")
                             .push_tag(b"ID", &"bc-bam-chunk")
                             .push_tag(b"PN", &"bc-bam-chunk"));
    chunk_header.push_comment(format!("Chunk {:03}", chunk_no).as_bytes());
    Ok(bam::Writer::from_path(output_filename(outbase, chunk_no)?, &chunk_header)?)
}

fn output_filename(outbase: &str, chunk_no: usize) -> Result<PathBuf, failure::Error> {
    let out_path: &Path = outbase.as_ref();
    let stem = out_path.file_stem().ok_or_else(|| format_err!("Bad output base {:?}", outbase))?;
    let ext = out_path.extension().map_or_else(String::new, |ext| format!(".{}", ext.to_string_lossy()));
    Ok(out_path.with_file_name(format!("{}_{:03}{}", stem.to_string_lossy(), chunk_no, ext)))
}
