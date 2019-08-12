extern crate clap;
#[macro_use]
extern crate failure;
extern crate rust_htslib;

use std::path::{Path, PathBuf};
use std::str::FromStr;

use clap::{App, Arg};

use rust_htslib::bam;
use rust_htslib::bam::Read;

fn main() {
    let matches = App::new("bc-bam-chunk")
        .version("1.0")
        .author("Nick Ingolia <ingolia@berkeley.edu>")
        .about("Split BAM file into chunks of a specified size")
        .arg(
            Arg::with_name("bam_input")
                .short("i")
                .long("input")
                .value_name("INPUT.BAM")
                .help("Large BAM file")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("outbase")
                .short("o")
                .long("outbase")
                .value_name("OUTBASE.BAM")
                .help("Output filename base (sequences in OUTBASE_001.bam, OUTBASE_002.bam, &c.)")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("nreads")
                .short("n")
                .long("nreads")
                .value_name("#READS")
                .help("Number of reads per chunk")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("verbose")
                .short("v")
                .long("verbose")
                .help("Verbose report on chunking progress"),
        )
        .get_matches();

    let input = matches.value_of("bam_input").unwrap();
    let outbase = matches.value_of("outbase").unwrap();
    let nreads_str = matches.value_of("nreads").unwrap();
    let verbose = matches.occurrences_of("verbose") > 0;

    match bc_bam_chunk(&input, &outbase, &nreads_str, verbose) {
        Ok(_) => (),
        Err(e) => panic!(e),
    }
}

fn bc_bam_chunk(
    input_filename: &str,
    outbase: &str,
    nreads_str: &str,
    verbose: bool,
) -> Result<(), failure::Error> {
    let nreads = usize::from_str(nreads_str)?;

    let mut input = bam::Reader::from_path(input_filename)?;
    let header = bam::Header::from_template(input.header());

    let mut chunk_no = 0;
    let mut read_no = 0;
    let mut subread_no = 0;

    let mut output = output_writer(outbase, &header, chunk_no)?;
    let mut current_zmw = Vec::new();

    for rr in input.records() {
        let r = rr?;

        let r_zmw = get_zmw(r.qname())?;

        if r_zmw != current_zmw.as_slice() {
            read_no += 1;

            if read_no >= nreads {
                if verbose {
                    println!(
                        "Finished chunk {} with {} reads and {} subreads",
                        chunk_no, read_no, subread_no
                    );
                }

                chunk_no += 1;
                read_no = 0;
                subread_no = 0;
                output = output_writer(outbase, &header, chunk_no)?;

                if verbose {
                    print!(
                        "Starting chunk {} with {:?}...",
                        chunk_no,
                        String::from_utf8_lossy(r_zmw)
                    );
                }
            }

            current_zmw = r_zmw.to_vec();

            //            println!("zmw {}\tread {}\tchunk {}\tqname {}", String::from_utf8_lossy(current_zmw.as_slice()), read_no, chunk_no, String::from_utf8_lossy(r.qname()));
        }

        output.write(&r)?;
        subread_no += 1;
    }

    Ok(())
}

fn output_writer(
    outbase: &str,
    input_header: &bam::header::Header,
    chunk_no: usize,
) -> Result<bam::Writer, failure::Error> {
    let mut chunk_header = input_header.clone();
    chunk_header.push_record(
        bam::header::HeaderRecord::new(b"PG")
            .push_tag(b"ID", &"bc-bam-chunk")
            .push_tag(b"PN", &"bc-bam-chunk"),
    );
    chunk_header.push_comment(format!("Chunk {:03}", chunk_no).as_bytes());
    Ok(bam::Writer::from_path(
        output_filename(outbase, chunk_no)?,
        &chunk_header,
    )?)
}

fn output_filename(outbase: &str, chunk_no: usize) -> Result<PathBuf, failure::Error> {
    let out_path: &Path = outbase.as_ref();
    let stem = out_path
        .file_stem()
        .ok_or_else(|| format_err!("Bad output base {:?}", outbase))?;
    let ext = out_path
        .extension()
        .map_or_else(String::new, |ext| format!(".{}", ext.to_string_lossy()));
    Ok(out_path.with_file_name(format!("{}_{:03}{}", stem.to_string_lossy(), chunk_no, ext)))
}

fn get_zmw(qname: &[u8]) -> Result<&[u8], failure::Error> {
    let mut iter = qname.split(|ch| *ch == b'/');
    let _discard = iter
        .next()
        .ok_or_else(|| format_err!("Bad query name {:?}", qname))?;
    iter.next()
        .ok_or_else(|| format_err!("Bad query name {:?}", qname))
}
