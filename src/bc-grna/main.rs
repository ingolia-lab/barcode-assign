extern crate bio;
#[macro_use]
extern crate clap;
#[macro_use]
extern crate error_chain;
extern crate rust_htslib;
extern crate barcode_assign;

use std::fs;

use std::io::{Write};
use std::path::{Path,PathBuf};
use std::rc::Rc;
use std::str;

use bio::io::fasta;
use clap::{Arg, App};
use rust_htslib::bam;
use rust_htslib::prelude::*;

use barcode_assign::barcode_group::*;

mod errors {
    error_chain!{
        links {
            BarcodeGroup(::barcode_assign::errors::Error, ::barcode_assign::errors::ErrorKind);
        }
        foreign_links {
            IO(::std::io::Error);
            FromUtf8(::std::string::FromUtf8Error);
            Utf8(::std::str::Utf8Error);
            BamRead(::rust_htslib::bam::ReadError);
            BamReaderPath(::rust_htslib::bam::ReaderPathError);
            BamWrite(::rust_htslib::bam::WriteError);
            BamWriterPath(::rust_htslib::bam::WriterPathError);
            Pileup(::rust_htslib::bam::pileup::PileupError);
        }
    }
}

use errors::*;

#[derive(Debug)]
struct Config {
    bowtie_bam: PathBuf,
    align_len: usize,
    align_start: usize,
    is_reverse: bool,
    out_base: PathBuf,
    min_reads: usize,
    min_qual: u8,
}

impl Config {
    pub fn outfile<T>(&self, file_suffix: T) -> PathBuf 
        where T: std::convert::AsRef<std::ffi::OsStr>
    {
        let mut filename = self.out_base.file_name().map_or_else(|| std::ffi::OsString::new(), |str| str.to_os_string());
        filename.push(file_suffix);
        self.out_base.with_file_name(filename)
    }
}

fn main() {
    let matches = App::new("bc-grna")
        .version("1.0")
        .author("Nick Ingolia <ingolia@berkeley.edu>")
        .about("Barcode / guide RNA assignment")
        .arg(Arg::with_name("bambyname")
             .short("b")
             .long("bam-by-name")
             .value_name("SORTED-BY-NAME-BAM")
             .help("BAM format alignment sorted by name")
             .takes_value(true)
             .required(true))
        .arg(Arg::with_name("outbase")
             .short("o")
             .long("outbase")
             .value_name("OUTBASE")
             .help("Output base filename")
             .takes_value(true)
             .required(true))
        .arg(Arg::with_name("alignlen")
             .long("align-len")
             .value_name("LEN")
             .help("Length of aligned sequence")
             .takes_value(true)
             .required(true))
        .arg(Arg::with_name("alignstart")
             .long("align-start")
             .value_name("START")
             .help("Alignment starting position (0-based)")
             .takes_value(true)
             .required(true))
        .arg(Arg::with_name("minreads")
             .long("min-reads")
             .value_name("MINREADS")
             .help("Minimum number of reads for good assignment")
             .takes_value(true)
             .default_value("3"))
        .arg(Arg::with_name("minqual")
             .long("min-qual")
             .value_name("MINQUAL")
             .help("Minimum (median) quality score to include read")
             .takes_value(true)
             .default_value("30"))
        .get_matches();
               
    let config = Config {
        bowtie_bam: PathBuf::from(matches.value_of("bambyname").unwrap()),
        align_len: value_t!(matches.value_of("alignlen"), usize).unwrap_or_else(|e| e.exit()),
        align_start: value_t!(matches.value_of("alignstart"), usize).unwrap_or_else(|e| e.exit()),

        is_reverse: false,
        
        out_base: PathBuf::from(matches.value_of("outbase").unwrap()),

        min_reads: value_t!(matches.value_of("minreads"), usize).unwrap_or_else(|e| e.exit()),

        min_qual: value_t!(matches.value_of("minqual"), u8).unwrap_or_else(|e| e.exit()),
    };

    if let Err(ref e) = run(&config) {
        println!("error: {}", e);
        
        for e in e.iter().skip(1) {
            println!("caused by: {}", e);
        }
        
        if let Some(backtrace) = e.backtrace() {
            println!("backtrace: {:?}", backtrace);
        }
        
        ::std::process::exit(1);
    }
}

fn run(config: &Config) -> Result<()> {
    barcode_to_grna(config)
}

fn barcode_to_grna(config: &Config) -> Result<()> {
    let mut good_assign = fs::File::create(config.outfile("barcode-grna-good.txt"))?;
    write!(good_assign, "barcode\tguide\n")?;

    let mut depth_out = fs::File::create(config.outfile("barcode-depth.txt"))?;
    write!(depth_out, "barcode\tnread\tnlowqual\n")?;

    let mut bam_reader = bam::Reader::from_path(&config.bowtie_bam)?;
    let header = bam::Header::from_template(bam_reader.header());
    let header_view = bam::HeaderView::from_header(&header);
    
    let barcode_groups = BarcodeGroups::new(&mut bam_reader)?;

    for barcode_group in barcode_groups {
        let (bc, qall) = barcode_group?;
        let bc_str = str::from_utf8(bc.as_slice()).unwrap_or("???");

        let nread_all = qall.len();
        
        let mut qvec = qall.into_iter()
            .filter(|r| (median_qual(r) >= config.min_qual))
            .collect::<Vec<bam::Record>>();
        qvec.sort_by_key(|r| (r.tid(), r.pos()));

        let nread = qvec.len();
        let nlowqual = nread_all - nread;
        
        write!(depth_out, "{}\t{}\t{}\n", bc_str, nread, nlowqual)?;
        
        if nread >= config.min_reads {
            
        }
    }

    Ok( () )
}

#[derive(Debug,Clone,Hash,PartialEq,Eq)]
enum ReadAssign {
    NoMatch,
    Match(AssignMatch)
}

impl ReadAssign {
    pub fn new(r: &bam::Record) -> Result<Self> {
        if r.tid() < 0 {
            Ok( ReadAssign::NoMatch )
        } else {
            let md_aux = r.aux(b"MD").ok_or_else(|| "No MD tag")?;
            let md: Vec<u8> = match md_aux {
                bam::record::Aux::String(md) => Ok( md.to_vec() ),
                _ => Err("MD tag not a string")
            }?;

            let cigar_view = r.cigar();
            let cigar: Vec<bam::record::Cigar> = cigar_view.iter().map(|&c| c.clone()).collect();
            
            let assign_match = AssignMatch { tid: r.tid() as u32, pos: r.pos(),
                                             is_reverse: r.is_reverse(),
                                             cigar: cigar, md: md };
            Ok( ReadAssign::Match(assign_match) )
        }
    }
}

#[derive(Debug,Clone,Hash,PartialEq,Eq)]
struct AssignMatch {
    tid: u32,
    pos: i32,
    is_reverse: bool,
    cigar: Vec<bam::record::Cigar>,
    md: Vec<u8>
}

fn median_qual(r: &bam::Record) -> u8 {
    let mut quals = r.qual().to_vec();
    quals.sort();
    quals.get(quals.len() / 2).map_or(0, |q| *q)
}

