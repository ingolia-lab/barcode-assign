extern crate bio;
#[macro_use]
extern crate clap;
#[macro_use]
extern crate error_chain;
extern crate rust_htslib;
extern crate barcode_assign;

use std::fs;

use std::collections::HashMap;
use std::io::{Write};
use std::path::{Path,PathBuf};
use std::str;

use clap::{Arg, App};
use rust_htslib::bam;
use rust_htslib::prelude::*;

use barcode_assign::barcode_group::*;

mod assign;
mod depth;
mod purity;

use assign::*;
use depth::*;
use purity::*;

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
    min_purity: f64
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

        // ZZZ
        min_purity: 0.89
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
    let mut fates: HashMap<Fate,usize> = HashMap::new();

    let mut good_assign = fs::File::create(config.outfile("barcode-grna-good.txt"))?;
    write!(good_assign, "barcode\tguide\n")?;

    let mut depth_out = fs::File::create(config.outfile("barcode-depth.txt"))?;
    write!(depth_out, "{}\n", Depth::header())?;

    let mut purity_out = fs::File::create(config.outfile("barcode-purity.txt"))?;
    write!(purity_out, "{}\n", Purity::header())?;
    
    let mut fidelity_out = fs::File::create(config.outfile("barcode-fidelity.txt"))?;
    write!(fidelity_out, "barcode\ttid\tpos\tcigar\tmd\n")?;
    
    let mut bam_reader = bam::Reader::from_path(&config.bowtie_bam)?;
    let header = bam::Header::from_template(bam_reader.header());
    let header_view = bam::HeaderView::from_header(&header);

    let targets_result: std::result::Result<Vec<&str>,std::str::Utf8Error> = header_view.target_names().iter()
        .map(|name_u8| str::from_utf8(name_u8))
        .collect();
    let targets = targets_result?;
    
    let barcode_groups = BarcodeGroups::new(&mut bam_reader)?;

    for barcode_group in barcode_groups {
        let (bc, qall) = barcode_group?;
        let bc_str = str::from_utf8(bc.as_slice()).unwrap_or("???");

        let (qvec, depth) = filter_by_quality(qall, config.min_qual);
        write!(depth_out, "{}\n", depth.line(bc_str))?;
        
        let fate = if depth.n_good() < config.min_reads {
            Fate::NoDepth
        } else {
            let purity = Purity::new(qvec.iter())?;
            write!(purity_out, "{}\n", purity.line(bc_str))?;

            if purity.purity() < config.min_purity {
                Fate::NoPurity
            } else {
                if let ReadAssign::Match(assign) = purity.primary_assign() {
                    write!(fidelity_out, "{}\t{}\t{}\t{:?}\t{}\n",
                           bc_str, assign.target(&targets), assign.pos(),
                           assign.cigar(), str::from_utf8(assign.md())?)?;

                    if (assign.pos() == config.align_start as i32)
                        && assign.is_cigar_perfect(config.align_len as u32)
                        && assign.is_md_perfect(config.align_len as u32)
                    {
                        write!(good_assign, "{}\t{}\n",
                               bc_str, assign.target(&targets))?;
                        Fate::Good
                    } else {
                        Fate::NoFidelity
                    }
                } else {
                    Fate::NoMatch
                }
            }
        };

        *fates.entry(fate).or_insert(0) += 1;
    }

    let mut fates_out = fs::File::create(config.outfile("barcode-fates.txt"))?;
    write!(fates_out, "Good\t{}\n", fates.get(&Fate::Good).unwrap_or(&0))?;
    write!(fates_out, "Bad Match\t{}\n", fates.get(&Fate::NoFidelity).unwrap_or(&0))?;
    write!(fates_out, "No Match\t{}\n", fates.get(&Fate::NoMatch).unwrap_or(&0))?;
    write!(fates_out, "Mixed\t{}\n", fates.get(&Fate::NoPurity).unwrap_or(&0))?;
    write!(fates_out, "Too Few\t{}\n", fates.get(&Fate::NoDepth).unwrap_or(&0))?;
    
    Ok( () )
}

#[derive(Debug,Clone,Hash,PartialEq,Eq)]
enum Fate {
    Good,
    NoFidelity,
    NoMatch,
    NoPurity,
    NoDepth
}
