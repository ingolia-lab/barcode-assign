//use std::fs;
//use std::io::Write;
use std::path::{Path, PathBuf};
use std::str::FromStr;

use failure;

use rust_htslib::bam;
use rust_htslib::bam::Read as BamRead;

// mod record_class;
// mod record_group;
// mod stats;

// use bam_suppress_duplicates::record_class::*;
// use bam_suppress_duplicates::record_group::*;
// use bam_suppress_duplicates::stats::*;

pub struct CLI {
    pub bambyname: String,
    pub outbase: String,
    pub minreads: String,
    pub minqual: String,
    pub minpurity: String,
}

pub struct Config {
    input: bam::Reader,
    min_reads: usize,
    min_qual: u8,
    min_purity: f64,
}

const DEFAULT_NLIM: usize = 100; // ZZZ

impl Config {
    pub fn new(cli: &CLI) -> Result<Self, failure::Error> {
        let input = if cli.bambyname == "-" {
            bam::Reader::from_stdin()?
        } else {
            bam::Reader::from_path(Path::new(&cli.bambyname))?
        };

        Ok(Config {
            input: input,
            min_reads: usize::from_str(&cli.minreads)?,
            min_qual: u8::from_str(&cli.minqual)?,
            min_purity: f64::from_str(&cli.minpurity)?,
        })
    }
}

pub fn read_tag(r1: &bam::Record) -> Option<&[u8]> {
    if let Some(delim_pos) = r1.qname().iter().position(|&ch| ch == b'_') {
        Some(r1.qname().split_at(delim_pos).0)
    } else {
        None
    }
}

pub fn bc_frag(mut config: Config) -> Result<(), failure::Error> {
    Ok(())
}
