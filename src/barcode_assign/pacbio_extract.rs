use std::boxed::Box;
use std::cmp::*;
use std::fs;
use std::io::{BufRead, BufReader, Write};
use std::io::Read as _;
use std::path::{Path, PathBuf};
use std::str::FromStr;

use anyhow::{Context, Result, bail};
use bio::alphabets::dna;
use bio::io::fasta;
use bio::io::fastq;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use serde::{Deserialize, Serialize};
use toml;

use flank_match::*;

#[derive(Debug)]
pub struct CLI {
    pub input_bam: String,
    pub config_file: String,
}
        
impl CLI {
    pub fn run(&self) -> Result<()> {
        let mut bam_in = if self.input_bam == "-" {
            bam::Reader::from_stdin()?
        } else {
            bam::Reader::from_path(&self.input_bam)?
        };

        let mut file = fs::File::open(&self.config_file)
            .with_context(|| format!("opening config file {:?}", self.config_file))?;
        let mut contents = String::new();
        file.read_to_string(&mut contents)
            .with_context(|| format!("reading config file {:?}", self.config_file))?;
        let config_toml: ConfigTOML = toml::from_str(&contents)
            .context("ConfigTOML")?;

        let lib_spec = LibSpec::new(&config_toml)?;
        
        println!("{:?}", config_toml);
        Ok(())
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConfigTOML
{
    output_base: String,
    fates: Option<String>,
    matching: Option<String>,
    max_errors: Option<u8>,
    inserts: Vec<InsertTOML>,
}

impl ConfigTOML {
    pub fn fates_filename(&self) -> PathBuf {
        self.fates.as_ref().map_or_else(
            || self.output_filename("-read-fates.txt"),
            |f| PathBuf::from(f))
    }
    
    pub fn output_filename(&self, name: &str) -> PathBuf {
        let base_ref: &Path = self.output_base.as_ref();
        let mut namebase = base_ref
            .file_name()
            .map_or(std::ffi::OsString::new(), std::ffi::OsStr::to_os_string);
        namebase.push(name);
        base_ref.with_file_name(namebase)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct InsertTOML
{
    name: String,
    before: String,
    after: String,
    min_len: Option<usize>,
    max_len: Option<usize>,
    max_errors: Option<u8>,
    reverse: Option<bool>,
    fastq: Option<String>,
}

pub struct LibSpec
{
    fates_out: Box<dyn Write>,
    matching_out: Box<dyn Write>,
    insert_specs: Vec<InsertSpec>
}

impl LibSpec {
    pub fn new(config: &ConfigTOML) -> Result<Self> {
        let fates_filename = config.fates_filename();

        let mut insert_specs = Vec::new();
        for insert_config in config.inserts.iter() {
            insert_specs.push(InsertSpec::new(config, insert_config)?);
        }
        
        Ok( LibSpec { fates_out: Box::new(std::fs::File::create(fates_filename)?),
                      matching_out: Box::new(std::io::sink()),
                      insert_specs: insert_specs,
        })

    }
    
}

pub struct InsertSpec
{
    name: String,
    match_spec: FlankMatchSpec,
    min_len: Option<usize>,
    max_len: Option<usize>,
    reverse: bool,
    fastq_writer: fastq::Writer<Box<dyn Write>>,
}

impl InsertSpec {
    pub fn new(config: &ConfigTOML, insert_config: &InsertTOML) -> Result<Self> {
        let max_err = insert_config.max_errors.or(config.max_errors).unwrap_or(Self::DEFAULT_MAX_ERRORS);
            
        let matcher = FlankMatchSpec::new(&Self::make_seq(&insert_config.before)?,
                                          &Self::make_seq(&insert_config.after)?,
                                          max_err);

        let fastq_out: Box<dyn Write> = match &insert_config.fastq {
            Some(filename) => Box::new(std::fs::File::create(filename)?),
            None => Box::new(std::io::sink()),
        };
        let fastq_writer = fastq::Writer::new(fastq_out);
        
        Ok(InsertSpec { name: insert_config.name.clone(),
                        match_spec: matcher,
                        min_len: insert_config.min_len,
                        max_len: insert_config.max_len,
                        reverse: insert_config.reverse.unwrap_or(false),
                        fastq_writer: fastq_writer,
        })
    }
    
    fn make_seq(raw: &str) -> Result<Vec<u8>> {
        let uc = raw.as_bytes().to_ascii_uppercase();
        if !uc.iter().all(|ch| b"ACGTN".contains(&ch)) {
            bail!("Bad sequence string {:?}", raw);
        }
        Ok(uc.to_vec())
    }

    const DEFAULT_MAX_ERRORS: u8 = 3;
}
