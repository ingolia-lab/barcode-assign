use std::boxed::Box;
use std::fs;
use std::io::Read as _;
use std::io::Write;
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use bio::alphabets::dna;
use bio::io::fastq;
use bio_types::strand::ReqStrand;
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
        let config_toml: ConfigTOML = toml::from_str(&contents).context("ConfigTOML")?;

        let mut lib_spec = LibSpec::new(&config_toml)?;

        println!("{:?}", config_toml);

        pacbio_extract(&mut lib_spec, &mut bam_in)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConfigTOML {
    output_base: String,
    fates: Option<String>,
    matching: Option<String>,
    max_errors: Option<u8>,
    inserts: Vec<InsertTOML>,
}

impl ConfigTOML {
    pub fn inserts_filename(&self) -> PathBuf {
        self.output_filename("-read-inserts-good.txt")
    }

    pub fn fates_filename(&self) -> PathBuf {
        self.fates.as_ref().map_or_else(
            || self.output_filename("-read-fates.txt"),
            |f| PathBuf::from(f),
        )
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
pub struct InsertTOML {
    name: String,
    before: String,
    after: String,
    min_len: Option<usize>,
    max_len: Option<usize>,
    max_errors: Option<u8>,
    reverse: Option<bool>,
    fastq: Option<String>,
}

pub struct LibSpec {
    fates_out: Box<dyn Write>,
    matching_out: Box<dyn Write>,
    inserts_out: Box<dyn Write>,
    insert_specs: Vec<InsertSpec>,
}

impl LibSpec {
    pub fn new(config: &ConfigTOML) -> Result<Self> {
        let inserts_filename = config.inserts_filename();
        let fates_filename = config.fates_filename();

        let mut insert_specs = Vec::new();
        for insert_config in config.inserts.iter() {
            insert_specs.push(InsertSpec::new(config, insert_config)?);
        }

        Ok(LibSpec {
            fates_out: Box::new(std::fs::File::create(fates_filename)?),
            matching_out: Box::new(std::io::sink()),
            inserts_out: Box::new(std::fs::File::create(inserts_filename)?),
            insert_specs: insert_specs,
        })
    }

    pub fn best_match<'a>(
        &mut self,
        query: &'a [u8],
        query_qual: &'a [u8],
    ) -> Vec<FlankMatchOut<'a>> {
        let mut insert_match_out = Vec::new();
        for insert_spec in self.insert_specs.iter_mut() {
            insert_match_out.push(insert_spec.best_match(query, query_qual));
        }
        insert_match_out
    }
}

pub struct InsertSpec {
    name: String,
    match_spec: FlankMatchSpec,
    min_len: Option<usize>,
    max_len: Option<usize>,
    reverse: bool,
    fastq_writer: fastq::Writer<Box<dyn Write>>,
}

impl InsertSpec {
    pub fn new(config: &ConfigTOML, insert_config: &InsertTOML) -> Result<Self> {
        let max_err = insert_config
            .max_errors
            .or(config.max_errors)
            .unwrap_or(Self::DEFAULT_MAX_ERRORS);

        let matcher = FlankMatchSpec::new(
            &Self::make_seq(&insert_config.before)?,
            &Self::make_seq(&insert_config.after)?,
            max_err,
        );

        let fastq_out: Box<dyn Write> = match &insert_config.fastq {
            Some(filename) => Box::new(std::fs::File::create(filename)?),
            None => Box::new(std::io::sink()),
        };
        let fastq_writer = fastq::Writer::new(fastq_out);

        Ok(InsertSpec {
            name: insert_config.name.clone(),
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

    pub fn best_match<'a>(&mut self, query: &'a [u8], query_qual: &'a [u8]) -> FlankMatchOut<'a> {
        self.match_spec.best_match(query, query_qual)
    }
}

pub fn pacbio_extract(spec: &mut LibSpec, bam_in: &mut bam::Reader) -> Result<()> {
    let mut rec = bam::Record::new();

    'bam: loop {
        match bam_in.read(&mut rec) {
            Some(Ok(())) => (),
            Some(Err(e)) => bail!(e),
            None => return Ok(()),
        }

        let read_id = String::from_utf8_lossy(extract_read_id(rec.qname()));

        let sequ_fwd = rec.seq().as_bytes();
        let qual_fwd = rec.qual();

        let sequ_rev = dna::revcomp(&sequ_fwd);
        let mut qual_rev = qual_fwd.to_vec();
        qual_rev.reverse();

        let lib_matches: Vec<(ReqStrand, Vec<FlankMatchOut>)> = vec![
            (ReqStrand::Forward, spec.best_match(&sequ_fwd, &qual_fwd)),
            (ReqStrand::Reverse, spec.best_match(&sequ_rev, &qual_rev)),
        ];

        for (ref strand, ref match_out) in lib_matches.iter() {
            write!(spec.matching_out, "{}\t{}", read_id, strand)?;
            for insert_match_out in match_out.iter() {
                write!(
                    spec.matching_out,
                    "\t{}\t{}",
                    insert_match_out.before_match_desc(),
                    insert_match_out.after_match_desc()
                )?;
            }
            write!(spec.matching_out, "\n")?;
        }

        let good_matches: Vec<(ReqStrand, Vec<FlankMatch>)> = lib_matches
            .iter()
            .filter_map(|(strand, lib_match_out)| {
                lib_match_out
                    .iter()
                    .map(FlankMatchOut::flank_match)
                    .collect::<Option<Vec<_>>>()
                    .map(|lm| (*strand, lm))
            })
            .collect();

        if good_matches.len() == 0 {
            write!(spec.fates_out, "{}\tNone\n", read_id)?;
        } else if good_matches.len() == 1 {
            let (ref strand, ref lib_match) = good_matches[0];

            for (insert_spec, insert_match) in spec.insert_specs.iter().zip(lib_match.iter()) {
                let insert_len = insert_match.insert_seq().len();
                if insert_spec
                    .min_len
                    .map_or(false, |min_len| insert_len < min_len)
                {
                    write!(spec.fates_out, "{}\t{}-short\n", read_id, insert_spec.name)?;
                    continue 'bam;
                } else if insert_spec
                    .max_len
                    .map_or(false, |max_len| insert_len > max_len)
                {
                    write!(spec.fates_out, "{}\t{}-long\n", read_id, insert_spec.name)?;
                    continue 'bam;
                }
            }

            write!(spec.fates_out, "{}\t{}\n", read_id, strand)?;

            write!(spec.inserts_out, "{}\t{}", read_id, strand)?;

            for (insert_spec, insert_match) in spec.insert_specs.iter_mut().zip(lib_match.iter()) {
                let insert_seq = if insert_spec.reverse {
                    dna::revcomp(insert_match.insert_seq())
                } else {
                    insert_match.insert_seq().to_vec()
                };
                let insert_qual = if insert_spec.reverse {
                    insert_match
                        .insert_qual()
                        .iter()
                        .map(|q| q + 33)
                        .collect::<Vec<_>>()
                } else {
                    insert_match
                        .insert_qual()
                        .iter()
                        .map(|q| q + 33)
                        .rev()
                        .collect::<Vec<_>>()
                };
                insert_spec.fastq_writer.write(
                    &format!("{}/0_{}", read_id, insert_seq.len()),
                    None,
                    &insert_seq,
                    &insert_qual,
                )?;

                write!(
                    spec.inserts_out,
                    "\t{}\t{}",
                    insert_match.insert_start(),
                    String::from_utf8_lossy(&insert_seq)
                )?;
            }

            write!(spec.inserts_out, "\n")?;
        } else {
            write!(spec.fates_out, "{}\tMulti\n", read_id)?;
        }
    }
}

pub fn extract_read_id(qname: &[u8]) -> &[u8] {
    let mut slash_iter = qname.rsplitn(2, |&ch| ch == b'/');
    let split_final = slash_iter.next().unwrap();
    slash_iter.next().unwrap_or(split_final)
}
