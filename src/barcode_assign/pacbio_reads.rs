use std::boxed::Box;
use std::cmp::*;
use std::fs;
use std::io::{BufRead,BufReader,Read,Write};
use std::path::{Path,PathBuf};
use std::str::FromStr;

use bio::alphabets::dna;
use bio::io::fasta;
use bio::io::fastq;
use bio::io::fastq::FastqRead;
use failure;

use flank_match::*;

#[derive(Debug)]
pub struct CLI {
    pub input_fastq: String,
    pub input_specs_file: String,
    pub output_base: String,
    pub output_file_frags: Option<String>,
    pub output_file_inserts: Option<String>,
    pub output_file_fates: Option<String>,
    pub output_file_matching: Option<String>,
    pub output_matching: bool,
    pub max_errors_str: String,
}

impl CLI {
    pub fn output_filename(&self, name: &str) -> PathBuf {
        let base_ref: &Path = self.output_base.as_ref();
        let mut namebase = base_ref.file_name().map_or(std::ffi::OsString::new(), std::ffi::OsStr::to_os_string);
        namebase.push(name);
        base_ref.with_file_name(namebase)
    }

    pub fn parse_lib_spec(line: &str, max_errors: u8) -> Result<LibSpec, failure::Error> {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != 6 {
            bail!("Malformed specification, expecting 6 fields: {:?}");
        }

        fn make_seq(raw: &str) -> Result<Vec<u8>, failure::Error> {
            let uc = raw.as_bytes().to_ascii_uppercase();
            if !uc.iter().all(|ch| b"ACGTN".contains(&ch)) {
                bail!("Bad sequence string {:?}", raw);
            }
            Ok(uc.to_vec())
        }

        let frag_matcher = FlankMatchSpec::new(&make_seq(fields[1])?, 
                                               &make_seq(fields[2])?,
                                               max_errors);
        let barcode_matcher = FlankMatchSpec::new(&make_seq(fields[3])?, 
                                                  &make_seq(fields[4])?,
                                                  max_errors);
        let barcode_rev = bool::from_str(fields[5])?;
        Ok(LibSpec::new(fields[0], frag_matcher, barcode_matcher, barcode_rev))
    }

    pub fn max_errors(&self) -> Result<u8, failure::Error> {
        Ok(u8::from_str(&self.max_errors_str)?)
        // ZZZ annotate error
    }
    
    pub fn read_lib_specs(&self) -> Result<Vec<LibSpec>, failure::Error> {
        let max_errors = self.max_errors()?;
        let mut specs = Vec::new();
        for line_res in BufReader::new(fs::File::open(&self.input_specs_file)?).lines() {
            let line = line_res?;
            if line.len() > 0 && !line.starts_with("#") {
                specs.push(Self::parse_lib_spec(&line, max_errors)?);
            }
        }
        Ok(specs)
    }

    pub fn outputs(&self) -> Result<Outputs, failure::Error> {
        let fasta_writer: Box<dyn Write> = Box::new(std::fs::File::create(self.output_filename("-frags.fasta"))?);
        
        let frag_out = fasta::Writer::new(fasta_writer);
        let good_insert_out = std::fs::File::create(self.output_filename("-read-inserts-good.txt"))?;
        let fates_out = std::fs::File::create(self.output_filename("-read-fates.txt"))?;

        let all_match_out: Box<dyn Write> = if self.output_matching {
            Box::new(std::fs::File::create(self.output_filename("-read-matching-all.txt"))?)
        } else {
            Box::new(std::io::sink())
        };
        
        Ok(Outputs { frags: frag_out,
                     inserts: Box::new(good_insert_out),
                     fates: Box::new(fates_out),
                     matching: all_match_out })
    }
    
    pub fn run(&self) -> Result<(), failure::Error> {
        let fastq_reader: Box<dyn Read> = if self.input_fastq == "-" {
            Box::new(std::io::stdin())
        } else {
            Box::new(std::fs::File::open(&self.input_fastq)?)
        };
      
        {
            let mut fastq_in = fastq::Reader::new(fastq_reader);
            let specs = self.read_lib_specs()?;
            let mut outputs = self.outputs()?;
            
            pacbio_reads(specs.as_slice(), &mut fastq_in, &mut outputs)
        }
    }
}

pub struct Outputs {
    frags: fasta::Writer<Box<dyn Write>>,
    inserts: Box<dyn Write>,
    fates: Box<dyn Write>,
    matching: Box<dyn Write>
}

impl Outputs {
    pub fn frags(&mut self) -> &mut fasta::Writer<Box<dyn Write>> { &mut self.frags }
    pub fn inserts(&mut self) -> &mut Write { self.inserts.as_mut() }
    pub fn fates(&mut self) -> &mut Write { self.fates.as_mut() }
    pub fn matching(&mut self) -> &mut Write { self.matching.as_mut() }
}

pub fn pacbio_reads<R: Read>(specs: &[LibSpec], fastq_in: &mut fastq::Reader<R>, outputs: &mut Outputs) -> Result<(), failure::Error>
{
    let mut rec = fastq::Record::new();
    
    loop {
        fastq_in.read(&mut rec)?;

        if rec.is_empty() {
            return Ok(());
        }

        let read_id = &rec.id()[0..rec.id().rfind("/").unwrap_or(rec.id().len())];

        let sequ_fwd = rec.seq();
        let sequ_rev = dna::revcomp(sequ_fwd);
        let lib_matches: Vec<(String, String, LibMatchOut)> = specs
            .iter()
            .flat_map(|ref spec| {
                vec![
                    (
                        spec.name().to_string(),
                        "Fwd".to_string(),
                        spec.best_match(&sequ_fwd),
                    ),
                    (
                        spec.name().to_string(),
                        "Rev".to_string(),
                        spec.best_match(&sequ_rev),
                    ),
                ]
            })
            .collect();

        for (ref lib, ref strand, ref match_out) in lib_matches.iter() {
            write!(outputs.matching(), "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                   read_id, lib, strand,
                   match_out.barcode_match().before_match_desc(),
                   match_out.barcode_match().after_match_desc(),
                   match_out.frag_match().before_match_desc(),
                   match_out.frag_match().after_match_desc())?;                   
        }
        
        let good_matches: Vec<(String, String, LibMatch)> = lib_matches
            .iter()
            .filter_map(|(name, strand, lib_match_out)| {
                lib_match_out
                    .lib_match()
                    .map(|lib_match| (name.to_string(), strand.to_string(), lib_match))
            })
            .collect();

        if good_matches.len() == 0 {
            write!(outputs.fates(), "{}\tNone\n", read_id)?;
        } else if good_matches.len() == 1 {
            let (ref name, ref strand, ref lib_match) = good_matches[0];
            write!(outputs.fates(), "{}\t{}\t{}\n", read_id, name, strand)?;
            let frag_seq = lib_match.frag_match().insert_seq();
            outputs.frags()
                .write(&format!("{}/0_{}", read_id, frag_seq.len()), None, frag_seq)
                ?;
            write!(
                outputs.inserts(),
                "{}\t{}\t{}\t{}\n",
                read_id,
                name,
                strand,
                format_match(&lib_match)
            )
                ?;
        } else {
            write!(outputs.fates(), "{}\tMulti\n", read_id)?;
        }
    }
}

fn format_match<'a>(res: &'a LibMatch<'a>) -> String {
    let frag = res.frag_match().insert_seq();

    let frag_start = &frag[0..min(30, frag.len())];
    let frag_end = &frag[(max(30, frag.len()) - 30)..frag.len()];

    format!(
        "{}\t{}\t{}\t{}...{}",
        res.barcode_actual(),
        res.between(),
        frag.len(),
        String::from_utf8_lossy(frag_start),
        String::from_utf8_lossy(frag_end),
    )
}
