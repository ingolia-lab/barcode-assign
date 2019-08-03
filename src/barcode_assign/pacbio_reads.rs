use std::cmp::*;
use std::io::Write;
use std::path::{Path,PathBuf};

use bio::alphabets::dna;
use bio::io::fasta;
use bio::io::fastq;
use failure;

use flank_match::*;

pub fn pacbio_reads<P: AsRef<Path>, Q: AsRef<Path>>(specs: &[&LibSpec], fastq_filename: P, out_base: Q) -> Result<(), failure::Error>
{
    let fastq_in = fastq::Reader::from_file(fastq_filename)?;

    let mut frag_out = fasta::Writer::to_file(output_filename(&out_base, "-frags.fasta"))?;
    let mut good_insert_out = std::fs::File::create(output_filename(&out_base, "-read-inserts-good.txt"))?;
    let mut all_match_out = std::fs::File::create(output_filename(&out_base, "-read-matching-all.txt"))?;
    let mut fates_out = std::fs::File::create(output_filename(&out_base, "-read-fates.txt"))?;

    for recres in fastq_in.records() {
        let rec = recres?;

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
            write!(all_match_out, "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
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
            write!(fates_out, "{}\tNone\n", read_id)?;
        } else if good_matches.len() == 1 {
            let (ref name, ref strand, ref lib_match) = good_matches[0];
            write!(fates_out, "{}\t{}\t{}\n", read_id, name, strand)?;
            let frag_seq = lib_match.frag_match().insert_seq();
            frag_out
                .write(&format!("{}/0_{}", read_id, frag_seq.len()), None, frag_seq)
                ?;
            write!(
                good_insert_out,
                "{}\t{}\t{}\t{}\n",
                read_id,
                name,
                strand,
                format_match(&lib_match)
            )
                ?;
        } else {
            write!(fates_out, "{}\tMulti\n", read_id)?;
        }
    }

    Ok(())
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

pub fn output_filename<Q: AsRef<Path>>(out_base: Q, name: &str) -> PathBuf {
    let mut namebase = out_base.as_ref().file_name().map_or(std::ffi::OsString::new(), std::ffi::OsStr::to_os_string);
    namebase.push(name);
    out_base.as_ref().with_file_name(namebase)
}
