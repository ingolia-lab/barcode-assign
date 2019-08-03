extern crate barcode_assign;
extern crate bio;

use std::cmp::*;
use std::io::Write;

use bio::alphabets::dna;
use bio::io::fasta;
use bio::io::fastq;

use barcode_assign::pacbio::*;

fn main() {
    let facs_spec = LibSpec::new(
        FlankMatchSpec::new(
            b"GCTCGGAGATGTGTATAAGAGACAG",
            b"CTGTCTCTTATACACATCTGACGCC",
            3,
        ),
        FlankMatchSpec::new(
            b"CTATAGCACGACGCTCTTCCGATCT",
            b"GATCCTGTAGCCCTAGACTTGATAG",
            3,
        ),
        false,
    );

    let syntf_spec = LibSpec::new(
        FlankMatchSpec::new(b"CTCGGAGATGTGTATAAGAGACAG", b"CTGTCTCTTATACACATCTGACGC", 3),
        FlankMatchSpec::new(b"GCGCTCTGTTGATAACTCCGGATC", b"AGATCGGAAGAGCGTCGTGCTATA", 3),
        true,
    );

    let specs = vec![(&"FACS", &facs_spec), (&"SynTF", &syntf_spec)];

    let fastq_filename = "/mnt/ingolialab/ingolia/Prog/pool-analysis/PacBio190731/samples.fastq";
    let fqin = fastq::Reader::from_file(&fastq_filename).unwrap();

    let mut frag_out = fasta::Writer::to_file("frags.fasta").unwrap();
    let mut good_insert_out = std::fs::File::create("read-inserts-good.txt").unwrap();
    let mut all_match_out = std::fs::File::create("read-matching-all.txt").unwrap();
    let mut fates_out = std::fs::File::create("read-fates.txt").unwrap();

    for recres in fqin.records() {
        let rec = recres.unwrap();

        let read_id = &rec.id()[0..rec.id().rfind("/").unwrap_or(rec.id().len())];

        let sequ_fwd = rec.seq();
        let sequ_rev = dna::revcomp(sequ_fwd);

        let lib_matches: Vec<(String, String, LibMatchOut)> = specs
            .iter()
            .flat_map(|(&name, ref spec)| {
                vec![
                    (
                        name.to_string(),
                        "Fwd".to_string(),
                        spec.best_match(&sequ_fwd),
                    ),
                    (
                        name.to_string(),
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
                   match_out.frag_match().after_match_desc()).unwrap();                   
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
            write!(fates_out, "{}\tNone\n", read_id).unwrap();
        } else if good_matches.len() == 1 {
            let (ref name, ref strand, ref lib_match) = good_matches[0];
            write!(fates_out, "{}\t{}\t{}\n", read_id, name, strand).unwrap();
            let frag_seq = lib_match.frag_match().insert_seq();
            frag_out
                .write(&format!("{}/0_{}", read_id, frag_seq.len()), None, frag_seq)
                .unwrap();
            write!(
                good_insert_out,
                "{}\t{}\t{}\t{}\n",
                read_id,
                name,
                strand,
                format_match(&lib_match)
            )
            .unwrap();
        } else {
            write!(fates_out, "{}\tMulti\n", read_id).unwrap();
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
