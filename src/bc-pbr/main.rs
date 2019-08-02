extern crate barcode_assign;
extern crate bio;

use std::cmp::*;
use std::io;

use bio::alphabets::dna;
use bio::io::fastq;

use barcode_assign::pacbio::*;

fn main() {
    let facs_conf = LibConf::new(MatchConf::new(b"GCTCGGAGATGTGTATAAGAGACAG",
                                                b"CTGTCTCTTATACACATCTGACGCC",
                                                4),
                                 MatchConf::new(b"CTATAGCACGACGCTCTTCCGATCT",
                                                b"GATCCTGTAGCCCTAGACTTGATAG",
                                                4),
                                 false);

    let syntf_conf = LibConf::new(MatchConf::new(b"CTCGGAGATGTGTATAAGAGACAG",
                                                 b"CTGTCTCTTATACACATCTGACGC",
                                                 4),
                                  MatchConf::new(b"GCGCTCTGTTGATAACTCCGGATC",
                                                 b"AGATCGGAAGAGCGTCGTGCTATA",
                                                 4),
                                  true);
    
    let fastq_filename = "/mnt/ingolialab/ingolia/Prog/pool-analysis/PacBio190731/samples.fastq";
    let fqin = fastq::Reader::from_file(&fastq_filename).unwrap();
        
    for recres in fqin.records() {
        let rec = recres.unwrap();

        let sequ_fwd = rec.seq();
        let sequ_rev = dna::revcomp(sequ_fwd);
        
        if let Some(facs_match) = facs_conf.best_match(&sequ_fwd) {
            println!("{}\tFACS\t+\t{}", rec.id(), format_match(&facs_match));
        } else if let Some(facs_match) = facs_conf.best_match(&sequ_rev) {
                println!("{}\tFACS\t-\t{}", rec.id(), format_match(&facs_match));
        } else if let Some(syntf_match) = syntf_conf.best_match(&sequ_fwd) {
            println!("{}\tSYNTF\t+\t{}", rec.id(), format_match(&syntf_match));
        } else if let Some(syntf_match) = syntf_conf.best_match(&sequ_rev) {
                println!("{}\tSYNTF\t-\t{}", rec.id(), format_match(&syntf_match));
        } else {
            println!("{}\tN/A", rec.id());
        }
    }
}

fn format_match<'a>(res: &'a LibMatchResult<'a>) -> String {
    let frag = res.frag_match.insert;
    
    let frag_start = &frag[0..min(30,frag.len())];
    let frag_end = &frag[(max(30, frag.len())-30)..frag.len()];

    let earlier_after = min(res.frag_match.after_match.0, res.barcode_match.after_match.0) as isize;
    let later_before = max(res.frag_match.before_match.1, res.barcode_match.before_match.1) as isize;
    
    format!("{}\t{}\t{}\t{}...{}",
            String::from_utf8_lossy(res.barcode_match.insert),
            later_before - earlier_after,
            frag.len(),
            String::from_utf8_lossy(frag_start),
            String::from_utf8_lossy(frag_end),
    )
}
