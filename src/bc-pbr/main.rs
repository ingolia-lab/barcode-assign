extern crate barcode_assign;
extern crate bio;

use std::cmp::*;
use std::io;

use bio::alphabets::dna;
use bio::io::fasta;
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

    let mut frag_out = fasta::Writer::to_file("frags.fasta").unwrap();
    
    for recres in fqin.records() {
        let rec = recres.unwrap();

        let sequ_fwd = rec.seq();
        let sequ_rev = dna::revcomp(sequ_fwd);
        
        if let Some(facs_match) = facs_conf.best_match(&sequ_fwd) {
            println!("{}\tFACS\t+\t{}", rec.id(), format_match(&facs_match));
            frag_out.write(rec.id(), None, facs_match.frag_match().insert).unwrap();
        } else if let Some(facs_match) = facs_conf.best_match(&sequ_rev) {
            println!("{}\tFACS\t-\t{}", rec.id(), format_match(&facs_match));
            frag_out.write(rec.id(), None, facs_match.frag_match().insert).unwrap();
        } else if let Some(syntf_match) = syntf_conf.best_match(&sequ_fwd) {
            println!("{}\tSYNTF\t+\t{}", rec.id(), format_match(&syntf_match));
            frag_out.write(rec.id(), None, syntf_match.frag_match().insert).unwrap();
        } else if let Some(syntf_match) = syntf_conf.best_match(&sequ_rev) {
            println!("{}\tSYNTF\t-\t{}", rec.id(), format_match(&syntf_match));
            frag_out.write(rec.id(), None, syntf_match.frag_match().insert).unwrap();
        } else {
            println!("{}\tN/A", rec.id());
        }
    }
}

fn format_match<'a>(res: &'a LibMatch<'a>) -> String {
    let frag = res.frag_match().insert;
    
    let frag_start = &frag[0..min(30,frag.len())];
    let frag_end = &frag[(max(30, frag.len())-30)..frag.len()];

    let between = if res.frag_match().insert_end() < res.barcode_match().insert_start() {
        (res.frag_match().insert_end() as isize) - (res.frag_match().insert_start() as isize)
    } else {
        (res.frag_match().insert_start() as isize) - (res.barcode_match().insert_end() as isize)
    };

    format!("{}\t{}\t{}\t{}...{}",
            String::from_utf8_lossy(res.barcode_match().insert),
            between,
            frag.len(),
            String::from_utf8_lossy(frag_start),
            String::from_utf8_lossy(frag_end),
    )
}
