extern crate barcode_assign;

use std::io::Write;

use barcode_assign::flank_match::*;
use barcode_assign::pacbio_reads::*;

fn main() {
    let facs_spec = LibSpec::new(
        "FACS",
        FlankMatchSpec::new(b"GCTCGGAGATGTGTATAAGAGACAG", b"CTGTCTCTTATACACATCTGACGCC", 3),
        FlankMatchSpec::new(b"CTATAGCACGACGCTCTTCCGATCT", b"GATCCTGTAGCCCTAGACTTGATAG", 3),
        false,
    );

    let syntf_spec = LibSpec::new(
        "SynTF",
        FlankMatchSpec::new(b"CTCGGAGATGTGTATAAGAGACAG", b"CTGTCTCTTATACACATCTGACGC", 3),
        FlankMatchSpec::new(b"GCGCTCTGTTGATAACTCCGGATC", b"AGATCGGAAGAGCGTCGTGCTATA", 3),
        true,
    );

    let specs = vec![&facs_spec, &syntf_spec];

    let fastq_filename = "/mnt/ingolialab/ingolia/Prog/pool-analysis/PacBio190731/samples.fastq";
    let out_base = "./pacbio-190731";

    pacbio_reads(specs.as_slice(), fastq_filename, out_base).unwrap_or_else(|err| {
        std::io::stderr()
            .write(format!("{}\n", err).as_bytes())
            .unwrap();
        std::process::exit(1);
    });
}
