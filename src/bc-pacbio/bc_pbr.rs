extern crate barcode_assign;
extern crate clap;

use std::io::Write;

use clap::{App,Arg};

use barcode_assign::pacbio_reads::CLI;

fn main() {
    let matches = App::new("bc-pbr")
        .version("1.0")
        .author("Nick Ingolia <ingolia@berkeley.edu>")
        .about("Extract barcodes and inserts from PacBio CCS")
        .arg(
            Arg::with_name("input")
                .short("i")
                .long("input")
                .value_name("INPUT.FQ")
                .help("FastQ file of PacBio CCS")
                .takes_value(true)
                .required(true))
        .arg(Arg::with_name("outbase")
             .short("o")
             .long("outbase")
             .value_name("OUTBASE")
             .help("Base name for constructing output names")
             .takes_value(true)
             .required(true))
        .arg(Arg::with_name("libspecs")
             .short("l")
             .long("libspecs")
             .value_name("LIB_SPECS.TXT")
             .help("Filename of tab-delimited library specifications")
             .takes_value(true)
             .required(true))
        .arg(Arg::with_name("max_errors")
             .short("e")
             .long("max-errors")
             .value_name("NERRORS")
             .help("Maximum number of errors allowed per match")
             .takes_value(true)
             .default_value("3"))
        .arg(Arg::with_name("matches")
             .short("m")
             .long("matches")
             .help("Write verbose linker matching information"))
        .get_matches();

    let cli = CLI {
        input_fastq: matches.value_of("input").unwrap().to_string(),
        input_specs_file: matches.value_of("libspecs").unwrap().to_string(),
        output_base: matches.value_of("outbase").unwrap().to_string(),
        output_file_frags:None,
        output_file_inserts:None,
        output_file_fates:None,
        output_file_matching:None,
        output_matching: matches.occurrences_of("matches") > 0,
        max_errors_str: matches.value_of("max_errors").unwrap().to_string(),
    };

    match cli.run() {
        Ok(_) => (),
        Err(err) => {
            std::io::stderr()
                .write(format!("{}\n", err).as_bytes())
                .unwrap();
            std::process::exit(1);
        }
    }

    // let facs_spec = LibSpec::new(
    //     "FACS",
    //     FlankMatchSpec::new(b"GCTCGGAGATGTGTATAAGAGACAG", b"CTGTCTCTTATACACATCTGACGCC", 3),
    //     FlankMatchSpec::new(b"CTATAGCACGACGCTCTTCCGATCT", b"GATCCTGTAGCCCTAGACTTGATAG", 3),
    //     false,
    // );

    // let syntf_spec = LibSpec::new(
    //     "SynTF",
    //     FlankMatchSpec::new(b"CTCGGAGATGTGTATAAGAGACAG", b"CTGTCTCTTATACACATCTGACGC", 3),
    //     FlankMatchSpec::new(b"GCGCTCTGTTGATAACTCCGGATC", b"AGATCGGAAGAGCGTCGTGCTATA", 3),
    //     true,
    // );

    // let specs = vec![facs_spec, syntf_spec];

    // let fastq_filename = "/mnt/ingolialab/ingolia/Prog/pool-analysis/PacBio190731/samples.fastq";
    // let out_base = "./pacbio-190731";
}
