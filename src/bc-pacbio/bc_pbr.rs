extern crate barcode_assign;
extern crate clap;

use std::io::Write;

use clap::{App, Arg};

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
                .value_name("INPUT.BAM")
                .help("BAM file of PacBio CCS")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("outbase")
                .short("o")
                .long("outbase")
                .value_name("OUTBASE")
                .help("Base name for constructing output names")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("libspecs")
                .short("l")
                .long("libspecs")
                .value_name("LIB_SPECS.TXT")
                .help("Filename of tab-delimited library specifications")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("max_errors")
                .short("e")
                .long("max-errors")
                .value_name("NERRORS")
                .help("Maximum number of errors allowed per match")
                .takes_value(true)
                .default_value("3"),
        )
        .arg(
            Arg::with_name("matches")
                .short("m")
                .long("matches")
                .help("Write verbose linker matching information"),
        )
        .arg(
            Arg::with_name("frags")
                .long("frags")
                .value_name("FRAGS.FA")
                .help("Filename of fragment sequence output file")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("inserts")
                .long("inserts")
                .value_name("INSERTS.TXT")
                .help("Filename of insert sequence output file")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("fates")
                .long("fates")
                .value_name("FATES.TXT")
                .help("Filename of read fate output file")
                .takes_value(true),
        )
        .get_matches();

    let cli = CLI {
        input_bam: matches.value_of("input").unwrap().to_string(),
        input_specs_file: matches.value_of("libspecs").unwrap().to_string(),
        output_base: matches.value_of("outbase").unwrap().to_string(),
        output_file_frags: matches.value_of("frags").map(String::from),
        output_file_inserts: matches.value_of("inserts").map(String::from),
        output_file_fates: matches.value_of("fates").map(String::from),
        output_file_matching: None,
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
}
