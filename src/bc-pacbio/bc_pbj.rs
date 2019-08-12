extern crate barcode_assign;
extern crate clap;
extern crate failure;

use std::io::Write;

use clap::{App, Arg};

use barcode_assign::pacbio_join::CLI;

fn main() {
    let matches = App::new("bc-pbj")
        .version("1.0")
        .author("Nick Ingolia <ingolia@berkeley.edu>")
        .about("Collate PacBio alignments to assign barcodes to fragments")
        .arg(
            Arg::with_name("inbase")
                .short("i")
                .long("inbase")
                .value_name("INBASE")
                .help("Base name for input files")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("outbase")
                .short("o")
                .long("outbase")
                .value_name("OUTBASE")
                .help("Base name for output files")
                .takes_value(true)
                .required(true),
        )
        .get_matches();

    let cli = CLI {
        input_base: matches.value_of("inbase").unwrap().to_string(),
        inserts_good_file: None,
        frags_aligned_file: None,
        output_base: matches.value_of("outbase").unwrap().to_string(),
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
