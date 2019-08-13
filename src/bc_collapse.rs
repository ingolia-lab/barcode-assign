extern crate barcode_assign;
extern crate bio;
extern crate clap;
#[macro_use]
extern crate failure;

use clap::{App, Arg};

use barcode_assign::collapse::CLI;

fn main() {
    let matches = App::new("bc-collapse")
        .version("1.0")
        .author("Nick Ingolia <ingolia@berkeley.edu>")
        .about("Collapse barcode sequences")
        .arg(
            Arg::with_name("fastq_in")
                .short("f")
                .long("fastq")
                .value_name("BARCODE.FASTQ")
                .help("FastQ file of barcode sequences")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("output_base")
                .short("o")
                .long("outbase")
                .value_name("OUTPUT_BASE")
                .help("Base name for output files")
                .takes_value(true)
                .required(true),
        )
        .get_matches();

    let cli = CLI {
        fastq_in: matches.value_of("fastq_in").unwrap().to_string(),
        output_base: matches.value_of("output_base").unwrap().to_string(),
    };

    match cli.run() {
        Ok(_) => (),
        Err(e) => panic!(e),
    }
}

