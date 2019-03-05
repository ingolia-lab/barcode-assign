extern crate barcode_assign;
extern crate clap;

use clap::{Arg, App};
use barcode_assign::bc_count::*;

fn main() {
    let matches = App::new("bc-count")
        .version("1.0")
        .author("Nick Ingolia <ingolia@berkeley.edu>")
        .about("Count barcode sequences")
        .arg(Arg::with_name("fastq")
             .short("f")
             .long("fastq")
             .value_name("BARCODE-FQ")
             .help("FastQ file of barcode sequences")
             .takes_value(true)
             .required(true))
        .arg(Arg::with_name("output")
             .short("o")
             .long("output")
             .value_name("OUTPUT-TXT")
             .help("Tab-delimited text file of barcode counts")
             .takes_value(true)
             .required(true))
        .get_matches();

    let config = Config {
        barcode_fastq: matches.value_of("fastq").unwrap().to_string(),
        out_barcodes: matches.value_of("output").unwrap().to_string(),
        freq_filename: None,
    };

    match bc_count(config) {
        Ok(_) => (),
        Err(e) => panic!(e),
    }
}
