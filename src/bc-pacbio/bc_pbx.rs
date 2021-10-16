extern crate barcode_assign;
extern crate clap;

use std::io::Write;

use clap::{App, Arg};

use barcode_assign::pacbio_extract::CLI;

fn main() {
    let matches = App::new("bc-pbx")
        .version("1.0")
        .author("Nick Ingolia <ingolia@berkeley.edu>")
        .about("Extract multiple barcodes and inserts from PacBio CCS")
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
            Arg::with_name("config")
                .short("c")
                .long("config")
                .value_name("CONFIG.TXT")
                .help("Filename of tab-delimited library specifications")
                .takes_value(true)
                .required(true),
        )
        .get_matches();

    let cli = CLI {
        input_bam: matches.value_of("input").unwrap().to_string(),
        config_file: matches.value_of("config").unwrap().to_string(),
    };

    match cli.run() {
        Ok(_) => (),
        Err(err) => {
            std::io::stderr()
                .write(format!("{:?}\n", err).as_bytes())
                .unwrap();
            std::process::exit(1);
        }
    }
}
