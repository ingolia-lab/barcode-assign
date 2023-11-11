extern crate barcode_assign;
extern crate clap;

use barcode_assign::bc_umi::*;
use clap::{App, Arg};

fn main() {
    let matches = App::new("bc-umi")
        .version("1.0")
        .author("Nick Ingolia <ingolia@berkeley.edu>")
        .about("Combine barcodes by UMI")
        .arg(
            Arg::with_name("fastq")
                .short("f")
                .long("fastq")
                .value_name("BARCODE-FQ")
                .help("FastQ file of barcode sequences")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("umi")
                .short("u")
                .long("umi")
                .value_name("UMI-PREFIX")
                .help("Prefix for UMI in header")
                .takes_value(true)
                .default_value("umi="),
        )
        .arg(
            Arg::with_name("output")
                .short("o")
                .long("output")
                .value_name("OUTPUT-TXT")
                .help("Tab-delimited text file of barcode counts")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("neighborhood")
                .short("n")
                .long("neighborhood")
                .value_name("NBHD_BASE")
                .help("Analyze barcode mutation neighborhoods")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("dedup-stats")
                .short("d")
                .long("dedup-stats")
                .value_name("DEDUP_BASE")
                .help("Deduplication statistics")
                .takes_value(true),
        )
        .get_matches();

    let config = Config {
        barcode_fastq: matches.value_of("fastq").unwrap().to_string(),
        umi_prefix: matches.value_of("umi").unwrap().to_string(),
        out_barcodes: matches.value_of("output").unwrap().to_string(),
        dedup_stats: matches.value_of("dedup-stats").map(|s| String::from(s)),
        neighborhood: matches.value_of("neighborhood").map(|s| String::from(s)),
    };

    match bc_umi(config) {
        Ok(_) => (),
        Err(e) => {
            eprint!("bc_umi: failed due to error:\n{:#}\n", e);
        }
    }
}
