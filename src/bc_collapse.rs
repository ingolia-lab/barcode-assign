extern crate barcode_assign;
extern crate bio;
extern crate clap;
extern crate failure;

use std::io::Read;

use bio::io::fastq;
use barcode_assign::collapse::*;
use clap::{App, Arg};

fn main() {
    let matches = App::new("bc-collapse")
        .version("1.0")
        .author("Nick Ingolia <ingolia@berkeley.edu>")
        .about("Collapse barcode sequences")
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
            Arg::with_name("output")
                .short("o")
                .long("output")
                .value_name("OUTPUT-TXT")
                .help("Tab-delimited text file of barcode counts")
                .takes_value(true)
                .required(true),
        )
        .get_matches();

    let cli = CLI {
        barcode_fastq: matches.value_of("fastq").unwrap().to_string(),
        output: matches.value_of("output").unwrap().to_string(),
    };

    match cli.run() {
        Ok(_) => (),
        Err(e) => panic!(e),
    }
}

struct CLI {
    barcode_fastq: String,
    output: String,
}

impl CLI {
    fn run(&self) -> Result<(), failure::Error> {
        let mut nbhds = BarcodeNbhdMap::new(1);

        let reader: Box<Read> = if self.barcode_fastq == "-" {
            Box::new(std::io::stdin())
        } else {
            Box::new(std::fs::File::open(&self.barcode_fastq)?)
        };
        let barcode_reader = fastq::Reader::new(reader);

        let mut last_report = std::time::Instant::now();
        let mut last_recno = 0;
        
        for (recno, rec_res) in barcode_reader.records().enumerate() {
            let rec = rec_res?;

            if rec.seq().len() < 24 || rec.seq().len() > 26 {
//                println!("Skipping {} because of length",
//                         String::from_utf8_lossy(rec.seq()));
                continue;
            }
            
            nbhds.insert(rec.seq());

            let last_duration = std::time::Instant::now().duration_since(last_report);
            if last_duration.as_secs() > 0
            {
                let recs = recno - last_recno;
                println!("{} usec, {} records, {} records / second ({} total)",
                         last_duration.as_micros(),
                         recs,
                         1000000.0 * (recs as f64) / (last_duration.as_micros() as f64),
                         recno);
                last_report = std::time::Instant::now();
                last_recno = recno;
            }

            if recno > 10000000 {
                break;
            }
        }
        
        nbhds.write_nbhds(std::fs::File::create(&self.output)?)
    }
}
