use std::collections::HashMap;
use std::fs::File;
use std::io::{self, Read, Write};

use bio::io::fastq;

//use counts::SampleCounts;
//use neighborhood::*;

#[derive(Debug)]
pub struct Config {
    pub barcode_fastq: String,
    pub out_barcodes: String,
//    pub neighborhood: Option<String>,
}

pub fn bc_umi(config: Config) -> Result<(), failure::Error> {
    let reader: Box<dyn Read> = if config.barcode_fastq == "-" {
        Box::new(io::stdin())
    } else {
        Box::new(File::open(&config.barcode_fastq)?)
    };
    let barcode_reader = fastq::Reader::new(reader);

    let mut umi_counts = UmiCounts::new();

    for recres in barcode_reader.records() {
        let rec = recres?;
        let umi = UmiCounts::find_umi(rec.desc().unwrap()).unwrap();
        umi_counts.count_one(rec.seq(), umi.as_bytes());
    }

    let writer: Box<dyn Write> = if config.out_barcodes == "-" {
        Box::new(io::stdout())
    } else {
        Box::new(File::create(&config.out_barcodes)?)
    };
    umi_counts.write(writer)?;

    Ok(())
}

/// Tabulation of barcode counts in a sample
#[derive(Debug, Clone)]
pub struct UmiCounts(HashMap<Vec<u8>, HashMap<Vec<u8>, usize>>);

impl UmiCounts {
    pub fn new() -> Self
    {
        UmiCounts(HashMap::new())
    }
    
    pub fn count_one(&mut self, barcode: &[u8], umi: &[u8]) -> ()
    {
        let barcode_umis = if let Some(existing) = self.0.get_mut(barcode) {
            existing
        } else {
            self.0.insert(barcode.to_owned(), HashMap::new());
            self.0.get_mut(barcode).unwrap()
        };

        let umi_count = if let Some(existing) = barcode_umis.get_mut(umi) {
            existing
        } else {
            barcode_umis.insert(umi.to_owned(), 0);
            barcode_umis.get_mut(umi).unwrap()
        };

        *umi_count += 1;        
    }

    pub fn find_umi(desc: &str) -> Option<&str> {
        let (_, rest) = desc.split_once("umi=")?;
        rest.split_whitespace().next()
    }

    pub fn write<W: Write>(&self, umi_out: W) -> Result<(), failure::Error> {
        let mut out = std::io::BufWriter::new(umi_out);

        for (barcode, umis) in self.0.iter() {
            write!(out, "{}", String::from_utf8_lossy(barcode))?;
            let mut umi_counts = umis.iter().map(|uc| *uc.1).collect::<Vec<usize>>();
            umi_counts.sort();
            umi_counts.reverse();

            write!(out, "\t{}\t{}\t",
                   umi_counts.iter().copied().sum::<usize>(), umi_counts.len())?;

            write!(out, "{}\t",
                   umi_counts[umi_counts.len() / 2])?;
            
            for count in umi_counts.iter() {
                write!(out, "{},", count)?;
            }
            write!(out, "\n")?;
        }
        
        Ok(())
    }
}
