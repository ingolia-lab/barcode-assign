extern crate barcode_assign;
extern crate failure;
//extern crate rust_htslib;

use std::cmp::*;
use std::collections::HashMap;
use std::io::{BufRead,BufReader,Write};

use barcode_assign::blasr::*;
//use rust_htslib::bam;
//use rust_htslib::prelude::*;

fn main() {
    pacbio_join().unwrap_or_else(|err| {
        std::io::stderr().write(format!("{}\n", err).as_bytes()).unwrap();
        std::process::exit(1);
    });
}

fn pacbio_join() -> Result<(), failure::Error> {
    let mut read_to_barcode = HashMap::new();
    let mut read_to_aligns = HashMap::new();
    let mut barcode_to_reads = HashMap::new();

    let barcodes_in = BufReader::new(std::fs::File::open("barcodes.txt")?);
    for line_res in barcodes_in.lines() {
        let line = line_res?;
        let fields: Vec<&str> = line.split("\t").collect();
        if fields.len() < 4 {
            println!("Bad barcode line {:?}", line);
        } else {
            if let Some(_prev) = read_to_barcode.insert(fields[0].to_string(), (fields[1].to_string(), fields[3].to_string())) {
                println!("Duplicated read name {:?}", line);
            }

            let reads = barcode_to_reads.entry((fields[1].to_string(), fields[3].to_string())).or_insert(Vec::new());
            reads.push(fields[0].to_string());
        }
    }

    let aligns_in = BufReader::new(std::fs::File::open("frags_aligned.txt")?);
    for line_res in aligns_in.lines() {
        let line = line_res?;
        let blasr = BlasrAlign::new(&line)?;
        let read = blasr.query_name_trimmed().to_string();
        let aligns = read_to_aligns.entry(read.to_string()).or_insert(Vec::new());
        aligns.push(blasr);
    }

    let mut good_out = std::fs::File::create("barcode-aligns-good.txt")?;
    let mut all_out = std::fs::File::create("barcode-aligns-all.txt")?;

    for (&ref read, &(ref library, ref barcode)) in read_to_barcode.iter() {
        if let Some(aligns) = read_to_aligns.get(read) {
            if aligns.len() == 1 {
                write!(all_out, "{}\t{}\tUnique\n", barcode, library)?;
                write!(good_out, "{}\t{}\t{}\n", barcode, library, aligns[0])?;
            } else {
                write!(all_out, "{}\t{}\tMulti\t{}\n", barcode, library, aligns.len())?;
            }
        } else {
            write!(all_out, "{}\t{}\tNone\n", barcode, library)?;
        }
    }

    let mut assign_out = std::fs::File::create("barcode-assign-good.txt")?;
    let mut barcode_out = std::fs::File::create("barcode-assign-all.txt")?;

    for (&(ref library, ref barcode), &ref reads) in barcode_to_reads.iter() {
        let aligns: Vec<Vec<BlasrAlign>> = reads.iter().map(|r| read_to_aligns.get(r).map_or(Vec::new(), |alns| alns.clone())).collect();

        if aligns.iter().any(|read_align| read_align.len() > 1) {
            write!(barcode_out, "{}\t{}\t{}\tMulti\n", library, barcode, reads.len())?;
            continue;
        }

        let read_aligns: Vec<&BlasrAlign> = aligns.iter().filter_map(|read_align| read_align.first()).collect();
        
        if read_aligns.len() == 0 {
            write!(barcode_out, "{}\t{}\t{}\tNone\n", library, barcode, reads.len())?;
            continue;            
        }

        if read_aligns.iter().all(|ra| ra.target_name() == read_aligns[0].target_name()
                                  && ra.target_start() == read_aligns[0].target_start()
                                  && ra.target_end() == read_aligns[0].target_end()) {
            write!(assign_out, "{}\t{}\t{}\t{}\n", library, barcode, reads.len(), read_aligns[0])?;
            write!(barcode_out, "{}\t{}\t{}\tGood\n", library, barcode, reads.len())?;
        } else {
            write!(barcode_out, "{}\t{}\t{}\tMixed\n", library, barcode, reads.len())?;
        }
        
    }

    Ok(())
}
