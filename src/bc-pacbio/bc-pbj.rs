extern crate barcode_assign;

use std::cmp::*;
use std::collections::HashMap;

use std::io::{BufRead,BufReader,Write};

fn main() {
    let mut read_to_barcode = HashMap::new();
    let mut read_to_aligns = HashMap::new();

    let barcodes_in = BufReader::new(std::fs::File::open("barcodes.txt").unwrap());
    for line_res in barcodes_in.lines() {
        let line = line_res.unwrap();
        let fields: Vec<&str> = line.split("\t").collect();
        if fields.len() < 4 {
            println!("Bad barcode line {:?}", line);
        } else {
            if let Some(_prev) = read_to_barcode.insert(fields[0].to_string(), (fields[1].to_string(), fields[3].to_string())) {
                println!("Duplicated read name {:?}", line);
            }
        }
    }

    let aligns_in = BufReader::new(std::fs::File::open("frags_aligned.txt").unwrap());
    for line_res in aligns_in.lines() {
        let line = line_res.unwrap();
        let fields: Vec<&str> = line.split(" ").collect();
        if fields.len() != 13 {
            println!("Bad alignment line {:?}", line);
        } else {
            let read = fields[0].rsplitn(2, "/").skip(1).next().unwrap();
            let aligns = read_to_aligns.entry(read.to_string()).or_insert(Vec::new());
            aligns.push(fields[1..fields.len()].join("\t"));
        }
    }

    let mut good_out = std::fs::File::create("barcode-aligns-good.txt").unwrap();
    let mut all_out = std::fs::File::create("barcode-aligns-all.txt").unwrap();
    for (&ref read, &(ref library, ref barcode)) in read_to_barcode.iter() {
        if let Some(aligns) = read_to_aligns.get(read) {
            if aligns.len() == 1 {
                write!(all_out, "{}\t{}\tUnique\n", barcode, library).unwrap();
                write!(good_out, "{}\t{}\t{}\n", barcode, library, aligns[0]).unwrap();
            } else {
                write!(all_out, "{}\t{}\tMulti\t{}\n", barcode, library, aligns.len()).unwrap();
            }
        } else {
            write!(all_out, "{}\t{}\tNone\n", barcode, library).unwrap();
        }
    }
}
