extern crate barcode_assign;
extern crate failure;
extern crate rust_htslib;

use std::collections::HashMap;
use std::io::{BufRead, BufReader, Write};

use rust_htslib::bam;
use rust_htslib::prelude::*;

fn main() {
    pacbio_join().unwrap_or_else(|err| {
        std::io::stderr()
            .write(format!("{}\n", err).as_bytes())
            .unwrap();
        std::process::exit(1);
    });
}

fn pacbio_join() -> Result<(), failure::Error> {
    let mut read_to_barcode = HashMap::new();
    let mut read_to_aligns = HashMap::new();
    let mut barcode_to_reads = HashMap::new();

    let barcodes_in = BufReader::new(std::fs::File::open("read-inserts-good.txt")?);
    for line_res in barcodes_in.lines() {
        let line = line_res?;
        let fields: Vec<&str> = line.split("\t").collect();
        if fields.len() < 4 {
            println!("Bad barcode line {:?}", line);
        } else {
            if let Some(_prev) = read_to_barcode.insert(
                fields[0].to_string(),
                (fields[3].to_string(), fields[1].to_string()),
            ) {
                println!("Duplicated read name {:?}", line);
            }

            let reads = barcode_to_reads
                .entry((fields[3].to_string(), fields[1].to_string()))
                .or_insert(Vec::new());
            reads.push(fields[0].to_string());
        }
    }

    let mut aligns_in = bam::Reader::from_path(&"frags_aligned.bam")?;
    let target_names_res: Result<Vec<_>, _> = aligns_in
        .header()
        .target_names()
        .into_iter()
        .map(|t| String::from_utf8(t.to_vec()))
        .collect();
    let target_names = target_names_res?;
    for res in aligns_in.records() {
        let r = res?;
        let qname = String::from_utf8(r.qname().to_vec())?;

        let read = &qname[0..qname.rfind("/").unwrap_or(qname.len())];
        let aligns = read_to_aligns.entry(read.to_string()).or_insert(Vec::new());
        aligns.push(r);
    }

    for ref mut aligns in read_to_aligns.values_mut() {
        aligns.sort_by_key(&align_sort_key);
    }
    
    let mut align_all_out = std::fs::File::create("read-aligns-all.txt")?;
    let mut align_unique_out = std::fs::File::create("read-aligns-unique.txt")?;
    
    for (&ref read, &(ref barcode, ref library)) in read_to_barcode.iter() {
        write!(align_all_out, "{}\t{}\t", barcode, library)?;
        if let Some(aligns) = read_to_aligns.get(read) {
            let status = if aligns.len() == 1 { "Unique" } else { "Multi" };
            let terse: Result<Vec<_>,_> = aligns.iter().map(|a| terse_align(&target_names, a)).collect();
            write!(align_all_out, "{}\t{}\n", status, terse?.join("\t"))?;

            if aligns.len() == 1 {
                write!(align_unique_out,
                       "{}\t{}\t{}\n",
                       barcode, library,
                       format_align(&target_names, &aligns[0])?)?;
            }
        } else {
            write!(align_all_out, "None\n")?;
        }
    }

    let mut barcode_all_out = std::fs::File::create("barcode-assign-all.txt")?;
    let mut barcode_unambig_out = std::fs::File::create("barcode-assign-unambig.txt")?;
    let mut barcode_unique_out = std::fs::File::create("barcode-assign-unique.txt")?;

    for (&(ref barcode, ref library), &ref reads) in barcode_to_reads.iter() {
        let aligns: Vec<Vec<bam::Record>> = reads
            .iter()
            .map(|r| {
                read_to_aligns
                    .get(r)
                    .map_or(Vec::new(), |alns| alns.clone())
            })
            .collect();

        write!(
            barcode_all_out,
            "{}\t{}\t{}\t",
            barcode, library, reads.len())?;
        
        if is_ambiguous(&aligns) {
            write!(barcode_all_out,
                   "Ambig\t{}\n",
                   format_ambiguous(&target_names, &aligns)?)?;
            continue;
        }

        let unambig = aligns.first().ok_or(failure::err_msg("Empty read set"))?;

        let status = if unambig.len() == 0 { "None" } else if unambig.len() > 1 { "Multi" } else { "Unique" };

        write!(barcode_all_out, "{}\n", status)?;
        
        let terse_res: Result<Vec<_>,_> = unambig.iter().map(|aln| terse_align(&target_names, aln)).collect();
        
        write!(
            barcode_unambig_out,
            "{}\t{}\t{}\t{}\t{}\n",
            barcode, library, 
            reads.len(), status,
            terse_res?.join("\t"))?;

        if unambig.len() == 1 {
            write!(
                barcode_unique_out,
                "{}\t{}\t{}\t{}\n",
                barcode, library,
                reads.len(),
                format_align(&target_names, &unambig[0])?)?;
        }
    }

    Ok(())
}

fn format_align(names: &Vec<String>, r: &bam::Record) -> Result<String, failure::Error> {
    Ok(format!(
        "{}\t{}\t{}\t{}",
        names[r.tid() as usize],
        r.pos(),
        r.cigar().end_pos()?,
        if r.is_reverse() { "-" } else { "+" }
    ))
}

fn terse_align(names: &Vec<String>, r: &bam::Record) -> Result<String, failure::Error> {
    Ok(format!("{}:{}-{}({})",
               names[r.tid() as usize],
               r.pos(),
               r.cigar().end_pos()?,
               if r.is_reverse() { "-" } else { "+" }
    ))               
}

fn align_sort_key(r: &bam::Record) -> (i32, i32) {
    (r.tid(), r.pos())
}

const POS_TOL: i32 = 3;

fn is_ambiguous(aligns: &Vec<Vec<bam::Record>>) -> bool {
    fn align_equivalent(r0: &bam::Record, r1: &bam::Record) -> bool {
        r0.tid() == r1.tid() && (r0.pos() - r1.pos()).abs() < POS_TOL && r0.is_reverse() == r1.is_reverse()
    }

    fn aligns_equivalent(alns0: &Vec<bam::Record>, alns1: &Vec<bam::Record>) -> bool {
        alns0.len() == alns1.len() && alns0.iter().zip(alns1.iter()).all(|(r0, r1)| align_equivalent(r0, r1))
    }
    
    let mut aligns_iter = aligns.iter();

    if let Some(aligns0) = aligns_iter.next() {
        aligns_iter.any(|aligns1| !aligns_equivalent(aligns0, aligns1))
    } else {
        false
    }
}

fn format_ambiguous(names: &Vec<String>, aligns: &Vec<Vec<bam::Record>>) -> Result<String, failure::Error> {
    fn format_read(names: &Vec<String>, alns: &Vec<bam::Record>) -> Result<String, failure::Error> {
        let terse_res: Result<Vec<_>,_> = alns.iter().map(|aln| { terse_align(names, aln) }).collect();
        Ok(terse_res?.join(";"))
    }

    let read_res: Result<Vec<_>,_> = aligns.iter().map(|alns| format_read(names, alns)).collect();
    Ok(read_res?.join("\t"))
}
