use std::collections::HashMap;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

use rust_htslib::bam;
use rust_htslib::bam::Read;

#[derive(Debug)]
pub struct CLI {
    pub input_base: String,
    pub inserts_good_file: Option<String>,
    pub frags_aligned_file: Option<String>,
    pub output_base: String,
}

impl CLI {
    pub fn inserts_good_file(&self) -> PathBuf {
        self.inserts_good_file.as_ref().map_or_else(
            || self.input_filename("-read-inserts-good.txt"),
            |f| PathBuf::from(f),
        )
    }

    pub fn frags_aligned_file(&self) -> PathBuf {
        self.inserts_good_file.as_ref().map_or_else(
            || self.input_filename("-frags-aligned.bam"),
            |f| PathBuf::from(f),
        )
    }

    pub fn input_filename(&self, name: &str) -> PathBuf {
        let base_ref: &Path = self.input_base.as_ref();
        let mut namebase = base_ref
            .file_name()
            .map_or(std::ffi::OsString::new(), std::ffi::OsStr::to_os_string);
        namebase.push(name);
        base_ref.with_file_name(namebase)
    }

    pub fn output_filename(&self, name: &str) -> PathBuf {
        let base_ref: &Path = self.output_base.as_ref();
        let mut namebase = base_ref
            .file_name()
            .map_or(std::ffi::OsString::new(), std::ffi::OsStr::to_os_string);
        namebase.push(name);
        base_ref.with_file_name(namebase)
    }

    pub fn outputs(&self) -> Result<Outputs, failure::Error> {
        Ok(Outputs {
            read_aligns_all: Outputs::output(self.output_filename("-read-aligns-all.txt"))?,
            read_aligns_unique: Outputs::output(self.output_filename("-read-aligns-unique.txt"))?,
            barcode_assign_all: Outputs::output(self.output_filename("-barcode-assign-all.txt"))?,
            barcode_assign_unambig: Outputs::output(
                self.output_filename("-barcode-assign-umabig.txt"),
            )?,
            barcode_assign_unique: Outputs::output(
                self.output_filename("-barcode-assign-unique.txt"),
            )?,
            barcode_assign_bed: Outputs::output(self.output_filename("-barcode-assign.bed"))?,
        })
    }

    pub fn run(&self) -> Result<(), failure::Error> {
        let mut read_inserts_good = std::fs::File::open(self.inserts_good_file())?;
        let frags_aligned = bam::Reader::from_path(self.frags_aligned_file())?;

        let mut outputs = self.outputs()?;

        pacbio_join(&mut read_inserts_good, frags_aligned, &mut outputs)
    }
}

pub struct Outputs {
    read_aligns_all: Box<dyn Write>,
    read_aligns_unique: Box<dyn Write>,

    barcode_assign_all: Box<dyn Write>,
    barcode_assign_unambig: Box<dyn Write>,
    barcode_assign_unique: Box<dyn Write>,
    barcode_assign_bed: Box<dyn Write>,
}

impl Outputs {
    pub fn read_aligns_all(&mut self) -> &mut dyn Write {
        &mut self.read_aligns_all
    }
    pub fn read_aligns_unique(&mut self) -> &mut dyn Write {
        &mut self.read_aligns_unique
    }

    pub fn barcode_assign_all(&mut self) -> &mut dyn Write {
        &mut self.barcode_assign_all
    }
    pub fn barcode_assign_unambig(&mut self) -> &mut dyn Write {
        &mut self.barcode_assign_unambig
    }
    pub fn barcode_assign_unique(&mut self) -> &mut dyn Write {
        &mut self.barcode_assign_unique
    }
    pub fn barcode_assign_bed(&mut self) -> &mut dyn Write {
        &mut self.barcode_assign_bed
    }

    pub fn output<P: AsRef<Path>>(filename: P) -> Result<Box<dyn Write>, failure::Error> {
        Ok(Box::new(std::fs::File::create(filename)?))
    }
}

pub fn pacbio_join<R: std::io::Read>(
    read_inserts_good: R,
    mut frags_aligned: bam::Reader,
    outputs: &mut Outputs,
) -> Result<(), failure::Error> {
    let mut read_to_barcode = HashMap::new();
    let mut read_to_aligns = HashMap::new();
    let mut barcode_to_reads = HashMap::new();

    let barcodes_in = BufReader::new(read_inserts_good);
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

    let target_names_res: Result<Vec<_>, _> = frags_aligned
        .header()
        .target_names()
        .into_iter()
        .map(|t| String::from_utf8(t.to_vec()))
        .collect();
    let target_names = target_names_res?;
    for res in frags_aligned.records() {
        let r = res?;
        let qname = String::from_utf8(r.qname().to_vec())?;

        let read = &qname[0..qname.rfind("/").unwrap_or(qname.len())];
        let aligns = read_to_aligns.entry(read.to_string()).or_insert(Vec::new());
        aligns.push(r);
    }

    for ref mut aligns in read_to_aligns.values_mut() {
        aligns.sort_by_key(&align_sort_key);
    }

    for (&ref read, &(ref barcode, ref library)) in read_to_barcode.iter() {
        write!(outputs.read_aligns_all(), "{}\t{}\t", barcode, library)?;
        if let Some(aligns) = read_to_aligns.get(read) {
            let status = if aligns.len() == 1 { "Unique" } else { "Multi" };
            let terse: Result<Vec<_>, _> = aligns
                .iter()
                .map(|a| terse_align(&target_names, a))
                .collect();
            write!(
                outputs.read_aligns_all(),
                "{}\t{}\n",
                status,
                terse?.join("\t")
            )?;

            if aligns.len() == 1 {
                write!(
                    outputs.read_aligns_unique(),
                    "{}\t{}\t{}\n",
                    barcode,
                    library,
                    format_align(&target_names, &aligns[0])?
                )?;
            }
        } else {
            write!(outputs.read_aligns_all(), "None\n")?;
        }
    }

    let empty = Vec::new();

    for (&(ref barcode, ref library), &ref reads) in barcode_to_reads.iter() {
        let aligns: Vec<&Vec<bam::Record>> = reads
            .iter()
            .map(|r| read_to_aligns.get(r).unwrap_or(&empty))
            .collect();

        write!(
            outputs.barcode_assign_all(),
            "{}\t{}\t{}\t",
            barcode,
            library,
            reads.len()
        )?;

        if is_ambiguous(&aligns) {
            write!(
                outputs.barcode_assign_all(),
                "Ambig\t{}\n",
                format_ambiguous(&target_names, &aligns)?
            )?;
            continue;
        }

        let unambig = aligns.first().ok_or(failure::err_msg("Empty read set"))?;

        let status = if unambig.len() == 0 {
            "None"
        } else if unambig.len() > 1 {
            "Multi"
        } else {
            "Unique"
        };

        write!(outputs.barcode_assign_all(), "{}\n", status)?;

        let terse_res: Result<Vec<_>, _> = unambig
            .iter()
            .map(|aln| terse_align(&target_names, aln))
            .collect();

        write!(
            outputs.barcode_assign_unambig(),
            "{}\t{}\t{}\t{}\t{}\n",
            barcode,
            library,
            reads.len(),
            status,
            terse_res?.join("\t")
        )?;

        if unambig.len() == 1 {
            let frag = &unambig[0];

            write!(
                outputs.barcode_assign_unique(),
                "{}\t{}\t{}\t{}\n",
                barcode,
                library,
                reads.len(),
                format_align(&target_names, frag)?
            )?;

            write!(
                outputs.barcode_assign_bed(),
                "{}\t{}\t{}\t{}_{}\t{}\t{}\n",
                target_names[frag.tid() as usize],
                frag.pos(),
                frag.cigar().end_pos(),
                barcode,
                library,
                reads.len(),
                if frag.is_reverse() { "-" } else { "+" }
            )?;
        }
    }

    Ok(())
}

fn format_align(names: &Vec<String>, r: &bam::Record) -> Result<String, failure::Error> {
    Ok(format!(
        "{}\t{}\t{}\t{}",
        names[r.tid() as usize],
        r.pos(),
        r.cigar().end_pos(),
        if r.is_reverse() { "-" } else { "+" }
    ))
}

fn terse_align(names: &Vec<String>, r: &bam::Record) -> Result<String, failure::Error> {
    Ok(format!(
        "{}:{}-{}({})",
        names[r.tid() as usize],
        r.pos(),
        r.cigar().end_pos(),
        if r.is_reverse() { "-" } else { "+" }
    ))
}

fn align_sort_key(r: &bam::Record) -> (i32, i64) {
    (r.tid(), r.pos())
}

const FRAC_TOL: f64 = 0.76;
const POS_TOL: i64 = 12;

fn is_ambiguous(aligns: &Vec<&Vec<bam::Record>>) -> bool {
    fn align_equivalent(r0: &bam::Record, r1: &bam::Record) -> bool {
        r0.tid() == r1.tid()
            && (r0.pos() - r1.pos()).abs() < POS_TOL
            && r0.is_reverse() == r1.is_reverse()
    }

    fn aligns_equivalent(alns0: &Vec<bam::Record>, alns1: &Vec<bam::Record>) -> bool {
        alns0.len() == alns1.len()
            && alns0
                .iter()
                .zip(alns1.iter())
                .all(|(r0, r1)| align_equivalent(r0, r1))
    }

    fn aligns_groups_insert<'a>(mut groups: Vec<(&'a Vec<bam::Record>, usize)>, aligns: &'a Vec<bam::Record>)
        -> Vec<(&'a Vec<bam::Record>, usize)> {
        for (group_aligns, group_count) in groups.iter_mut() {
            if aligns_equivalent(group_aligns, aligns) {
                *group_count += 1;
                return groups;
            }
        }
        groups.push((aligns, 1));
        groups
    }

    let groups = aligns.iter().fold(Vec::new(), |groups, ref aligns| aligns_groups_insert(groups, aligns));
    
    let count_first = groups.first().map_or(0, |(_aligns, count)| *count);
    let count_total: usize = groups.iter().map(|(_aligns, count)| *count).sum();

    if count_total > 0 {
        (count_first as f64) / (count_total as f64) < FRAC_TOL
    } else {
        false
    }
}

fn format_ambiguous(
    names: &Vec<String>,
    aligns: &Vec<&Vec<bam::Record>>,
) -> Result<String, failure::Error> {
    fn format_read(names: &Vec<String>, alns: &Vec<bam::Record>) -> Result<String, failure::Error> {
        if alns.len() == 0 {
            return Ok("none".to_string());
        }
        let terse_res: Result<Vec<_>, _> = alns.iter().map(|aln| terse_align(names, aln)).collect();
        Ok(terse_res?.join(";"))
    }

    let read_res: Result<Vec<_>, _> = aligns.iter().map(|alns| format_read(names, alns)).collect();
    Ok(read_res?.join("\t"))
}

pub fn input_filename<P: AsRef<Path>>(in_base: P, name: &str) -> PathBuf {
    let mut namebase = in_base
        .as_ref()
        .file_name()
        .map_or(std::ffi::OsString::new(), std::ffi::OsStr::to_os_string);
    namebase.push(name);
    in_base.as_ref().with_file_name(namebase)
}

pub fn output_filename<Q: AsRef<Path>>(out_base: Q, name: &str) -> PathBuf {
    let mut namebase = out_base
        .as_ref()
        .file_name()
        .map_or(std::ffi::OsString::new(), std::ffi::OsStr::to_os_string);
    namebase.push(name);
    out_base.as_ref().with_file_name(namebase)
}
