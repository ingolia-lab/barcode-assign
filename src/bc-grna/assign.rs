use std::str;

use rust_htslib::bam;
use rust_htslib::bam::record::{Cigar,CigarString};

use errors::*;

#[derive(Debug,Clone,Hash,PartialEq,Eq)]
pub enum ReadAssign {
    NoMatch,
    Match(AssignMatch)
}

impl ReadAssign {
    pub fn new(r: &bam::Record) -> Result<Self> {
        if r.tid() < 0 {
            Ok( ReadAssign::NoMatch )
        } else {
            let md_aux = r.aux(b"MD").ok_or_else(|| "No MD tag")?;
            let md: Vec<u8> = match md_aux {
                bam::record::Aux::String(md) => Ok( md.to_vec() ),
                _ => Err("MD tag not a string")
            }?;

            let cigar = (*r.cigar()).clone();
            
            let assign_match = AssignMatch { tid: r.tid() as u32, pos: r.pos(),
                                             is_reverse: r.is_reverse(),
                                             cigar: cigar, md: md };
            Ok( ReadAssign::Match(assign_match) )
        }
    }

    pub fn is_no_match(&self) -> bool {
        match self {
            ReadAssign::NoMatch => true,
            _ => false
        }
    }
}

#[derive(Debug,Clone,Hash,PartialEq,Eq)]
pub struct AssignMatch {
    tid: u32,
    pos: i32,
    is_reverse: bool,
    cigar: CigarString,
    md: Vec<u8>
}

impl AssignMatch {
    #[allow(dead_code)]
    pub fn tid(&self) -> u32 { self.tid }
    
    pub fn target(&self, targets: &[&str]) -> String {
        targets.get(self.tid as usize).unwrap_or(&"???").to_string()
    }

    pub fn pos(&self) -> i32 { self.pos }
    pub fn is_reverse(&self) -> bool { self.is_reverse }

    #[allow(dead_code)]
    pub fn cigar(&self) -> CigarString { self.cigar.clone() }
    pub fn cigar_string(&self) -> String { self.cigar.to_string() }
    pub fn md(&self) -> &[u8] { self.md.as_slice() }

    pub fn is_cigar_perfect(&self, len: u32) -> bool {
        let CigarString(ref v) = self.cigar;
        (v.len() == 1) && (v[0] == Cigar::Match(len))
    }

    pub fn is_md_perfect(&self, len: u32) -> bool {
        len.to_string().as_str().as_bytes() == self.md.as_slice()
    }

    pub fn header() -> String {
        "barcode\ttid\tpos\tcigar\tmd".to_string()
    }
    
    pub fn line(&self, bc_str: &str, targets: &[&str]) -> Result<String> {
        Ok( format!("{}\t{}\t{}\t{:?}\t{}",
                    bc_str, self.target(&targets), self.pos(),
                    self.cigar_string(), str::from_utf8(self.md())?) )
    }
}


