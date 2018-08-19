use rust_htslib::bam;

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

            let cigar_view = r.cigar();
            let cigar: Vec<bam::record::Cigar> = cigar_view.iter().map(|&c| c.clone()).collect();
            
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
    cigar: Vec<bam::record::Cigar>,
    md: Vec<u8>
}

impl AssignMatch {
    pub fn tid(&self) -> u32 { self.tid }
    
    pub fn target(&self, targets: &[&str]) -> String {
        targets.get(self.tid as usize).unwrap_or(&"???").to_string()
    }

    pub fn pos(&self) -> i32 { self.pos }
    pub fn is_reverse(&self) -> bool { self.is_reverse }
    pub fn cigar(&self) -> &[bam::record::Cigar] { self.cigar.as_slice() }
    pub fn md(&self) -> &[u8] { self.md.as_slice() }
}


