use rust_htslib::bam;
use rust_htslib::prelude::*;

use errors::*;

pub fn barcode(qname: &[u8]) -> Result<&[u8]> {
    qname.split(|&ch| ch == b'_').next()
        .ok_or_else(|| ErrorKind::NoBarcode(qname.to_vec()).into())
}

pub struct BarcodeGroups<'a> {
    bam_reader: &'a mut bam::Reader,
    next_record: Option<bam::Record>,
}

impl <'a> BarcodeGroups<'a> {
    pub fn new(bam_reader: &'a mut bam::Reader) -> Result<Self> {
        let mut bg = BarcodeGroups{ bam_reader: bam_reader, next_record: None };
        bg.next_record = bg.read_next_record()?;
        Ok( bg )
    }

    fn read_next_record(&mut self) -> Result<Option<bam::Record>> {
        let mut rec = bam::Record::new();
        match self.bam_reader.read(&mut rec) {
            Ok( () ) => Ok( Some(rec) ),
            Err( bam::ReadError::NoMoreRecord ) => Ok( None ),
            Err( e ) => Err( e.into() ),
        }
    }

    fn barcode_group(&mut self, curr: bam::Record) -> Result<(Vec<u8>, Vec<bam::Record>)>
    {
        let curr_bc = barcode(curr.qname())?.to_vec();
        let mut bc_group = Vec::new();
        bc_group.push(curr);
        
        loop {
            let next = self.read_next_record()?;
            if let Some(rec) = next {
                if rec.qname().starts_with(&curr_bc) {
                    bc_group.push(rec);
                } else {
                    self.next_record = Some(rec);
                    break;
                }
            } else {
                self.next_record = None;
                break;
            }
        }
        
        Ok( (curr_bc, bc_group) )
    }
}

impl <'a> Iterator for BarcodeGroups<'a> {
    type Item = Result<(Vec<u8>, Vec<bam::Record>)>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(curr) = self.next_record.take() {
            Some(self.barcode_group(curr))
        } else {
            None
        }
    }
}


