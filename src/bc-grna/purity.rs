use rust_htslib::bam;

use assign::ReadAssign;
use errors::*;

#[derive(Debug,Clone,Hash,PartialEq,Eq)]
pub struct Purity {
    assign_counts: Vec<(ReadAssign,usize)>
}    
    
impl Purity {
    pub fn new<'a, I>(r_iter: I) -> Result<Self>
        where I: Iterator<Item = &'a bam::Record>
    {
        let mut cts = Vec::new();
        
        for r in r_iter {
            let asn = ReadAssign::new(r)?;
            
            let mut extant = false;
            for &mut(ref mut a, ref mut n) in cts.iter_mut() {
                if asn == *a {
                    *n += 1;
                    extant = true;
                }
            }
            if !extant {
                cts.push( (asn, 1) );
            }
        }
        
        cts.sort_by_key(|&(ref _a, ref n)| - (*n as isize));
        
        Ok( Purity { assign_counts: cts } )
    }

    pub fn primary_assign(&self) -> ReadAssign {
        self.assign_counts.first()
            .map_or(ReadAssign::NoMatch, |&(ref a, ref _n)| a.clone())
    }
    
    pub fn purity(&self) -> f64 {
        self.assign_counts.first()
            .map_or(0.0, |&(ref _a, ref nprim)| (*nprim as f64) / (self.n_total() as f64))
    }
    
    pub fn n_total(&self) -> usize {
        self.assign_counts.iter().map(|&(ref _a, ref n)| *n).sum()
    }

    pub fn n_primary(&self) -> usize {
        self.assign_counts.first()
            .map_or(0, |&(ref _a, ref n)| *n)
    }
    
    pub fn n_no_match(&self) -> usize {
        self.assign_counts.iter()
            .filter(|&(ref a, ref _n)| a.is_no_match())
            .map(|&(ref _a, ref n)| *n).sum()
    }

    pub fn n_other_match(&self) -> usize {
        self.assign_counts.iter().skip(1)
            .filter(|&(ref a, ref _n)| !a.is_no_match())
            .map(|&(ref _a, ref n)| *n).sum()
    }

    pub fn header() -> String {
        "barcode\tprimary\tother_match\tno_match".to_string()
    }
    
    pub fn line(&self, bc_str: &str) -> String {
        format!("{}\t{}\t{}\t{}\n", bc_str,
                self.n_primary(), self.n_other_match(), self.n_no_match())
    }
}

