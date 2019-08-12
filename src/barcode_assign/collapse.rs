use std::cell::RefCell;
use std::collections::HashMap;
use std::io::Write;
use std::rc::Rc;

use bio::pattern_matching::myers::Myers;

// New barcode -- check first for exact match, then scan all existing barcodes
//   for near-matches
// Insert later barcodes with links to earlier barcodes?
//   Or, create barcode equivalency sets?

// No need for complexity during counting -- post-process key set in
// map to find neighborhoods. No need to handle insertions, can find
// all edges with substitutions and deletions.

// 1. Pick a node
//    a. remove from key set
//    b. initialize work stack with node
// 2. Handle a node from work stack
//    a. check key set for all near-neighbors
//       i. remove near-neighbor from key set
//       ii. push near-neighbor onto work stack
//    b. add node to neighborhood

pub struct BarcodeNbhdMap {
    nbhds: Vec<Rc<RefCell<BarcodeNbhd>>>,
    barcode_nbhds: HashMap<Vec<u8>, Rc<RefCell<BarcodeNbhd>>>,
    max_dist: u8,
}

impl BarcodeNbhdMap {
    pub fn new(max_dist: usize) -> Self {
        BarcodeNbhdMap {
            nbhds: Vec::new(),
            barcode_nbhds: HashMap::new(),
            max_dist: max_dist as u8,
        }
    }

    pub fn insert(&mut self, bc: &[u8]) -> () {
        if self.barcode_nbhds.contains_key(bc) {
            self.barcode_nbhds.get(bc).unwrap().borrow_mut().insert(bc);
        } else {
            let mut subst_iter = Substitutions::new(bc).chain(Deletions::new(bc));
            let mut nbhd_found = false;
            
            for subst in subst_iter {
                if self.barcode_nbhds.contains_key(&subst) {
                    let nbhd = self.barcode_nbhds.get(&subst).unwrap();
                    nbhd.borrow_mut().insert(bc);
                    self.barcode_nbhds.insert(bc.to_vec(), Rc::clone(nbhd));
                    return;
                }
            }
            
            let mut nbhd = BarcodeNbhd::new();
            nbhd.insert(bc);
            let mut nbhd_rc = Rc::new(RefCell::new(nbhd));
            self.nbhds.push(Rc::clone(&nbhd_rc));
            self.barcode_nbhds.insert(bc.to_vec(), nbhd_rc);
        }
    }

    pub fn write_nbhds<W: Write>(&self, mut out: W) -> Result<(), failure::Error> {
        for nbhd in self.nbhds.iter() {
            let mut barcodes = nbhd.borrow().barcodes.clone();
            barcodes.sort_by_key(|(_, ct)| -(*ct as isize));
            let (keybc, keyct) = barcodes.first().unwrap();
            let total: usize = barcodes.iter().map(|(_, ct)| *ct).sum();
            write!(out, "{}\t{}\t{}\t{:0.3}",
                   String::from_utf8_lossy(keybc),
                   barcodes.len(), total,
                   (*keyct as f64) / (total as f64))?;
            for (bc, ct) in barcodes.iter() {
                write!(out, "\t{}\t{}",
                       String::from_utf8_lossy(bc), ct)?;
            }
            write!(out, "\n")?;
        }
        Ok(())
    }
    
    pub fn write_barcode_nbhds<W: Write>(&self, mut out: W) -> Result<(), failure::Error> {
        for (key, val) in self.barcode_nbhds.iter() {
            write!(out, "{}",
                   String::from_utf8_lossy(key))?;
            
            for (bc, ct) in val.borrow().barcodes.iter() {
                write!(out, "\t{}\t{}",
                       String::from_utf8_lossy(bc), ct)?;
            }

            write!(out, "\n")?;
        }

        Ok(())
    }
}

struct BarcodeNbhd {
    barcodes: Vec<(Vec<u8>, usize)>,
}

impl BarcodeNbhd {
    fn new() -> Self {
        BarcodeNbhd { barcodes: Vec::new() }
    }

    fn insert(&mut self, bc_new: &[u8]) -> () {
        for &mut (ref bc, ref mut ct) in self.barcodes.iter_mut() {
            if bc_new == bc.as_slice() {
                *ct += 1;
                return;
            }
        }

        self.barcodes.push((bc_new.to_vec(), 1));
    }
}

// Switch to an interface where mutations (acting on a slice buffer)
// are returned to avoid allocation.

const NTS_LEN: usize = 4;
static NTS: [u8; NTS_LEN] = [b'A', b'C', b'G', b'T'];

struct Substitutions<'a> {
    original: &'a [u8],
    position: usize,
    nt: usize,
}

impl <'a> Substitutions<'a> {
    pub fn new(original: &'a [u8]) -> Self {
        Substitutions { original: original, position: 0, nt: 0 }
    }
}

impl <'a> Iterator for Substitutions<'a> {
    type Item = Vec<u8>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.position >= self.original.len() {
            return None;
        }

        if self.nt >= NTS_LEN {
            self.position += 1;
            self.nt = 0;
            return self.next();
        }

        if self.original[self.position] == NTS[self.nt] {
            self.nt += 1;
            return self.next();
        }

        let mut variant = self.original.to_vec();
        variant[self.position] = NTS[self.nt];
        self.nt += 1;
        return Some(variant);
    }
}

struct Deletions<'a> {
    original: &'a [u8],
    position: usize,
}

impl <'a> Deletions<'a> {
    pub fn new(original: &'a [u8]) -> Self {
        Deletions { original: original, position: 0 }
    }
}

impl <'a> Iterator for Deletions<'a> {
    type Item = Vec<u8>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.position >= self.original.len() {
            return None;
        }

        let mut variant = self.original.to_vec();
        variant.remove(self.position);
        self.position += 1;
        return Some(variant);
    }
}

