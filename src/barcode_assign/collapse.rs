use std::cell::RefCell;
use std::collections::HashMap;
use std::io::Write;
use std::rc::Rc;

use bio::pattern_matching::myers::Myers;

// New barcode -- check first for exact match, then scan all existing barcodes
//   for near-matches
// Insert later barcodes with links to earlier barcodes?
//   Or, create barcode equivalency sets?

pub struct BarcodeNbhdMap {
    barcode_nbhds: HashMap<Vec<u8>, Rc<RefCell<BarcodeNbhd>>>,
    max_dist: u8,
}

impl BarcodeNbhdMap {
    pub fn new(max_dist: usize) -> Self {
        BarcodeNbhdMap {
            barcode_nbhds: HashMap::new(),
            max_dist: max_dist as u8,
        }
    }

    pub fn insert(&mut self, bc: &[u8]) -> () {
        if self.barcode_nbhds.contains_key(bc) {
            self.barcode_nbhds.get(bc).unwrap().borrow_mut().insert(bc);
        } else {
            // let myers = Myers::<u64>::new(bc);

            // for key in self.barcode_nbhds.keys() {
            //     let (_, dist) = myers.find_best_end(key);
            //     if dist <= self.max_dist {
            //         println!("Found {} vs {} at {}",
            //                  String::from_utf8_lossy(bc),
            //                  String::from_utf8_lossy(key),
            //                  dist);
            //     }
            // }

            let mut subst_iter = Substitutions::new(bc).chain(Deletions::new(bc));
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
            self.barcode_nbhds.insert(bc.to_vec(), Rc::new(RefCell::new(nbhd)));
        }
    }

    pub fn write_nbhds<W: Write>(&self, mut out: W) -> Result<(), failure::Error> {
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

