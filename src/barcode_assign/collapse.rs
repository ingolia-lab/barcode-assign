use std::cell::RefCell;
use std::collections::HashMap;
use std::rc::Rc;

pub struct BarcodeEquivMap {
    barcode_equivs: HashMap<Vec<u8>, Rc<RefCell<BarcodeEquiv>>>,
    max_dist: usize,
}

impl BarcodeEquivMap {
    pub fn new(max_dist: usize) -> Self {
        BarcodeEquivMap {
            barcode_equivs: HashMap::new(),
            max_dist: max_dist,
        }
    }

    pub fn insert(&mut self, bc: &[u8]) -> () {
        unimplemented!();
    }
}

struct BarcodeEquiv {
    barcodes: Vec<(Vec<u8>, usize)>,
}

impl BarcodeEquiv {
    fn new() -> Self {
        BarcodeEquiv { barcodes: Vec::new() }
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
        return Some(variant);
    }
}
