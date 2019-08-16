use std::collections::HashMap;
use std::io::Write;

pub struct Neighborhood<T> {
    barcodes: Vec<(Vec<u8>, T)>
}

impl <T> Neighborhood<T> {
    fn new() -> Self {
        Neighborhood {
            barcodes: Vec::new()
        }
    }
    
    fn insert(&mut self, barcode: Vec<u8>, value: T) -> () {
        self.barcodes.push((barcode, value));
    }

    pub fn barcodes(&self) -> impl Iterator<Item = &(Vec<u8>, T)> {
        self.barcodes.iter()
    }

    pub fn len(&self) -> usize { self.barcodes.len() }

    pub fn key_barcode(&self) -> (&[u8], &T) {
        let (keybc, keyct) = self.barcodes.first().unwrap();
        (keybc, keyct)
    }
    
    // Collecting a neighborhood:
    // 1. Pick a node arbitrarily
    //    a. remove from key set
    //    b. initialize work stack with node
    // 2. Handle a node from work stack
    //    a. check key set for all near-neighbors
    //       i. remove near-neighbor from key set
    //       ii. push near-neighbor onto work stack
    //    b. add node to neighborhood
    // 3. Repeat handling nodes from work stack until empty
    
    pub fn gather_neighborhoods(mut bc_map: HashMap<Vec<u8>, T>) -> Vec<Neighborhood<T>>
    {
        let mut neighborhoods = Vec::new();
        
        loop {
            let mut work_stack = Vec::new();
            let mut neighborhood = Neighborhood::new();
            
            let start_ref = match bc_map.iter().next() {
                Some(start_ref) => start_ref,
                None => { break; }
            };
            
            let start = start_ref.0.to_vec();
            let value = bc_map.remove(&start).unwrap();
            work_stack.push((start, value));
            
            while work_stack.len() > 0 {
                let (curr, curr_value) = work_stack.pop().unwrap();
                
                let neighbors = Substitutions::new(&curr).chain(Deletions::new(&curr)).chain(Insertions::new(&curr));
                for neighbor in neighbors {
                    if bc_map.contains_key(&neighbor) {
                        let neighbor_value = bc_map.remove(&neighbor).unwrap();
                        work_stack.push((neighbor, neighbor_value));
                    }
                }
                
                neighborhood.insert(curr, curr_value);
            }
            
            neighborhoods.push(neighborhood);
        }
        
        neighborhoods
    }
}

impl Neighborhood<usize> {
    pub fn sort_by_counts(&mut self) -> () {
        self.barcodes.sort_by_key(|(_, ct)| -(*ct as isize));
    }

    pub fn total(&self) -> usize {
        self.barcodes().map(|(_, ct)| *ct).sum()
    }

    pub fn write_total_counts<W: Write>(&self, out: &mut W) -> Result<(), std::io::Error> {
        write!(out, "{}\t{}\n", String::from_utf8_lossy(self.key_barcode().0), self.total())
    }
    
    pub fn barcode_counts_header() -> &'static str { "barcode\tneighborhood\tcount\ttotal\tfraction" }

    pub fn write_barcode_counts<W: Write>(&self, out: &mut W) -> Result<(), std::io::Error> {
        let (keybc, _keyct) = self.key_barcode();
        let total = self.total();
        
        for (bc, ct) in self.barcodes() {
            write!(out, "{}\t{}\t{}\t{}\t{:0.3}\n",
                   String::from_utf8_lossy(bc),
                   String::from_utf8_lossy(keybc),
                   ct, total, (*ct as f64) / (total as f64))?;
        }

        Ok(())
    }

    pub fn nbhd_counts_header() -> &'static str { "neighborhood\tnum_barcodes\ttotal\tfract_nbhd" }
    
    pub fn write_nbhd_counts<W: Write>(&self, out: &mut W) -> Result<(), std::io::Error> {
        let (keybc, keyct) = self.key_barcode();
        let total = self.total();

        write!(out, "{}\t{}\t{}\t{:0.3}",
               String::from_utf8_lossy(keybc),
               self.len(), total,
               (*keyct as f64) / (total as f64))?;
        for (bc, ct) in self.barcodes() {
            write!(out, "\t{}\t{}",
                   String::from_utf8_lossy(bc), ct)?;
        }
        write!(out, "\n")        
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

struct Insertions<'a> {
    original: &'a [u8],
    position: usize,
    nt: usize,
}

impl <'a> Insertions<'a> {
    pub fn new(original: &'a [u8]) -> Self {
        Insertions { original: original, position: 0, nt: 0 }
    }
}

impl <'a> Iterator for Insertions<'a> {
    type Item = Vec<u8>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.position > self.original.len() {
            return None;
        }

        if self.nt >= NTS_LEN {
            self.position += 1;
            self.nt = 0;
            return self.next();
        }

        let mut variant = Vec::with_capacity(self.original.len() + 1);
        variant.extend_from_slice(&self.original[..self.position]);
        variant.push(NTS[self.nt]);
        variant.extend_from_slice(&self.original[self.position..]);
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

        let mut variant = Vec::with_capacity(self.original.len() - 1);
        variant.extend_from_slice(&self.original[..self.position]);
        variant.extend_from_slice(&self.original[(self.position+1)..]);
        self.position += 1;
        return Some(variant);
    }
}

