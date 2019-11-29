use std::collections::HashMap;
use std::io::Write;

#[derive(Debug)]
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

impl <T> Neighborhood<Vec<T>> {
    pub fn sort_by_counts(&mut self) -> () {
        fn cmp_entries<T>((bcl, entl): &(Vec<u8>, Vec<T>), (bcr, entr): &(Vec<u8>, Vec<T>)) -> std::cmp::Ordering {
            match entl.len().cmp(&entr.len()) {
                std::cmp::Ordering::Less => std::cmp::Ordering::Greater,
                std::cmp::Ordering::Greater => std::cmp::Ordering::Less,
                std::cmp::Ordering::Equal => {
                    bcl.cmp(&bcr)
                }
            }
        }
        self.barcodes.sort_unstable_by(cmp_entries);
    }

    pub fn to_counts(&self) -> Neighborhood<usize> {
        let mut count_nbhd = Neighborhood::new();
        for (bc, ents) in self.barcodes.iter() {
            count_nbhd.insert(bc.to_vec(), ents.len());
        }
        count_nbhd
    }
}

impl Neighborhood<usize> {
    pub fn sort_by_counts(&mut self) -> () {
        fn cmp_entries((bcl, ctl): &(Vec<u8>, usize), (bcr, ctr): &(Vec<u8>, usize)) -> std::cmp::Ordering {
            match ctl.cmp(&ctr) {
                std::cmp::Ordering::Less => std::cmp::Ordering::Greater,
                std::cmp::Ordering::Greater => std::cmp::Ordering::Less,
                std::cmp::Ordering::Equal => {
                    bcl.cmp(&bcr)
                }
            }
        }
        self.barcodes.sort_unstable_by(cmp_entries);
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::borrow::Borrow;
    use counts::*;
    
    fn vec_count_map<V: AsRef<[u8]>, I: IntoIterator<Item = (V, usize)>>(bc_counts: I) -> HashMap<Vec<u8>, usize> {
        let mut ctmap = HashMap::new();
        for (bc, ct) in bc_counts {
            if ctmap.insert(bc.as_ref().to_vec(), ct).is_some() {
                panic!("Duplicate barcode {}",
                       String::from_utf8_lossy(bc.as_ref()));
            }
        }
        ctmap
    }

    fn nbhd_map<T: Clone, N: Borrow<Neighborhood<T>>>(nbhd: N) -> HashMap<Vec<u8>, T> {
        nbhd.borrow().barcodes().map(|pair| pair.clone()).collect()
    }
    
    #[test]
    fn single_nbhd() {
        let count_vec = vec![(b"ACGTACGT", 5),
                             (b"ACGTTCGT", 3),
                             (b"ACATACGT", 2),
                             (b"ACATATGT", 1)];
        
        let nbhds = Neighborhood::gather_neighborhoods(vec_count_map(count_vec.clone()));
        assert_eq!(nbhds.len(), 1);
        let exp_nbhd: HashMap<Vec<u8>, usize> = vec_count_map(count_vec.clone());
        let act_nbhd: HashMap<Vec<u8>, usize> = nbhd_map(&nbhds[0]);
        assert_eq!(exp_nbhd, act_nbhd);
    }

    #[test]
    fn two_nbhds() {
        let count_vec_1 = vec![(b"ACGTACGT", 5),
                               (b"ACGTTCGT", 3),
                               (b"ACATACGT", 2)];
        let count_vec_2 = vec![(b"CGTACGTA", 8),
                               (b"CGTACGAA", 3)];

        let mut count_vec = count_vec_1.clone();
        count_vec.extend(count_vec_2.clone());

        let nbhds = Neighborhood::gather_neighborhoods(vec_count_map(count_vec));
        assert_eq!(nbhds.len(), 2);
        let exp_nbhd_1 = vec_count_map(count_vec_1.clone());
        let exp_nbhd_2 = vec_count_map(count_vec_2.clone());

        if nbhd_map(&nbhds[0]) == exp_nbhd_1 {
            assert_eq!(nbhd_map(&nbhds[1]), exp_nbhd_2);
        } else {
            assert_eq!(nbhd_map(&nbhds[0]), exp_nbhd_2);
            assert_eq!(nbhd_map(&nbhds[1]), exp_nbhd_1);
        }
    }

    #[test]
    fn three_nbhds() {
        let count_table = r#"ACGTACGT	5
ACGTTCGT	3
ACATACGT	2
CGTACGTA	8
CGTACGAA	4
GTACGTACG	7
GTACGTCG	1
GTACGCACG	9
GTACGCATCG	6"#;
        let count_map = SampleCounts::read(count_table.as_bytes()).unwrap().count_map();
        
        let nbhds = Neighborhood::gather_neighborhoods(count_map);
        assert_eq!(nbhds.len(), 3);

        let mut exp_a = vec![b"ACGTACGT".to_vec(),
                             b"ACGTTCGT".to_vec(),
                             b"ACATACGT".to_vec()];
        exp_a.sort();
        let mut exp_c = vec![b"CGTACGTA".to_vec(),
                             b"CGTACGAA".to_vec()];
        exp_c.sort();
        let mut exp_g = vec![b"GTACGTACG".to_vec(),
                             b"GTACGTCG".to_vec(),
                             b"GTACGCACG".to_vec(),
                             b"GTACGCATCG".to_vec()];
        exp_g.sort();
        let mut exp = vec![exp_a, exp_c, exp_g];
        exp.sort();

        let mut act: Vec<Vec<Vec<u8>>> = nbhds.iter().map(|n| {
            let mut n_act: Vec<Vec<u8>> = n.barcodes().map(|(bc, _ct)| bc.clone()).collect();
            n_act.sort();
            n_act }).collect();
        act.sort();

        assert_eq!(act, exp);
    }
}
