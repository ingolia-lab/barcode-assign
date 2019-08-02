use bio::alphabets::dna;
use bio::pattern_matching::myers::Myers;

#[derive(Debug,Clone,Hash,PartialEq,Eq)]
pub struct LibConf {
    frag_matcher: MatchConf,
    barcode_matcher: MatchConf,
    barcode_rev: bool,
}

impl LibConf {
    pub fn new(frag_matcher: MatchConf, barcode_matcher: MatchConf, barcode_rev: bool) -> Self {
        LibConf { frag_matcher: frag_matcher,
                  barcode_matcher: barcode_matcher,
                  barcode_rev: barcode_rev }
    }

    pub fn best_match<'a>(&self, target: &'a [u8]) -> Option<LibMatch<'a>> {
        if let Some(barcode_match) = self.barcode_matcher.best_match(target) {
            if let Some(frag_match) = self.frag_matcher.best_match(target) {
                Some(LibMatch{frag_match: frag_match,
                                    barcode_match: barcode_match,
                                    barcode_rev: self.barcode_rev})
            } else {
                None
            }
        } else {
            None
        }
    }
}

#[derive(Debug,Clone,Hash,PartialEq,Eq)]
pub struct LibMatch<'a> {
    frag_match: MatchResult<'a>,
    barcode_match: MatchResult<'a>,
    barcode_rev: bool,
}
    
impl <'a> LibMatch<'a> {
    pub fn frag_match(&self) -> &MatchResult<'a> { &self.frag_match }
    pub fn barcode_match(&self) -> &MatchResult<'a> { &self.barcode_match }
    pub fn barcode_rev(&self) -> bool { self.barcode_rev }

    pub fn barcode_actual(&self) -> String {
        if self.barcode_rev {
            String::from_utf8_lossy(&dna::revcomp(self.barcode_match.insert)).to_string()
        } else {
            String::from_utf8_lossy(self.barcode_match.insert).to_string()
        }
    }
}

#[derive(Debug,Clone,Hash,PartialEq,Eq)]
pub struct MatchConf {
    before: Vec<u8>,
    after: Vec<u8>,
    max_errors: u8
}

impl MatchConf {
    pub fn new(before: &[u8], after: &[u8], max_errors: u8) -> Self {
        MatchConf { before: before.to_vec(), after: after.to_vec(), max_errors: max_errors }
    }

    pub fn best_match<'a>(&self, target: &'a [u8]) -> Option<MatchResult<'a>> {
        let mut myers_before = Myers::<u64>::new(&self.before);
        let mut myers_after = Myers::<u64>::new(&self.after);

        // N.B. end coordinate is not included in match
        let hit_before = myers_before.find_all(target, self.max_errors).by_ref().min_by_key(|&(_, _, dist)| dist);
        if let Some((before_start, before_end, before_dist)) = hit_before {
            let rest = target.split_at(before_end).1;
            
            let hit_after = myers_after.find_all(rest, self.max_errors - before_dist).by_ref().min_by_key(|&(_, _, dist)| dist);
            if let Some((rest_after_start, rest_after_end, after_dist)) = hit_after {
                Some(MatchResult{before_match: (before_start, before_end, before_dist),
                                 after_match: (before_end + rest_after_start, before_end + rest_after_end, after_dist),
                                 insert: rest.split_at(rest_after_start).0})
            } else {
                None
            }
        } else {
            None
        }
    }
}

#[derive(Debug,Clone,Hash,PartialEq,Eq)]
pub struct MatchResult<'a> {
    pub before_match: (usize, usize, u8),
    pub after_match: (usize, usize, u8),
    pub insert: &'a [u8],
}

impl <'a> MatchResult<'a> {
    pub fn score(&self) -> u8 { self.before_match.2 + self.after_match.2 }

    pub fn insert_start(&self) -> usize { self.before_match.1 }
    pub fn insert_end(&self) -> usize { self.after_match.0 }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn perfect_match() {
        let match_one = MatchConf::new(b"ACGTACGT", b"CAGTCAGT", 0);

        //              01234567890123456789012345678
        let query_a = b"CTAACGTACGTTTAATTAACAGTCAGTAA";
        let res = match_one.best_match(query_a);
        assert_eq!(res, Some(MatchResult{before_match:(3, 11, 0),
                                         after_match: (19, 27, 0),
                                         insert: &query_a[11..19]}));
    }
}
