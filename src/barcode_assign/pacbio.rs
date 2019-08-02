use bio::pattern_matching::myers::Myers;

#[derive(Debug,Clone,Hash,PartialEq,Eq)]
pub struct LibConf {
    frag_match: MatchConf,
    barcode_match: MatchConf,
    barcode_rev: bool,
}

impl LibConf {
    pub fn new(frag_match: MatchConf, barcode_match: MatchConf, barcode_rev: bool) -> Self {
        LibConf { frag_match: frag_match,
                  barcode_match: barcode_match,
                  barcode_rev: barcode_rev }
    }
}

#[derive(Debug,Clone,Hash,PartialEq,Eq)]
pub struct MatchConf {
    before: Vec<u8>,
    after: Vec<u8>,
    max_errors: u8
}

pub struct MatchResult<'a> {
    pub before_match: (usize, usize, u8),
    pub after_match: (usize, usize, u8),
    pub insert: &'a [u8],
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
    


