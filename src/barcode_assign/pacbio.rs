use bio::alphabets::dna;
use bio::pattern_matching::myers::Myers;

#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub struct LibSpec {
    frag_matcher: FlankMatchSpec,
    barcode_matcher: FlankMatchSpec,
    barcode_rev: bool,
}

impl LibSpec {
    pub fn new(
        frag_matcher: FlankMatchSpec,
        barcode_matcher: FlankMatchSpec,
        barcode_rev: bool,
    ) -> Self {
        LibSpec {
            frag_matcher: frag_matcher,
            barcode_matcher: barcode_matcher,
            barcode_rev: barcode_rev,
        }
    }

    pub fn best_match<'a>(&self, query: &'a [u8]) -> LibMatchOut<'a> {
        LibMatchOut {
            frag: self.frag_matcher.best_match(query),
            barcode: self.barcode_matcher.best_match(query),
            barcode_rev: self.barcode_rev,
        }
    }
}

#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub struct LibMatchOut<'a> {
    frag: FlankMatchOut<'a>,
    barcode: FlankMatchOut<'a>,
    barcode_rev: bool,
}

impl<'a> LibMatchOut<'a> {
    pub fn frag_match(&self) -> &FlankMatchOut<'a> {
        &self.frag
    }
    pub fn barcode_match(&self) -> &FlankMatchOut<'a> {
        &self.barcode
    }
    pub fn barcode_rev(&self) -> bool {
        self.barcode_rev
    }

    pub fn lib_match(&self) -> Option<LibMatch<'a>> {
        self.frag.flank_match().and_then(|frag| {
            self.barcode.flank_match().map(|barcode| LibMatch {
                frag: frag,
                barcode: barcode,
                barcode_rev: self.barcode_rev,
            })
        })
    }
}

#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub struct LibMatch<'a> {
    frag: FlankMatch<'a>,
    barcode: FlankMatch<'a>,
    barcode_rev: bool,
}

impl<'a> LibMatch<'a> {
    pub fn frag_match(&self) -> &FlankMatch<'a> {
        &self.frag
    }
    pub fn barcode_match(&self) -> &FlankMatch<'a> {
        &self.barcode
    }
    pub fn barcode_rev(&self) -> bool {
        self.barcode_rev
    }

    pub fn between(&self) -> isize {
        if self.frag.insert_start() > self.barcode.insert_end() {
            self.frag.insert_start() as isize - self.barcode.insert_end() as isize
        } else {
            self.barcode.insert_start() as isize - self.frag.insert_end() as isize
        }
    }

    pub fn barcode_actual(&self) -> String {
        if self.barcode_rev {
            String::from_utf8_lossy(&dna::revcomp(self.barcode.insert_seq())).to_string()
        } else {
            String::from_utf8_lossy(self.barcode.insert_seq()).to_string()
        }
    }
}

#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub struct FlankMatchSpec {
    before: Vec<u8>,
    after: Vec<u8>,
    max_errors: u8,
}

impl FlankMatchSpec {
    pub fn new(before: &[u8], after: &[u8], max_errors: u8) -> Self {
        FlankMatchSpec {
            before: before.to_vec(),
            after: after.to_vec(),
            max_errors: max_errors,
        }
    }

    pub fn best_match<'a>(&self, query: &'a [u8]) -> FlankMatchOut<'a> {
        let mut myers_before = Myers::<u64>::new(&self.before);
        let mut myers_after = Myers::<u64>::new(&self.after);

        // N.B. end coordinate is not included in match
        let best_before = myers_before
            .find_all(query, self.max_errors)
            .by_ref()
            .min_by_key(|&(_, _, dist)| dist);
        let best_after = myers_after
            .find_all(query, self.max_errors)
            .by_ref()
            .min_by_key(|&(_, _, dist)| dist);

        FlankMatchOut {
            before: best_before,
            after: best_after,
            query: query,
        }
    }
}

/// Attempted alignment between flanking constant sequences and a
/// target query. Either of the flanking sequence matches may fail.
#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub struct FlankMatchOut<'a> {
    before: Option<(usize, usize, u8)>,
    after: Option<(usize, usize, u8)>,
    query: &'a [u8],
}

impl<'a> FlankMatchOut<'a> {
    /// Returns the successful match, or `None` if either the before
    /// or after match failed.
    pub fn flank_match(&self) -> Option<FlankMatch<'a>> {
        if let Some(before) = self.before {
            if let Some(after) = self.after {
                if before.1 <= after.0 {
                    Some(FlankMatch {
                        before: before,
                        after: after,
                        query: self.query,
                    })
                } else {
                    None
                }
            } else {
                None
            }
        } else {
            None
        }
    }
}

/// Successful alignment between flanking constant sequences and a
/// target query.
#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub struct FlankMatch<'a> {
    before: (usize, usize, u8),
    after: (usize, usize, u8),
    query: &'a [u8],
}

impl<'a> FlankMatch<'a> {
    /// Returns the total score (edit distance) of the alignments
    /// between the before and after flanking sequences and the query.
    pub fn score(&self) -> u8 {
        self.before.2 + self.after.2
    }

    /// Returns the starting position of the insert sequence in the
    /// query.
    pub fn insert_start(&self) -> usize {
        self.before.1
    }

    /// Returns the ending position of the insert sequence in the
    /// query, *non*-inclusive.
    pub fn insert_end(&self) -> usize {
        self.after.0
    }

    /// Returns the insert sequence.
    pub fn insert_seq(&self) -> &[u8] {
        &self.query[self.before.1..self.after.0]
    }

    /// Returns the query sequence that matched the constant sequence
    /// before the insert.
    pub fn before_seq(&self) -> &[u8] {
        &self.query[self.before.0..self.before.1]
    }

    /// Returns the query sequence that matched the constant sequence
    /// after the insert.
    pub fn after_seq(&self) -> &[u8] {
        &self.query[self.after.0..self.after.1]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn perfect_match() {
        let match_one = FlankMatchSpec::new(b"ACGTACGT", b"CAGTCAGT", 0);

        //              01234567890123456789012345678
        let query_a = b"CTAACGTACGTTTAATTAACAGTCAGTAA";
        let m = match_one.best_match(query_a).flank_match().unwrap();
        assert_eq!(m.score(), 0);
        assert_eq!(m.insert_start(), 11);
        assert_eq!(m.insert_end(), 19);
        assert_eq!(m.insert_seq(), &query_a[11..19]);
    }
}
