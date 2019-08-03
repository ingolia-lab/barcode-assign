use std::fmt;

use failure;

#[derive(Debug,Clone,Hash,PartialEq,Eq)]
pub struct BlasrAlign {
    query_name: String,
    target_name: String,
    query_strand: bool,
    target_strand: bool,
    score: isize,
    percent_similarity: String,
    target_start: usize,
    target_end: usize,
    target_length: usize,
    query_start: usize,
    query_end: usize,
    query_length: usize,
    n_cells: usize,
}

impl BlasrAlign {
    pub fn new(line: &str) -> Result<Self, failure::Error> {
        let fields: Vec<&str> = line.split(" ").collect();
        if fields.len() != 13 {
            return Err(format_err!("Wrong number of fields in {:?}", line));
        }

        Ok( BlasrAlign {
            query_name: fields[0].to_string(),
            target_name: fields[1].to_string(),
            query_strand: fields[2].parse::<u8>()? > 0,
            target_strand: fields[3].parse::<u8>()? > 0,
            score: fields[4].parse::<isize>()?,
            percent_similarity: fields[5].to_string(),
            target_start: fields[6].parse::<usize>()?,
            target_end: fields[7].parse::<usize>()?,
            target_length: fields[8].parse::<usize>()?,
            query_start: fields[9].parse::<usize>()?,
            query_end: fields[10].parse::<usize>()?,
            query_length: fields[11].parse::<usize>()?,
            n_cells: fields[12].parse::<usize>()?,
        })
    }

    pub fn query_name(&self) -> &str { &self.query_name }
    pub fn target_name(&self) -> &str { &self.target_name }
    pub fn query_strand(&self) -> bool { self.query_strand }
    pub fn target_strand(&self) -> bool { self.target_strand }
    pub fn score(&self) -> isize { self.score }
    pub fn percent_similarity(&self) -> &str { &self.percent_similarity }
    pub fn target_start(&self) -> usize { self.target_start }
    pub fn target_end(&self) -> usize { self.target_end }
    pub fn target_length(&self) -> usize { self.target_length }
    pub fn query_start(&self) -> usize { self.query_start }
    pub fn query_end(&self) -> usize { self.query_end }
    pub fn query_length(&self) -> usize { self.query_length }
    pub fn n_cells(&self) -> usize { self.n_cells }

    pub fn query_name_trimmed(&self) -> &str {
        &self.query_name[0..self.query_name.rfind("/").unwrap_or(self.query_name.len())]
    }
}


impl fmt::Display for BlasrAlign {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} {} {} {} {} {} {} {} {} {} {} {} {}",
               self.query_name, self.target_name, self.query_strand as u8, self.target_strand as u8,
               self.score, self.percent_similarity, self.target_start, self.target_end, self.target_length,
               self.query_start, self.query_end, self.query_length, self.n_cells)
    }
}
