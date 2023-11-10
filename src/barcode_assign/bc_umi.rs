use std::collections::HashMap;
use std::fs::File;
use std::io::{self, Read, Write};

use bio::io::fastq;

use neighborhood::*;

#[derive(Debug)]
pub struct Config {
    pub barcode_fastq: String,
    pub umi_prefix: String,
    pub out_barcodes: String,
    pub dedup_stats: Option<String>,
    pub neighborhood: Option<String>,
}

pub fn bc_umi(config: Config) -> Result<(), failure::Error> {
    let reader: Box<dyn Read> = if config.barcode_fastq == "-" {
        Box::new(io::stdin())
    } else {
        Box::new(File::open(&config.barcode_fastq)?)
    };
    let barcode_reader = fastq::Reader::new(reader);

    let mut barcode_umis = BarcodeUmis::new();

    for recres in barcode_reader.records() {
        let rec = recres?;
        let desc = rec.desc().ok_or_else(|| format_err!("No header for {:?}", rec.id()))?;
        let umi = BarcodeUmis::find_umi(&config.umi_prefix, desc)
            .ok_or_else(|| format_err!("No UMI in {:?} for {:?}", desc, rec.id()))?;
        barcode_umis.count_one(rec.seq(), umi.as_bytes());
    }

    let final_counts = if let Some(nbhd_filename) = config.neighborhood {
        neighborhood_counts(barcode_umis, &nbhd_filename)?
    } else {
        barcode_umis
    };

    if let Some(dedup_base) = config.dedup_stats {
        final_counts.write_tables(&dedup_base)?;
    }
    
    let writer: Box<dyn Write> = if config.out_barcodes == "-" {
        Box::new(io::stdout())
    } else {
        Box::new(File::create(&config.out_barcodes)?)
    };
    final_counts.write(writer)?;

    Ok(())
}

fn neighborhood_counts(
    raw_umis: BarcodeUmis,
    nbhd_filename: &str,
) -> Result<BarcodeUmis, failure::Error> {
    let nbhds_raw = Neighborhood::gather_neighborhoods(raw_umis.barcode_map());
    let nbhds: Vec<_> = nbhds_raw.into_iter().map(|n| n.into_sorted()).collect();

    let nbhd_sizes: Vec<_> = nbhds.iter().map(|n| n.with_mapped_values(|u| u.total_counts())).collect::<Vec<SortedNeighborhood<usize>>>();
    
    SortedNeighborhood::write_tables(&nbhd_filename, nbhd_sizes.iter())?;

    Ok(std::iter::FromIterator::from_iter(nbhds.into_iter().map(|n| UmiCounts::merge(n))))
}


#[derive(Debug, Clone)]
pub struct UmiCounts {
    total: usize,
    umi_counts: HashMap<Vec<u8>, usize>
}

impl UmiCounts {
    pub fn new() -> Self {
        UmiCounts { total: 0, umi_counts: HashMap::new() }
    }
    
    pub fn count_one(&mut self, umi: &[u8]) -> () {
        self.count_add(umi, 1)
            
    }

    pub fn count_add(&mut self, umi: &[u8], count: usize) -> () {
        let umi_count = if let Some(existing) = self.umi_counts.get_mut(umi) {
            existing
        } else {
            self.umi_counts.insert(umi.to_owned(), 0);
            self.umi_counts.get_mut(umi).unwrap()
        };

        *umi_count += count;
        self.total += count;
    }

    pub fn total_counts(&self) -> usize { self.total }
    pub fn total_umis(&self) -> usize { self.umi_counts.len() }

    pub fn counts(&self) -> Vec<usize> {
        let mut cts = self.umi_counts.iter().map(|uc| *uc.1).collect::<Vec<usize>>();
        cts.sort();
        cts.reverse();
        cts
    }

    fn merge(sorted_nbhd: SortedNeighborhood<UmiCounts>) -> (Vec<u8>, UmiCounts) {
        let mut barcode_umi_iter = sorted_nbhd.into_barcodes();

        let (key, mut umi_counts) = barcode_umi_iter.next().unwrap();
        for (_, more_counts) in barcode_umi_iter {
            umi_counts += more_counts;
        }
        
        (key, umi_counts)        
    }
}

impl std::ops::AddAssign for UmiCounts {
    fn add_assign(&mut self, other: Self) {
        for (umi, count) in other.umi_counts.into_iter() {
            self.count_add(&umi, count);
        }
    }
}

impl OrdEntry for UmiCounts {
    fn entry_cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.total.cmp(&other.total)
    }
}

 /// Tabulation of barcode counts in a sample
#[derive(Debug, Clone)]
pub struct BarcodeUmis(HashMap<Vec<u8>, UmiCounts>);

impl BarcodeUmis {
    pub fn new() -> Self
    {
        BarcodeUmis(HashMap::new())
    }

    pub fn barcode_map(self) -> HashMap<Vec<u8>, UmiCounts> {
        self.0
    }
    
    pub fn count_one(&mut self, barcode: &[u8], umi: &[u8]) -> ()
    {
        let barcode_umis = if let Some(existing) = self.0.get_mut(barcode) {
            existing
        } else {
            self.0.insert(barcode.to_owned(), UmiCounts::new());
            self.0.get_mut(barcode).unwrap()
        };

        barcode_umis.count_one(umi);
    }

    pub fn find_umi<'a, 'b>(prefix: &'b str, desc: &'a str) -> Option<&'a str> {
        let (_, rest) = desc.split_once(prefix)?;
        rest.split_whitespace().next()
    }

    pub fn write<W: Write>(&self, umi_out: W) -> Result<(), failure::Error> {
        let mut out = std::io::BufWriter::new(umi_out);

        for (barcode, umis) in self.0.iter() {
            write!(out, "{}\t{}\n", String::from_utf8_lossy(barcode), umis.total_umis())?;
        }
        
        Ok(())
    }
    
    pub fn write_tables(&self, filebase: &str) -> Result<(), std::io::Error> {
        // UMI deduplication statistics
        let mut umis_out_file =
            std::fs::File::create(SortedNeighborhood::output_filename(filebase, "-umi-dedup.txt"))?;
        let mut umis_out = std::io::BufWriter::new(umis_out_file);

        for (barcode, umis) in self.0.iter() {
            write!(umis_out, "{}", String::from_utf8_lossy(barcode))?;
            let umi_counts = umis.counts();

            write!(umis_out, "\t{}\t{}\t", umis.total_counts(), umis.total_umis())?;

            write!(umis_out, "{}\t",
                   umi_counts[umi_counts.len() / 2])?;
            
            for count in umi_counts.iter() {
                write!(umis_out, "{},", count)?;
            }
            write!(umis_out, "\n")?;
        }
        
        Ok(())
    }
}

impl std::iter::FromIterator<(Vec<u8>, UmiCounts)> for BarcodeUmis {
    fn from_iter<I>(iter: I) -> Self
        where I: IntoIterator<Item = (Vec<u8>, UmiCounts)>
    {
        BarcodeUmis(HashMap::from_iter(iter))
    }
}
