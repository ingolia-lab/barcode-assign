use std::fs::File;
use std::io::{self, Read, Write};
use std::path::{Path,PathBuf};

use bio::io::fastq;

use neighborhood::Neighborhood;
use counts::SampleCounts;

#[derive(Debug)]
pub struct Config {
    pub barcode_fastq: String,
    pub out_barcodes: String,
    pub freq_filename: Option<String>,
    pub neighborhood: Option<String>,
}

pub fn bc_count(config: Config) -> Result<(), failure::Error> {
    let reader: Box<Read> = if config.barcode_fastq == "-" {
        Box::new(io::stdin())
    } else {
        Box::new(File::open(&config.barcode_fastq)?)
    };
    let barcode_reader = fastq::Reader::new(reader);

    let barcode_counts_res: Result<SampleCounts, std::io::Error>
        = barcode_reader.records().collect();
    let barcode_counts = barcode_counts_res?;
    
    let writer: Box<Write> = if config.out_barcodes == "-" {
        Box::new(io::stdout())
    } else {
        Box::new(File::create(&config.out_barcodes)?)
    };
    barcode_counts.write(writer)?;
    
    if let Some(freq_filename) = config.freq_filename {
        barcode_counts.write_freq_table(File::create(freq_filename)?)?;
    }

    if let Some(nbhd_filename) = config.neighborhood {
        neighborhood_counts(barcode_counts, &nbhd_filename)?;
    }
    
    Ok(())
}

fn neighborhood_counts(barcode_counts: SampleCounts, nbhd_filename: &str) -> Result<(), failure::Error>
{
    let mut barcode_to_nbhd_out = std::fs::File::create(output_filename(nbhd_filename, "-barcode-to-nbhd.txt"))?;
    let mut nbhd_count_out = std::fs::File::create(output_filename(nbhd_filename, "-nbhd-count.txt"))?;
    let mut nbhds_out = std::fs::File::create(output_filename(nbhd_filename, "-nbhds.txt"))?;

    let mut nbhds = Neighborhood::gather_neighborhoods(barcode_counts.count_map());
    for nbhd in nbhds.iter_mut() {
        nbhd.sort_by_counts();
    }

    writeln!(barcode_to_nbhd_out, "{}", Neighborhood::barcode_counts_header())?;
    writeln!(nbhds_out, "{}", Neighborhood::nbhd_counts_header())?;
    
    for nbhd in nbhds.iter() {
        nbhd.write_total_counts(&mut nbhd_count_out)?;
        nbhd.write_barcode_counts(&mut barcode_to_nbhd_out)?;
        nbhd.write_nbhd_counts(&mut nbhds_out)?;
    }

    Ok(())
}

fn output_filename(output_base: &str, name: &str) -> PathBuf {
    let base_ref: &Path = output_base.as_ref();
    let mut namebase = base_ref
        .file_name()
        .map_or(std::ffi::OsString::new(), std::ffi::OsStr::to_os_string);
    namebase.push(name);
    base_ref.with_file_name(namebase)
}

#[cfg(test)]
mod tests {
    extern crate tempfile;
    extern crate rand;
    
    use self::rand::thread_rng;
    use self::rand::Rng;
    
    use super::*;

    #[test]
    fn count_literal() {
        let barcode_fq = r#"@one
ACGTTGCA
+
~~~~~~~~
@two
CGTAATGC
+
~~~~~~~~
@three
ACGTTGCA
+
~~~~~~~~
@four
GTACCATG
+
~~~~~~~~
@five
CGTAATGC
+
~~~~~~~~
@six
ACGTTGCA
+
~~~~~~~~
"#;

        let mut fastq_file = tempfile::NamedTempFile::new().unwrap();
        fastq_file.write_all(barcode_fq.as_bytes()).unwrap();

        let count_path = tempfile::NamedTempFile::new().unwrap().into_temp_path();
        
        let config = Config { barcode_fastq: fastq_file.path().to_string_lossy().into_owned(),
                              out_barcodes: count_path.to_string_lossy().into_owned(),
                              freq_filename: None,
                              neighborhood: None
        };

        bc_count(config).unwrap();

        let counts = SampleCounts::from_file(count_path).unwrap();
        let mut cvec: Vec<(Vec<u8>, usize)> = counts.into_iter().collect();
        cvec.sort();
        assert_eq!(cvec, vec![(b"ACGTTGCA".to_vec(), 3),
                              (b"CGTAATGC".to_vec(), 2),
                              (b"GTACCATG".to_vec(), 1)]);
    }

    fn barcode_records(bc: &[u8], ct: usize) -> impl Iterator<Item = fastq::Record> {
        let barcode = bc.to_vec();
        let qual = vec![b'~'; barcode.len()];
        std::iter::repeat_with(move || fastq::Record::with_attrs("barcode", None, &barcode, &qual)).take(ct)
    }
    
    #[test]
    fn count_many() {
        let mut counts = vec![(b"ACGTACGT".to_vec(), 5),
                              (b"ACGTTCGT".to_vec(), 3),
                              (b"ACATACGT".to_vec(), 2),
                              (b"CGTACGTA".to_vec(), 8),
                              (b"CGTACGAA".to_vec(), 4),
                              (b"GTACGTACG".to_vec(), 7),
                              (b"GTACGTCG".to_vec(), 1),
                              (b"GTACGCACG".to_vec(), 9),
                              (b"GTACGCATCG".to_vec(), 6)];

        let mut records: Vec<fastq::Record> = counts.iter().flat_map(|(bc, ct)| barcode_records(bc, *ct)).collect();
        thread_rng().shuffle(&mut records);
        
        let fastq_file = tempfile::NamedTempFile::new().unwrap();
        {
            let mut fastq_writer = fastq::Writer::new(&fastq_file);
            for rec in records.iter() {
                fastq_writer.write_record(rec).unwrap();
            }
        }
        let fastq_path = fastq_file.into_temp_path();
        
        let count_path = tempfile::NamedTempFile::new().unwrap().into_temp_path();
        
        let config = Config { barcode_fastq: fastq_path.to_string_lossy().into_owned(),
                              out_barcodes: count_path.to_string_lossy().into_owned(),
                              freq_filename: None,
                              neighborhood: None
        };

        bc_count(config).unwrap();

        let out_counts = SampleCounts::from_file(count_path).unwrap();
        let mut cvec: Vec<(Vec<u8>, usize)> = out_counts.count_map().into_iter().collect();
        cvec.sort();
        counts.sort();
        assert_eq!(counts, cvec);
    }

    #[test]
    fn count_frequencies() {
        let mut counts = vec![(b"ACGTACGT".to_vec(), 5),
                              (b"ACGTTCGT".to_vec(), 3),
                              (b"ACATACGT".to_vec(), 2),
                              (b"CGTACGTA".to_vec(), 5),
                              (b"CGTACGAA".to_vec(), 1),
                              (b"GTACGTACG".to_vec(), 4),
                              (b"GTACGTCG".to_vec(), 3),
                              (b"GTACGCACG".to_vec(), 7),
                              (b"GTACGCATCG".to_vec(), 5)];

        let mut records: Vec<fastq::Record> = counts.iter().flat_map(|(bc, ct)| barcode_records(bc, *ct)).collect();
        thread_rng().shuffle(&mut records);
        
        let fastq_file = tempfile::NamedTempFile::new().unwrap();
        {
            let mut fastq_writer = fastq::Writer::new(&fastq_file);
            for rec in records.iter() {
                fastq_writer.write_record(rec).unwrap();
            }
        }
        let fastq_path = fastq_file.into_temp_path();
        
        let count_path = tempfile::NamedTempFile::new().unwrap().into_temp_path();
        let freq_path = tempfile::NamedTempFile::new().unwrap().into_temp_path();
        
        let config = Config { barcode_fastq: fastq_path.to_string_lossy().into_owned(),
                              out_barcodes: count_path.to_string_lossy().into_owned(),
                              freq_filename: Some(freq_path.to_string_lossy().into_owned()),
                              neighborhood: None
        };

        bc_count(config).unwrap();

        let out_counts = SampleCounts::from_file(count_path).unwrap();
        let mut cvec: Vec<(Vec<u8>, usize)> = out_counts.count_map().into_iter().collect();
        cvec.sort();
        counts.sort();
        assert_eq!(counts, cvec);

        let out_freq = std::fs::read(freq_path).unwrap();
        assert_eq!(String::from_utf8_lossy(&out_freq), r#"1	1
2	1
3	2
4	1
5	3
7	1
"#);
    }

    #[test]
    fn count_nbhd() {
        let counts = vec![(b"ACGTACGT".to_vec(), 5),
                          (b"ACGTTCGT".to_vec(), 3),
                          (b"ACATACGT".to_vec(), 2),
                          (b"CGTACGTA".to_vec(), 8),
                          (b"CGTACGAA".to_vec(), 4),
                          (b"GTACGTACG".to_vec(), 7),
                          (b"GTACGTCG".to_vec(), 1),
                          (b"GTACGCACG".to_vec(), 9),
                          (b"GTACGCATCG".to_vec(), 6)];
        
        let mut nbhds = vec![(b"ACGTACGT".to_vec(), 5 + 3 + 2),
                             (b"CGTACGTA".to_vec(), 8 + 4),
                             (b"GTACGCACG".to_vec(), 7 + 1 + 9 + 6)];
        
        let mut records: Vec<fastq::Record> = counts.iter().flat_map(|(bc, ct)| barcode_records(bc, *ct)).collect();
        thread_rng().shuffle(&mut records);
        
        let fastq_file = tempfile::NamedTempFile::new().unwrap();
        {
            let mut fastq_writer = fastq::Writer::new(&fastq_file);
            for rec in records.iter() {
                fastq_writer.write_record(rec).unwrap();
            }
        }
        let fastq_path = fastq_file.into_temp_path();
        
        let count_path = tempfile::NamedTempFile::new().unwrap().into_temp_path();
        let nbhd_path = tempfile::NamedTempFile::new().unwrap().into_temp_path();
        
        let config = Config { barcode_fastq: fastq_path.to_string_lossy().into_owned(),
                              out_barcodes: count_path.to_string_lossy().into_owned(),
                              freq_filename: None,
                              neighborhood: Some(nbhd_path.to_string_lossy().into_owned()),
        };

        bc_count(config).unwrap();

        let mut nbhd_count_filename = nbhd_path.to_string_lossy().into_owned();
        nbhd_count_filename += "-nbhd-count.txt";
        
        let out_counts = SampleCounts::from_file(&nbhd_count_filename).unwrap();
        let mut cvec: Vec<(Vec<u8>, usize)> = out_counts.count_map().into_iter().collect();
        cvec.sort();
        nbhds.sort();
        assert_eq!(nbhds, cvec);
    }
}
