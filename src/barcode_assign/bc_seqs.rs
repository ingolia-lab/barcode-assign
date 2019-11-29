use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use std::path::{Path,PathBuf};

use bio::io::fastq;

use fastq_pair;
use neighborhood::*;

#[derive(Debug)]
pub struct Config {
    pub barcode_fastq: String,
    pub sequ_fastq: String,
    pub out_fastq: String,
    pub out_barcodes: Option<String>,
    pub out_barcode_freqs: Option<String>,
    pub neighborhood: Option<String>,
}

pub fn bc_seqs(config: Config) -> Result<(), failure::Error> {
    let barcode_reader = fastq::Reader::from_file(&config.barcode_fastq)?;
    let sequ_reader = fastq::Reader::from_file(&config.sequ_fastq)?;

    let mut fastq_writer = fastq::Writer::to_file(&config.out_fastq)?;

    let pair_records =
        fastq_pair::PairRecords::new(barcode_reader.records(), sequ_reader.records());

    if let Some(nbhd_filename) = config.neighborhood {
        let mut barcode_recs = HashMap::new();

        for pair_result in pair_records {
            let (barcode_record, sequ_record) = pair_result?;
            let barcode = barcode_record.seq().to_vec();
            let recs = barcode_recs.entry(barcode).or_insert_with(|| Vec::new());
            recs.push(sequ_record);
        }

        let mut nbhds = Neighborhood::gather_neighborhoods(barcode_recs);
        for nbhd in nbhds.iter_mut() {
            nbhd.sort_by_counts();
        }
        nbhds.sort_by_cached_key(|nbhd| nbhd.key_barcode().0.to_vec());

        let mut nbhd_counts = Vec::new();
        for nbhd in nbhds {
            let barcode = String::from_utf8(nbhd.key_barcode().0.to_vec())?;
            let mut idx = 1;

            for (_, sequ_recs) in nbhd.barcodes() {
                for sequ_rec in sequ_recs.iter() {
                    let name = format!("{}_{}", barcode, idx);
                    idx += 1;
                    let named_record =
                        fastq::Record::with_attrs(&name, None, sequ_rec.seq(), sequ_rec.qual());
                    fastq_writer.write_record(&named_record)?;
                }
            }
            
            nbhd_counts.push(nbhd.to_counts());
        }

        // Neighborhood grouping statistics
        let mut barcode_to_nbhd_out = std::fs::File::create(output_filename(&nbhd_filename, "-barcode-to-nbhd.txt"))?;
        let mut nbhd_count_out = std::fs::File::create(output_filename(&nbhd_filename, "-nbhd-count.txt"))?;
        let mut nbhds_out = std::fs::File::create(output_filename(&nbhd_filename, "-nbhds.txt"))?;

        writeln!(barcode_to_nbhd_out, "{}", Neighborhood::barcode_counts_header())?;
        writeln!(nbhds_out, "{}", Neighborhood::nbhd_counts_header())?;
        
        for nbhd in nbhd_counts {
            nbhd.write_total_counts(&mut nbhd_count_out)?;
            nbhd.write_barcode_counts(&mut barcode_to_nbhd_out)?;
            nbhd.write_nbhd_counts(&mut nbhds_out)?;
        }        
    } else {
        let mut barcode_counts = HashMap::new();
        
        for pair_result in pair_records {
            let (barcode_record, sequ_record) = pair_result?;
            
            let name = sequ_barcoded_name(&config, &mut barcode_counts, &barcode_record)?;
            
            let named_record =
                fastq::Record::with_attrs(&name, None, sequ_record.seq(), sequ_record.qual());
            
            fastq_writer.write_record(&named_record)?;
        }
        
        if let Some(barcode_filename) = config.out_barcodes {
            let barcode_writer = File::create(barcode_filename)?;
            write_barcode_table(barcode_writer, &barcode_counts)?;
        }
        
        if let Some(freq_filename) = config.out_barcode_freqs {
            let freq_writer = File::create(freq_filename)?;
            write_freq_table(freq_writer, &barcode_counts)?;
        }
    }

    Ok(())
}

fn write_barcode_table<W>(
    barcode_out: W,
    barcode_counts: &HashMap<String, usize>,
) -> Result<(), failure::Error>
where
    W: std::io::Write,
{
    let mut bcout = std::io::BufWriter::new(barcode_out);

    for (barcode, count) in barcode_counts.iter() {
        bcout.write(barcode.as_bytes())?;
        bcout.write("\t".as_bytes())?;
        bcout.write(count.to_string().as_bytes())?;
        bcout.write("\n".as_bytes())?;
    }

    Ok(())
}

fn write_freq_table<W>(
    freq_out: W,
    barcode_counts: &HashMap<String, usize>,
) -> Result<(), failure::Error>
where
    W: std::io::Write,
{
    let mut fout = std::io::BufWriter::new(freq_out);

    let mut freq_counts = HashMap::new();

    for freq in barcode_counts.values() {
        let freq_count = freq_counts.entry(freq).or_insert(0);
        *freq_count += 1;
    }

    let mut freqs: Vec<usize> = freq_counts.keys().map(|&&k| k).collect();
    freqs.sort();

    for freq in freqs {
        write!(fout, "{}\t{}\n", freq, freq_counts.get(&freq).unwrap_or(&0))?;
    }

    Ok(())
}

fn sequ_barcoded_name(
    _config: &Config,
    barcode_counts: &mut HashMap<String, usize>,
    barcode_record: &fastq::Record,
) -> Result<String, failure::Error> {
    let barcode = String::from_utf8(barcode_record.seq().to_vec())?;
    let barcode_count = barcode_counts.entry(barcode.to_string()).or_insert(0);
    *barcode_count += 1;

    Ok(barcode + "_" + &barcode_count.to_string())
}

fn output_filename(output_base: &str, name: &str) -> PathBuf {
    let base_ref: &Path = output_base.as_ref();
    let mut namebase = base_ref
        .file_name()
        .map_or(std::ffi::OsString::new(), std::ffi::OsStr::to_os_string);
    namebase.push(name);
    base_ref.with_file_name(namebase)
}
