use std::collections::HashMap;
use std::fs::File;
use std::io::Write;

use bio::io::fastq;

use fastq_pair;

#[derive(Debug)]
pub struct Config {
    pub barcode_fastq: String,
    pub sequ_fastq: String,
    pub out_fastq: String,
    pub out_barcodes: Option<String>,
    pub out_barcode_freqs: Option<String>,
}

pub fn bc_seqs(config: Config) -> Result<(), failure::Error> {
    let barcode_reader = fastq::Reader::from_file(&config.barcode_fastq)?;
    let sequ_reader = fastq::Reader::from_file(&config.sequ_fastq)?;

    let mut fastq_writer = fastq::Writer::to_file(&config.out_fastq)?;

    let mut barcode_counts = HashMap::new();

    let pair_records = fastq_pair::PairRecords::new(barcode_reader.records(), sequ_reader.records());

    for pair_result in pair_records {
        let (barcode_record, sequ_record) = pair_result?;

        let name = sequ_barcoded_name(&config, &mut barcode_counts, &barcode_record)?;

        let named_record = fastq::Record::with_attrs(&name, None, sequ_record.seq(), sequ_record.qual());

        fastq_writer.write_record(&named_record)?;
    }

    if let Some(barcode_filename) = config.out_barcodes {
        let mut barcode_writer = File::create(barcode_filename)?;
        write_barcode_table(barcode_writer, &barcode_counts)?;
    }

    if let Some(freq_filename) = config.out_barcode_freqs {
        let mut freq_writer = File::create(freq_filename)?;
        write_freq_table(freq_writer, &barcode_counts)?;
    }
    
    Ok(())
}

fn write_barcode_table<W>(barcode_out: W, barcode_counts: &HashMap<String, usize>) -> Result<(), failure::Error>
    where W: std::io::Write
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

fn write_freq_table<W>(freq_out: W, barcode_counts: &HashMap<String, usize>) -> Result<(), failure::Error>
    where W: std::io::Write
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

fn sequ_barcoded_name(_config: &Config, barcode_counts: &mut HashMap<String, usize>, barcode_record: &fastq::Record) -> Result<String, failure::Error>
{
    let barcode = String::from_utf8(barcode_record.seq().to_vec())?;
    let barcode_count = barcode_counts.entry(barcode.to_string()).or_insert(0);
    *barcode_count += 1;

    Ok(barcode + "_" + &barcode_count.to_string())
}

