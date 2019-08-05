extern crate bio;
#[macro_use]
extern crate failure;
extern crate rust_htslib;

pub mod assign;
pub mod barcode_group;
pub mod bc_count;
pub mod bc_frag;
pub mod bc_seqs;
pub mod blasr;
pub mod depth;
pub mod fastq_pair;
pub mod flank_match;
pub mod frag_purity;
pub mod pacbio_join;
pub mod pacbio_reads;
pub mod purity;
