extern crate barcode_assign;

use barcode_assign::pacbio::*;

fn main() {
    let facs_conf = LibConf::new(MatchConf::new(b"GCTCGGAGATGTGTATAAGAGACAG",
                                                b"CTGTCTCTTATACACATCTGACGCC",
                                                4),
                                 MatchConf::new(b"CTATAGCACGACGCTCTTCCGATCT",
                                                b"GATCCTGTAGCCCTAGACTTGATAG",
                                                4),
                                 false);

    
}
