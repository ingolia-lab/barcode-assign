extern crate barcode_assign;
extern crate failure;

use std::io::Write;

use barcode_assign::pacbio_join::*;

fn main() {
    pacbio_join("./pacbio-190731", "./pacbio-190731").unwrap_or_else(|err| {
        std::io::stderr()
            .write(format!("{}\n", err).as_bytes())
            .unwrap();
        std::process::exit(1);
    });
}
