#[macro_use]
extern crate error_chain;
extern crate rust_htslib;

pub mod barcode_group;

pub mod errors {
    error_chain!{
        errors {
            NoBarcode(qname: Vec<u8>) {
                description("no barcode in read name")
                display("no barcode in read name {}",
                        ::std::str::from_utf8(&qname).unwrap_or("<Not utf8>"))
            }
        }
        foreign_links {
            BamRead(::rust_htslib::bam::ReadError);
        }
    }
}
