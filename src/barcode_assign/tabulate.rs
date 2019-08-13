use failure;

pub struct CLI {
    pub inputs: Vec<String>,
    pub output: String,
    pub mintotal: Option<usize>,
    pub minsamples: Option<usize>,
    pub mininsample: Option<usize>,
    pub omitfile: Option<String>
}

impl CLI {
    pub fn run(&self) -> Result<(), failure::Error> {
        unimplemented!()
    }
}
