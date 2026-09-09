use clap::Parser;
use fluosubtraction_rust::functions::{cakeget1d, fluosub_curvefit, readcake, save1d};

#[derive(Parser, Debug)]
#[command(version,about="program for correcting fluorescence from cake files", long_about=None)]
struct Params{
    /// file to correct fluorescence
    #[arg(short, long)]
    pub filename: String,
    
    /// polarisation factor
    #[arg(short, long, default_value_t=0.85)]
    pub pfactor: f64,

    /// 2theta index to correct on
    #[arg(short, long, default_value_t=4500)]
    pub tthindex: usize,

    /// initial fluo constant to use
    #[arg(short, long, default_value_t=1.)]
    pub k0: f64, 
}

fn main(){
    let ap = Params::parse();
    let filename = ap.filename;
    let pfactor = ap.pfactor;
    let tthindex = ap.tthindex;
    let k0 = ap.k0;

    let cake = readcake(&filename);
    let (newcake, _fluok) = fluosub_curvefit(k0, cake, pfactor, tthindex);
    let newfilename = &filename.replace(".edf", "_fluosub.edf");
    let vec1d = cakeget1d(&newcake.cake);
    let tth = newcake.radial_positions.to_vec();
    let newfile1d = filename.replace(".edf", "_fluosub.xy");
    save1d(newfile1d, &tth, &vec1d, None);
    println!("saving fluo sub cake to {newfilename}");
    newcake.store(newfilename, None).unwrap();
}