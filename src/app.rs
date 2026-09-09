use std::{fs::create_dir, path::{Path, PathBuf}, process::exit};

use clap::{Parser};
use fluosubtraction_rust::functions::{cakeget1d, fluosub_curvefit, readcake, save1d};

#[derive(Parser, Debug)]
#[command(version,about="program for correcting fluorescence from cake files", long_about=None)]
struct Params{
    /// cake file to correct fluorescence on
    pub filename: String,
    
    /// polarisation factor
    #[arg(short, long, default_value_t=0.85)]
    pub pfactor: f64,

    /// 2theta index to correct on (default 90% of number of bins)
    #[arg(short='i', long)]
    pub tthindex: Option<usize>,

    /// initial fluo constant to use
    #[arg(short, long, default_value_t=1.)]
    pub k0: f64,

    /// optional subdir for saved files
    #[arg(short, long)]
    pub subdir: Option<String>,
}

fn main(){
    let ap = Params::parse();
    let filename = ap.filename;
    let pfactor = ap.pfactor;
    let tthindex = ap.tthindex;
    let k0 = ap.k0;
    let subdir = ap.subdir;

    let cake = match readcake(&filename){
        Ok(c) => c,
        Err(_e) => {println!("couldn't read cake, exiting"); exit(1)}
    };
    let tth = cake.radial_positions.to_vec();
    let tthi = match tthindex{
        Some(i) => i,
        None => tth.len()*90/100,
    };
    println!("optimising fluorescence correction using tthindex: {tthi}, k0: {k0}, polarisation factor: {pfactor},\non file: {filename}");
    let (newcake, _fluok) = fluosub_curvefit(k0, cake, pfactor, tthi);
    let fpath = Path::new(&filename);
    let mut dirname = PathBuf::from(fpath.parent().unwrap());
    let filestem = fpath.file_stem().unwrap();
    if let Some(s) = subdir{
        dirname.push(s);
        if !dirname.exists(){
            let _ = create_dir(&dirname);
        }
    };
    let outfile = format!("{}_fluosub.edf", &filestem.to_str().unwrap());
    let outfile1d = format!("{}_fluosub.xy",&filestem.to_str().unwrap());
    let newfilename = dirname.join(outfile);
    let newfile1d = dirname.join(&outfile1d);

    let vec1d = cakeget1d(&newcake.cake);
    
    save1d(&newfile1d, &tth, &vec1d, None);
    println!("saving fluo sub cake to {newfilename:?}");
    newcake.store(newfilename, None).unwrap();
}