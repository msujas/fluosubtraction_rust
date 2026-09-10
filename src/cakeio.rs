use std::{ fs::File, io::Write, path::Path, sync::Arc };

use cryiorust::{edf::Edf, frame::{Array, Frame, HeaderEntry}};
use integrustio::integrator::Cake;

fn linspace(start:f64,stop:f64, npoints:usize)->Vec<f64>{
    let spacing = (stop-start)/((npoints-1) as f64);
    let mut outvec : Vec<f64> = Vec::new();
    for i in 0..npoints{
        outvec.push(start + (i as f64)*spacing);
    }
    outvec
}

struct Pattern1d{
    tth: Vec<f64>,
    intensity: Vec<f64>,
    sigma: Vec<f64>
}

fn parse_bubblepattern(bubblepattern_s:&String)-> Pattern1d{
    let mut ttharray : Vec<f64> = Vec::new();
    let mut intarray :Vec<f64> = Vec::new();
    let mut sigarray : Vec<f64> = Vec::new();
    for item in bubblepattern_s.split("\t"){

        for (count,val) in item.split(" ").enumerate(){
            let value = val.parse::<f64>().unwrap();
            match count {
                0 => ttharray.push(value),
                1 => intarray.push(value),
                2 => sigarray.push(value),
                _ => continue,
            }
        }
    }
    Pattern1d { tth: ttharray, intensity: intarray, sigma: sigarray }
}

struct IntegrationRange{
    tth0: f64,
    tthend: f64,
    chi0: f64,
    chiend: f64,
}

fn parse_bubblecake(bubblecake_s:&String)-> IntegrationRange{
    let mut bcsplit = bubblecake_s.split(" ");
    let tth0 = bcsplit.next().unwrap().parse::<f64>().unwrap();
    let tthend = bcsplit.next().unwrap().parse::<f64>().unwrap();
    let chi0 = bcsplit.next().unwrap().parse::<f64>().unwrap();
    let chiend = bcsplit.next().unwrap().parse::<f64>().unwrap();
    IntegrationRange { tth0, tthend, chi0, chiend }
}

pub struct CakeReadError;

pub fn readcake(cakefile:&String)-> Result<Cake, CakeReadError>{
    if !cakefile.ends_with(".edf"){
        let cakeext = Path::new(cakefile).extension().unwrap();
        eprintln!("cake file must be in .edf format. Given file is .{cakeext:?}");
        return Err(CakeReadError)
    };
    let im = match Edf::open(cakefile) {
        Ok(e) => e,
        Err(_e) => {eprintln!("couldn't read file"); return Err(CakeReadError)}
    };
    let a = im.array().data().clone();
    let bubblepattern_s = match im.header().get("Bubble_pattern"){
        Some(HeaderEntry::String(s)) => s,
        _ => {println!("couldn't read Bubble_pattern"); return Err(CakeReadError)},
    };
    let pattern = parse_bubblepattern(bubblepattern_s);
    let bubblecake = match im.header().get("Bubble_cake"){
        Some(HeaderEntry::String(s)) => s,
        _ => {println!("couldn't read Bubble_cake"); return Err(CakeReadError)},
    };
    let chisize = im.dim1();
    let tthsize = im.dim2();
    let irange = parse_bubblecake(bubblecake);
    let _tth0 = irange.tth0;
    let _tthend = irange.tthend;
    let chi0 = irange.chi0;
    let chiend = irange.chiend;
    let chirange = linspace(chi0, chiend, chisize);
    let mut cake:Cake = Default::default();
    cake.azimuthal_positions = Arc::new( chirange);
    cake.radial_positions = Arc::new(pattern.tth);
    cake.cake = Array::with_data(chisize, tthsize, a);
    cake.radial.intensity = pattern.intensity;
    cake.radial.sigma = pattern.sigma;
    Ok(cake)

}

pub fn save1d(fname:&Path, tthrange: &Vec<f64>, vec1d: &Vec<f64>, sigma : Option<&Vec<f64>>){
    let mut outstring = String::new();
    let mut x:f64;
    let mut y:f64;
    let mut e:f64;
    let dosig:bool = match sigma  {
        None => false,
        Some(_s) => true
    };
    for i in 0..tthrange.len(){
        x = tthrange[i];
        y=vec1d[i];
        outstring = outstring + &String::from(format!("{x:.6} {y:.6}"));
        if dosig{
            e = sigma.unwrap()[i];
            outstring = outstring + &String::from(format!(" {e:.6}"));
            }
        outstring = outstring + &String::from("\n");
        }
    println!("saving 1d pattern to {:?}", &fname);
    
    let mut file = File::create(&fname).expect(&format!("error creating file {:?}",&fname));
    file.write(outstring.as_bytes()).unwrap();    
}