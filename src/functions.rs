use std::{ f64::consts::PI, fs::File, io::Write, path::Path, sync::Arc };

use cryiorust::{edf::Edf, frame::{Array, Frame, HeaderEntry}};
use integrustio::integrator::Cake;
use rmpfit::MPFitter;


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

fn polarization(tth:f64, chi:f64, pfactor:f64)->f64{
    //0.5*(1.0 + np.cos(tthr)**2 - pfactor * np.cos(2.0 *chir) * (1.0 - np.cos(tthr)**2))
    let tthr = tth*PI/180.;
    let chir = chi*PI/180.;
    0.5*(1.+tthr.cos().powi(2) - pfactor * (2.*chir).cos()*(1.- tthr.cos().powi(2)))
}

fn array_get_tthslice(a:&Array, tthindex:usize)-> Vec<f64>{
    let chilen = a.dim1();
    let tthlen = a.dim2();
    let mut slice: Vec<f64> = Vec::new();
    for chii in 0..chilen{
        let i = chii*tthlen + tthindex;
        let intensity = a[i];
        slice.push(intensity);
    }
    slice
}


pub fn getpolcake(tthrange:Vec<f64>, chirange:Vec<f64>, pfactor:f64, fname:Option<String>)->Vec<f64>{
    let mut pvec: Vec<f64> = Vec::new();
    for chi in chirange.iter(){
        for tth in tthrange.iter(){
            pvec.push(polarization(*tth, *chi, pfactor));
        }
    }
    if let Some(fname) = fname{
        let a = Array::with_data(chirange.len(),tthrange.len(), pvec.clone());
        let mut cake:Cake = Default::default();
        let (intensity, sigma) = intensity_from_array(&a);
        cake.cake = a;
        cake.azimuthal_positions = Arc::new(chirange);
        cake.radial_positions = Arc::new(tthrange.clone());
        cake.radial.intensity = intensity;
        cake.radial.sigma = sigma;
        cake.radial.positions = Arc::new(tthrange.clone());
        cake.store(fname,None).unwrap();
    }
    pvec
}


fn intensity_from_array(array:&Array)-> (Vec<f64>, Vec<f64>){
    let a = array.data();
    let chilen = array.dim1();
    let tthlen = array.dim2();
    let mut intensity:Vec<f64> = vec![0.; tthlen];
    let mut divvec: Vec<f64> = vec![0.; tthlen];
    let mut sigma: Vec<f64> = Vec::new();
    for chii in 0..chilen{
        for tthi in 0..tthlen{
            let aindex = chii * tthlen + tthi;
            if a[aindex] > 0.{
                intensity[tthi] += a[aindex];
                divvec[tthi] += 1.;
            }
        }
    }

    for (i,d) in intensity.iter_mut().zip(divvec.iter_mut()){
        if *d <= 0.{
            *d = 1.;
        }
        *i = *i/ *d;
        sigma.push(f64::powf(*i, 0.5)/ *d);
    }
    (intensity, sigma)
}

pub fn fluosub_cake(cake:Cake, pfactor:f64, fluo_k: f64)->Cake{
    let tthrangea = cake.radial_positions;
    let chirangea = cake.azimuthal_positions;
    let tthrange = tthrangea.to_vec();
    let chirange = chirangea.to_vec();
    let chisize = cake.cake.dim1();
    let tthsize = cake.cake.dim2();
    let pmap = getpolcake(tthrange, chirange, pfactor,None);
    let cake = cake.cake.data();
    let mut fluosubcake: Vec<f64> = Vec::new();
    for (i, p) in cake.iter().zip(pmap.iter()){
        let isub = i-fluo_k/p;
        if isub > 0.{
            fluosubcake.push(i - fluo_k/p);
        }
        else {
            fluosubcake.push(0.);
        }
    }
    let a = Array::with_data(chisize, tthsize, fluosubcake);
    let (intensity, sigma) = intensity_from_array(&a);
    let mut newcake:Cake = Default::default();
    newcake.cake = a;
    newcake.azimuthal_positions = chirangea.clone();
    newcake.radial_positions = tthrangea.clone();
    newcake.radial.intensity = intensity;
    newcake.radial.sigma = sigma;
    newcake.radial.positions = tthrangea.clone();
    newcake
}


pub fn fluosub_curvefit(fluo_k0:f64, cake:Cake, pfactor:f64, tthindex:usize)->(Cake, f64){
    let tthrange = cake.radial_positions;
    let chirange = cake.azimuthal_positions;
    let chilen = cake.cake.dim1();
    let tthlen=  cake.cake.dim2();
    let polcake = getpolcake(tthrange.to_vec(), chirange.to_vec(), pfactor,None);
    let polcakearray = Array::with_data( chilen, tthlen, polcake);
    let pslice = array_get_tthslice(&polcakearray, tthindex);
    let cakeslice = array_get_tthslice(&cake.cake, tthindex);
    assert!(&cakeslice.len() == &pslice.len());
    let mut cakeslicecut: Vec<f64> = Vec::new();
    let mut pslicecut :Vec<f64> = Vec::new();
    for (c,p) in cakeslice.iter().zip(pslice.iter()){
        if *c > 0.{
            cakeslicecut.push(*c);
            pslicecut.push(*p);
        }
    }
    let mut l = Linear{x:pslicecut,y: cakeslicecut};
    let mut init = [0., fluo_k0];
    let _res = l.mpfit(&mut init).unwrap();
    let newfluok = init[1];
    println!("new fluok: {newfluok}");
    let mut newcakevec = cake.cake.data().clone();
    for (c,p) in newcakevec.iter_mut().zip(polcakearray.data().iter()){
        if *c > 0.{
            *c = *c -  newfluok/ *p;
        }
    }
    let mut newcake : Cake = Default::default();
    let newarray = Array::with_data(chilen, tthlen, newcakevec);
    let (i, sig) = intensity_from_array(&newarray);
    newcake.cake = newarray;
    newcake.azimuthal_positions = chirange.clone();
    newcake.radial_positions = tthrange.clone();
    newcake.radial.intensity = i;
    newcake.radial.sigma = sig;
    newcake.radial.positions = tthrange.clone();
    (newcake, newfluok)
}

struct Linear{
    x: Vec<f64>, // x - polarisation slice, y - data slice
    y: Vec<f64>,
}

impl MPFitter for Linear{
    fn eval(&mut self, params: &[f64], deviates: &mut [f64]) -> rmpfit::MPResult<()> {
        for ((d, x), y) in deviates
            .iter_mut()
            .zip(self.x.iter())
            .zip(self.y.iter())
        {
            let f = params[0] + params[1] / *x; 
            *d = *y - f;
        }
        Ok(())
    }
    fn number_of_points(&self) -> usize {
        self.x.len()
    }
}


pub fn cakeget1d(cakearray: &Array)-> Vec<f64>{
    let chisize = cakearray.dim1();
    let tthsize = cakearray.dim2();
    let mut pattern1d: Vec<f64> = vec![0.;tthsize];
    for t in 0..tthsize{
        let mut tthslice = 0.;
        let mut div = 0.;
        for c in 0..chisize{
            let index = t + c*tthsize;
            let value = cakearray.data()[index];
            if value > 0.{
                tthslice += value;
                div += 1.
            }
        }
        if div > 0.{
            pattern1d[t] = tthslice/div;
        }
    }
    pattern1d
}

pub fn save1d(fname:String, tthrange: &Vec<f64>, vec1d: &Vec<f64>, sigma : Option<&Vec<f64>>){
    let mut outstring = String::new();
    //for (x,y ) in  tthrange.iter().zip(vec1d.iter()){
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
    println!("saving 1d pattern to {}", &fname);
    
    let mut file = File::create(&fname).expect(&format!("error creating file {:?}",&fname));
    file.write(outstring.as_bytes()).unwrap();    
}


/*
/// possible alternative minimisation
fn fluosub_lsquare(fluo_k: f64, cake:Cake, pfactor:f64)->f64{
    let newcake = fluosub_cake(cake, pfactor, fluo_k);
    let chilen = newcake.cake.dim1();
    let tthlen = newcake.cake.dim2();
    let index = tthlen*96/100;
    let mut slice : Vec<f64> = Vec::new();
    let mut mean: f64 = 0.;
    let mut div = 0.;
    let a = newcake.cake.data();
    for chii in 0..chilen{
        let i = chii*tthlen + index;
        let intensity = a[i];
        if intensity > 0.{
            slice.push(intensity);
            mean += intensity;
            div += 1.;
        }
    }
    mean = mean/div;

    let mut sum = 0.;
    for item in slice{
        sum += f64::powi(item-mean,2);
    }
    sum
}  */ 

#[cfg(test)]
mod tests{
    use super::*;
    #[test]
    fn polcaketest(){
        let tthrange:Vec<f64> = (1..60).map(f64::from).collect();
        let chirange:Vec<f64> = (0..359).map(f64::from).collect();
        getpolcake(tthrange, chirange, 0.85, Some(String::from("./polcake.edf")));
    }
}