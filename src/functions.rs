use std::{ f64::consts::PI };

use cryiorust::frame::Array;
use integrustio::integrator::Cake;
use rmpfit::MPFitter;


fn polarization(tth:f64, chi:f64, pfactor:f64)->f64{
    //0.5*(1.0 + np.cos(tthr)**2 - pfactor * np.cos(2.0 *chir) * (1.0 - np.cos(tthr)**2))
    let tthr = tth*PI/180.;
    let chir = chi*PI/180.;
    0.5*(1.+tthr.cos() - pfactor * (2.*chir).cos()*(1.- tthr.cos().powi(2)))
}

fn array_get_tthslice(a:&Array, tthindex:usize)-> Vec<f64>{
    let chilen = a.dim1();
    let tthlen = a.dim2();
    let mut slice: Vec<f64> = Vec::new();
    for chii in 0..chilen{
        let i = chii*tthlen + tthindex;
        let intensity = a[i];
        if intensity > 0.{
            slice.push(intensity);
        }
    }
    slice
}


fn getpolcake(tthrange:Vec<f64>, chirange:Vec<f64>, pfactor:f64)->Vec<f64>{
    let mut pvec: Vec<f64> = Vec::new();
    for chi in chirange.iter(){
        for tth in tthrange.iter(){
            pvec.push(polarization(*tth, *chi, pfactor));
        }
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
    let pmap = getpolcake(tthrange, chirange, pfactor);
    let cake = cake.cake.data();
    let mut fluosubcake: Vec<f64> = Vec::new();
    for (i, p) in cake.iter().zip(pmap.iter()){
        fluosubcake.push(i - fluo_k*p);
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


pub fn fluosub_curvefit(fluo_k:f64, cake:Cake, pfactor:f64, tthindex:usize)->Cake{
    let tthrange = cake.radial_positions;
    let chirange = cake.azimuthal_positions;
    let chilen = cake.cake.dim1();
    let tthlen=  cake.cake.dim2();
    let polcake = getpolcake(tthrange.to_vec(), chirange.to_vec(), pfactor);
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
    let mut init = [0., fluo_k];
    let _res = l.mpfit(&mut init).unwrap();
    let newfluok = init[1];
    let mut newcakevec = cake.cake.data().clone();
    for (c,p) in newcakevec.iter_mut().zip(polcakearray.data().iter()){
        *c = *c - *p * newfluok;
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
    newcake
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
            let f = params[0] + params[1] * *x; 
            *d = *y - f;
        }
        Ok(())
    }
    fn number_of_points(&self) -> usize {
        self.x.len()
    }
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