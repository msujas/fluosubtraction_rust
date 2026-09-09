A library and command line application for doing fluorescence correction on area detectors. The library contains functions for correcting the fluorescence using a least squares optimisation on cake arrays (2D azimuthally integrated images), and a function for reading cake files.

The command line application (cakefluosub) takes the filename as a positional argument, polarisation factor, 2theta index to optimise, and initial fluorescence constant to correct on a single file.

```
cakefluosub --help

program for correcting fluorescence from cake files

Usage: cakefluosub.exe [OPTIONS] <FILENAME>

Arguments:
  <FILENAME>  cake file to correct fluorescence on

Options:
  -p, --pfactor <PFACTOR>    polarisation factor [default: 0.85]
  -i, --tthindex <TTHINDEX>  2theta index to correct on (default 90% of number of bins)
  -k, --k0 <K0>              initial fluo constant to use [default: 1]
  -h, --help                 Print help
  -V, --version              Print version
```