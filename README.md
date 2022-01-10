<!-- <img src="https://user-images.githubusercontent.com/1958085/148283531-18c1de46-7709-434d-9944-d0bbe64d73ee.png" width="100"> -->
<!-- <img src="https://user-images.githubusercontent.com/1958085/147726000-0be6dc67-b849-4cfb-b589-22e3229041b5.png" width="100"> -->

<img src="https://user-images.githubusercontent.com/1958085/148284111-430555a5-22c3-4012-a71e-f51fb0eb88c0.png" width="100">

![Tests](https://github.com/mvaliev/cdftpy/actions/workflows/tests.yaml/badge.svg)
# CDFTPY: Python package for performing classical density functional theory calculations for molecular liquids <!-- omit in toc --> 
#### Marat Valiev and Gennady Chuev<!-- omit in toc --> 

- [Introduction](#introduction)
- [Running from command line](#running-from-command-line)
  - [First steps](#first-steps)
  - [General usage:](#general-usage)
- [References:](#references)
  
## Introduction

Classical Density Functional Theory (CDFT) is a 
rigrous theoretical framework that builds statistical mechanics 
model of molecular liquid system in terms of its average density.
The `cdftpy` python package provides implementation of two CDFT methods -
the Renormalized Site Density Functional Theory (RSDFT) and Reference 
Interaction Site Model (RISM). At this stage of development, the code
allows to investigate the problem of ion solvation providing free energy
of solvation and solvent density profiles around the ion.


## Running from command line

---
### First steps
___

Running CDFT calculaton from the command line requires an input file.
Simple example of what it may look like for Cl- solvation calculation 
is given below: 

```
<solute>
#name   sigma(A)        eps(kj/mol)    charge(e)  
  Cl       4.83         0.05349244     -1.0      
```

Once the input file has been prepared, CDFT calculation can be run
using `cdft1d` executable script as shown below 

```
cdft1d <input_file>
```
In this simple case, RSDFT calculation will be performed using `s2` water model. The output will contain solvation free energy as well as peak
analysis of solvent density around the ion, e.g.

```
....
-----------------------------
   Self-consistent cycle     
-----------------------------
iter  d_g         Free Energy 
0     4.21e+00   -8.2510119
10    9.02e-01   -195.0800613
20    1.17e-02   -297.0235445
30    2.71e-03   -297.1768922
40    2.87e-06   -296.7840683
41    8.80e-08   -296.7836460

Reached convergence, d_g < 1e-07


Total Free Energy               -296.783646
----------------------------------
Solvent density structure analysis
----------------------------------
O 1st peak position/height: 3.13 4.68  
O 2nd peak position/height: 5.93 1.38  
O 1st min position/height: 4.44 0.54  
  
H 1st peak position/height: 2.16 8.00  
H 2nd peak position/height: 4.93 1.29  
H 1st min position/height: 3.26 0.37 
```

Single ion solvation calculation with RSDFT flavor of CDFT
can be performed as

where example input file for Cl- ion is given by

```
<solute>
#name   sigma(A)        eps(kj/mol)    charge(e)  
  Cl       4.83         0.05349244     -1.0      
<simulation>
  solvent s2
  method rsdft
  tol 1.0E-7
  max_iter 200
  output_rate 10
  rcoul 1.25
  rmax 100
```
The output will contain solvation free energy as well as peak
analysis of solvent density around the ion, e.g.

```
...
-----------------------------
   Self-consistent cycle     
-----------------------------
iter  d_g         Free Energy 
0     1.13e+00   -258.8752239
10    4.14e+01   -4839.5648029
20    9.29e-02   -292.7447959
30    1.53e-03   -301.9588939
40    7.17e-05   -301.5734464
50    2.65e-07   -301.5653213
53    9.82e-08   -301.5653695

Reached convergence, d_g < 1e-07


Total Free Energy               -301.565369 kj/mol
----------------------------------
Solvent density structure analysis
----------------------------------
O 1st peak position/height: 3.12 4.785  
O 2nd peak position/height: 5.91 1.384  
O 1st min position/height: 4.43 0.543  
  
H 1st peak position/height: 2.15 8.291  
H 2nd peak position/height: 4.92 1.301  
H 1st min position/height: 3.24 0.359  

```
The same calculation but with RISM methodology
can be run as

    cdft1d -m rism <input_file>

___
### General usage:
___
```
Usage: cdft1d [OPTIONS] INPUT_FILE

Options:
  -m, --method [rism|rsdft]       
    Calculation method (default:rsdft)
  -s, --solvent <solvent model>   
    solvent model  [default: s2]
  -a, --adjust [charge|sigma|eps] <value>
    adjust solute parameters
  -r, --range [charge|sigma|eps] <values>
    Run calculation over the range of solute
    "charge","sigma","eps" values. Values could
    specified as comma separated sequence (e.g.
    0,0.5,..) or in triplets notation
    [start]:stop:nsteps. To avoid  issues with
    blank spaces, it is recommended that values
    are enclosed in double quotes.
  -o, --output
    generate data output (currently only site density)
  -d, --dashboard [filename]      
    Generate dashboard for analysis. 
    The dashboard will be saved as html file 
    under the name provided by optional argument. 
    In the absence of the latter dashboard will be
    open in browser
  --version                       
    display version and exit
  --help                         
    Show this message and exit
```
___
## References:
___
G.N. Chuev, M. V. Fedotova and M. Valiev,
 Renormalized site density functional theory,
 J. Stat. Mech. (2021) 033205, https://doi.org/10.1088/1742-5468/abdeb3

G.N. Chuev, M. V. Fedotova and M. Valiev,
Chemical bond effects in classical site density 
functional theory of inhomogeneous molecular liquids. 
J. Chem Phys. 2020 Jan 31;152(4):041101,
https://doi.org/10.1063/1.5139619

M. Valiev and G.N. Chuev,
 Site density models of inhomogeneous classical molecular liquids,
 J. Stat. Mech. (2018) 093201,
https://doi.org/10.1088/1742-5468/aad6bf

 
 
