# srev-to-mt


## Description

This program convert a SREV nucleosone conformation into a MT theory  readable conformation.
The program makes a nucleosume chain of N beads or nucleosomes in the 
3SPN-AIGC representation, by combing the 3SPN-AIGC representation of a 
Nuclear Core Paricle with SREV nucleosome chain of N beads. The outputted 3SPN-AIGC chain and auxilary files 
is can be read in by the molecular theory program ( nucleosome-chain) allowing the computation of the electratostatic enviroment. 
 Program base on n
Uses the molecular theory approach to compute the structure and charge of a nucleosome n-mer in 3D in 
cartesian and oblique (prism or hexagonal) coordinates. 
The aqueous solution is characterized by pH, and salt concentration of NaCl, KCl, CaCl2 and MgCl2.

The outputted nucleosome chains or conformation has the resolution of the 3SPN-AIGC model, 
inwhich DNA is represented by 3 sites. Namely P(phosphate), S(sugar) and an nucliotude (A, T, G, and C ). The amino-acid residue is represented by one site, namely the alpha carbon.
An refined representation of the core histone octomer is can be obtained  (in the nucleosome-chain code). 
Here the amino-acids are represented 
on atomistic level. It includes explicit the locations of the CA, CB, CD, NZ, NH1 etc. Positions are obtained for the 1kx5 pdb file. 

This code is based on nucleosome code with as input has input a set of  3SPN-AIGC nucleosome chain.
As such many input variables are inherited form the nuclesoem-cahin program and not strictly nessary  and will be removed


### Prerequisites

The program uses the lapack eigenvalue routine.

### Installing

Modify Makefile such that in contains the appropriate linker flags to the sundials library.
Linker flags can be obtained from the examples directory of the sundials installation.



### Running
The program requires a number of input files similar to the input files of nucloesome-chain ( branch serial) program
A central configuration file is called 'input.in' whose spefication is detailed in nucloesome-chain  code. Here only the added new keywords ar listed.


| Keyword (type)     	      | Description                       |Value                                                   |
| --------------------------|-----------------------------------|--------------------------------------------------------|
| nnucl_sph `int`           | Number of nucleosome of outputted chain    |                                               |  
| nseg_sph `int`            | Number of segement/unit of outputted chain | optional                                      |
|nnucl `int`                | Number of nucleosomes of inputted chain    | value needs to be one                         |
|nseg  `int`		            | Number of segments of inputted chain       |                                               |
|cuantas `int`              | Number of conformations                    |   value needs to be one                       |

Other input files:
 
* Input file '\<trajSREV.xyz\>' contains the computed com coordianates of a SREV chain of N=nnucl_sph nucleosomes.

* Output files:

* Output file '\<trajchain_sph_0.xyz\>' contains the 3SPN-AIGC representation output conformation
* Output file '\<types_NCP.in\>' contains the sequence of the outputted conformation. It list the segment identifier number and char for every segment number.
* Output file '\<segcom_NCP.in\>' contains segment id number of residue that is assigned the COM of each (nucl_sph) nucleosome.
* Output file '\<nuc_orient_NCP.in\>' contains nnucl_sph triplets of segment id used for orientation of nucleosome AA elements.


he
## Built With



## Versioning

version 0.1 2-november-2026

## Authors

* **Rikkert J Nap**
