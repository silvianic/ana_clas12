# ana_clas12

Analysis code on CLAS12 data for pDVCS and nDVCS channels. 
At its heart it is a TSelector reading the hipo4 file converted to ROOT by Nick Tyler's hipo2root. 
It has a main code (ana.cxx and ana_mc.cxx) which reads an input file (runListNUM) containing the files to be analyized, chains them and gives them to the TSelector. 
It produces an output ROOT files which contains 3 trees: ep, nDVCS, pDVCS. 
Right now there are TWO CODES, one for data, one for MC. 

Usage, on IPN machines:

To compile:

use root
make  (for the MC version: make -f Makefile_mc)

To run:
Create a file called runListNUM (where NUM is a number of your choice) containing the desired input files. The command to run is:

./ana NUM NUM2 (where NUM2 is how many of the files in the list we want the code run on) for data

./ana_mc NUM NUM2 (where NUM2 is how many of the files in the list we want the code run on) for MC




