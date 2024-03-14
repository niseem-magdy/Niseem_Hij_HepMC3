# Niseem_Hij_HepMC3

This Niseem Magdy code prepared on 12:15:2023

These are codes to run the Hijing model with hepmc3 and root Ttree options.

The "CMakeLists.txt" file will look up the C++ and root versions. You will need to set the HepMC3 dir manually in line 9.

In the same directory, use

cmake .

make

To run the code, use,

./Hij  Mode Decay FRAME N_EVENT IAP IZP IAT IZT EFRM BMIN BMAX

The code will take different arguments as,

Argument-01 is for 0-HepMC 1-ROOT-tree

Argument-02 is for 0-printing final particles only OR 1-for printing all particles; use 1 with HepMC.

Argument-03 is for 0-CMS 1 LAB

Argument-04 is for N-events 

Argument-05 is for IAP

Argument-06 is for IZP

Argument-07 is for IAT

Argument-08 is for IZT

Argument-09 is for EFRM-> Energy in GeV

Argument-10 is for the minimum impact parameter

Argument-11 is for the maximum impact parameter  

