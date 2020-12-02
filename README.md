# subroutine_saamd:

Subroutines used in self-adaptive accelerated molecular dynamics (SAAMD) 
code, which is based on mdcask code. 

The subroutines are used to (1) to set the boost atoms in AV, (2) calculate
the displacement of atoms, (3) calculate the boost energy and force for 
each atom in AV and (4) check the state to determine the maximum displacement
has been satisfied. 

Based on these results, the total energy and force can be calculated 
according to equations in our paper (https://iopscience.iop.org/article/10.1088/1361-648X/aa574b). 
The system evolves then according to the energy from the empirical potential
and biased potential.

These subroutines can be compiled with MOLDY code.
