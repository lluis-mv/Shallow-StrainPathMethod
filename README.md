 ## StrainPathMethod
Matlab/Octave code to compute displacements, velocities and strains due to the penetration of a tube using the Strain Path Method (Baligh, 1985).

The Shallow Strain Path Method (Sagaseta et al, 1991) is also included.

Both codes follow the approach proposed by Sagaseta et al (1991): the tube moves from the far field to the center of the domain.

Elliptical integrals are computed according to Hollos (2012) [http://www.exstrom.com/math/elliptic/ellipint.html] since it is by far more eficient than using the Matlab build-in elliptical integrals.

Data postprocessing is done in GiD.



### Octave compatibility with symbolic generation

(The code is fully compatible with Octave. However, the symbolic generation of files is slow and probably should be changed)

The symbolic toolbox of octave might be installed (in Linux) by typping in the octave command window:

> pkg install -forge symbolic

and loaded with:
   
> pkg load symbolic

It might be necessary to install sympy. In a terminal, enter the following line:

> sudo apt-get install python-sympy
