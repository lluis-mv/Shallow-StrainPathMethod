 ## Shallow-StrainPathMethod
 [![License][license-image]][license]
 
 [license-image]: https://img.shields.io/badge/license-BSD-green.svg?style=flat
[license]: https://github.com/lluis-mv/Shallow-StrainPathMethod/blob/master/LICENSE
 
Matlab/Octave code to compute displacements, velocities and strains due to the penetration of a tube using the Strain Path Method ([Baligh, 1985](http://doi.org/10.1061/(ASCE)0733-9410(1985)111:9(1108))).

The Shallow Strain Path Method ([Sagaseta et al, 1991](http://doi.org/10.1002/(SICI)1096-9853(199710)21:10<687::AID-NAG897>3.0.CO;2-3)) is also considered.

Both codes follow the approach proposed by [Sagaseta et al (1998)](http://doi.org/10.1002/(SICI)1096-9853(199710)21:10<687::AID-NAG897>3.0.CO;2-3): the tube moves from the far field to the center of the domain.

Elliptical integrals are computed according to [Hollos (2012)](http://www.exstrom.com/math/elliptic/ellipint.html) since it is by far more eficient than using the Matlab build-in elliptical integrals.

Data postprocessing is done in [GiD](http://www.gidhome.com/).



### Octave compatibility with symbolic generation

The code is fully compatible with Octave. However, the symbolic generation of these files is slow and probably should be changed. 

The symbolic toolbox of octave might be installed (in Linux) by typping in the octave command window:

> pkg install -forge symbolic

and loaded with:
   
> pkg load symbolic

It might be necessary to install sympy. In a terminal, enter the following line:

> sudo apt-get install python-sympy

## License and citation

This code is distributed under BSD 3-Clause License. 


If you find these codes useful for your work, please cite as follows:
- Monforte, L., Arroyo, M., Carbonell, J.M. & Gens, A (2020) Large strain analysis of undrained smooth tube samplig. Géotechnique. (Accepted)
