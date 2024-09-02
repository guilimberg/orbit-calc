# orbit-calc
My implementation of orbit calculation with AGAMA.

This code uses [Astropy](https://www.astropy.org/), [pyGaia](https://github.com/agabrown/PyGaia), and [AGAMA](http://agama.software/) ([Vasiliev 2019](https://academic.oup.com/mnras/article/482/2/1525/5114593)) to integrate orbits in an axisymmetric Galactic model potential and, then, extract several kinematic/dynamical quantities, including angular momenta, actions, and orbital energies.

The standard model potential I adopt is [McMillan 2017](https://academic.oup.com/mnras/article/465/1/76/2417479), which is described in the file "McMillan17.ini" in this repository. 

If you have any questions, get in contact with me: guilherme.limberg [at] usp.br.

**UPDATE Sep 2024**
Unfortunately, this requires running older versions of some packages. I will eventually update the code to run with contemporary Python package versions. For now, I was able to run it with Astropy v5.3.4, NumPy v1.23.4, and Scipy v1.11.1. 
