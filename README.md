# lbxflow
Modelling complex flow using lattice Boltzmann method

## lattice Boltzmann method
Lattice Boltzmann method is a computational technique that streams and
collides "psuedo particles", or particle frequency distributions, on a lattice
to approximate hydrodynamics.

## TODO list
* pass the inverse of the transformation matrix via function call or cache
  result in order to make MRT simulation more computationally efficient
* parallelize streaming and collision steps, multiscale mapping, and BCs
* profile simulation steps to identify bottlenecks
* implement MRT scheme for Bingham plastic flow
* implement free surface flow
