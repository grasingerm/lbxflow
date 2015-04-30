# lbxflow
Modelling complex flow using lattice Boltzmann method

## lattice Boltzmann method
Lattice Boltzmann method is a computational technique that streams and
collides "psuedo particles", or particle frequency distributions, on a lattice
to approximate hydrodynamics.

## TODO list
* [x] implement MRT scheme for Bingham plastic flow
* [ ] implement free surface flow
* [ ] clean up api
* [ ] clean up existing BCs / more general BCs
* [ ] develop post processing and animations
* [ ] consider "DSL" for mapping particle to macroscale and vice versa
* [ ] implement Jonas Latt BCs
* [ ] parallelize streaming and collision steps, multiscale mapping, and BCs
* [ ] profile simulation steps to identify bottlenecks
* [ ] write some tests?
