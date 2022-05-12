# lbxflow
Modelling two-dimensional complex flow using the lattice Boltzmann method

## lattice Boltzmann method
Lattice Boltzmann method is a computational technique that streams and
collides "psuedo particles", or particle frequency distributions, on a lattice
to approximate hydrodynamics.

## Dependencies

* Julia language (recently updated for julia 1.x.x)

### Julia language
* PyCall
* PyPlot
* HDF5
* JLD
* ArgParse
* ProfileView
* Roots
* YAML

## Quick Start

Install the both the [julia language](http://julialang.org/) and 
[python](https://www.python.org/). Note: the install script assumes an Ubuntu environment and that pip is installed.

    git clone https://github.com/grasingerm/lbxflow.git
    cd lbxflow
    julia scripts/install_dependencies.jl
    julia lbxflow.jl --help
    julia lbxflow.jl -vf example_sims/poiseuille/poiseuille_velocity-profile.yaml

## Documentation

An introductory [tutorial](https://github.com/grasingerm/lbxflow/tree/master/doc/tutorial/tutorial.pdf) 
is available in the `doc` directory. Currently documentation is sparse, but the
source files are all well documented and I am responsive to questions and feedback.
Example input files are available in the `example_sims` directory.

## Citations

There is a paper associated with lbxflow. If you found this code useful for your research, please be kind enough to cite it, using the DOI https://doi.org/10.1016/j.compfluid.2018.02.008, or the following BiBTeX entry:

    @article{grasinger2018numerical,
        title = {Numerical investigation of the accuracy, stability, and efficiency of lattice Boltzmann methods in simulating non-Newtonian flow},
        author = {Matthew Grasinger and Scott Overacker and John Brigham},
        journal = {Computers & Fluids},
        volume = {166},
        pages = {253-274},
        year = {2018},
        issn = {0045-7930},
        doi = {https://doi.org/10.1016/j.compfluid.2018.02.008}
    }

## Bugs

Please report all bugs and issues using the 
[github ticketing system](https://github.com/grasingerm/lbxflow/issues).

## TODO list
* [ ] add capability for 3D simulations
* [x] add animations and plots using gnuplot and Gadfly
* [x] implement free surface flow
* [x] break collision functions down into forcing method, constitutive eqn, etc.
* [x] create alt. stream/collide methods using grids of obstacle flags
* [x] write user documentation
* [ ] parallelize streaming and collision steps, multiscale mapping, and BCs
* [x] wrap simulation code in a module
* [ ] profile simulation steps to identify bottlenecks
* [x] write some tests?
* [x] write mrt collision functions with "hard coded" default inverse transformation matrix
* [x] verify Newtonian Poiseuille flow
* [x] debug BCs at inlet/outlet (bounceback, periodic)
* [x] implement MRT scheme for Bingham plastic flow
* [x] clean up api
* [x] make an autosave feature
* [x] should be able to "catch" from exceptions and/or keyboard interrupt
* [x] should be able to restart from interrupted model where left off
* [ ] clean up existing BCs / more general BCs
* [x] develop post processing and animations
* [ ] implement HB-Bingham constitutive model used in Fluent
