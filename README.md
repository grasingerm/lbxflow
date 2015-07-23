# lbxflow
Modelling complex flow using lattice Boltzmann method

## lattice Boltzmann method
Lattice Boltzmann method is a computational technique that streams and
collides "psuedo particles", or particle frequency distributions, on a lattice
to approximate hydrodynamics.

## Dependencies
* PyCall
* PyPlot
* HDF5
* ArrayViews
* ProfileView
* PyYaml       (python library)
* numpy        (python library)
* matplotlib   (python library)

## Quick Start
    git clone https://github.com/grasingerm/lbxflow.git
    cd lbxflow
    julia scripts/install_dependencies.jl
    julia lbxflow.jl --help
    julia lbxflow.jl -vf example_sims/ex/ex1.yaml

## Input Files

### Overview

Input files are created according to the [YAML](http://yaml.org/) format. If you
need to check if your YAML is wellformed, check out this 
[online parser](http://yaml-online-parser.appspot.com/) with helpful error
reporting.

### preamble
    preamble: >
    const ni = 128;
    const nj = 64;
    const pgrad = -5.6e-6;
    const F = [-pgrad; 0.0];
    const mu_p = 0.2;
    const tau_y = 4.0e-5;
    const m = 1.0e8;
    const max_iters = 11;
    const tol = 1e-3;
    const nsteps = 20000;
    const id = "poise-chen-e-000004";
    const datadir = joinpath("data","poise","chen","explicit","000004");
    const constit_rel_f = init_constit_mrt_bingham_explicit(mu_p, tau_y, m, 1.0e-9, 1.0);
    const forcing_kf = init_korner_Fk(F);

The `preamble` section of an input file is parsed and evaluated in the global
scope. This means the preamble can contain any julia code, and that anything
defined in the preamble will be available throughout the input file. Because
the preamble is evaluated globally, it should avoid *lbxflow* names, such
as `simulate!`, `col_srt!`, `north_bounce_back!`, etc. Also, the pound sign
(`#`) should not be used in the preamble. Because this starts a comment or
comment block in julia, including it in the preamble will lead to unpredictable
behavior.

### datadir

    datadir:    { value: datadir,   expr: true    }

The data directory is an output directory that will be created if it does not
exist. Note: `expr` flag determines whether or not an input is parsed and
evaluated, or if the value determined by the YAML specification is used.

### material properties
    # material init
    rho_0:      { value: 1.0,       expr: false   }
    nu:         { value: mu_p,      expr: true    }

`rho_0` specifies the ref/initial density. `nu` specifies the reference
kinematic viscosity.

### lattice configuration
    # lattice configuration
    dx:         { value: 1.0,       expr: false   }
    dt:         { value: 1.0,       expr: false   }
    ni:         { value: ni,        expr: true    }
    nj:         { value: nj,        expr: true    }

`dx` specifies the distance between lattice nodes. `dt` specifies the time
step. `ni` specifies the number of nodes in the x-direction. `nj` specifies
the number of nodes in the y-direction.

### simulation parameters
    # simulation parameters
    simtype:    default
    nsteps:     { value: nsteps,    expr: true    }
    col_f:      init_col_mrt(constit_rel_f, forcing_kf, S_luo)

`nsteps` specifies the number of steps to simulate. `col_f` is a function
that performes LBM collisions. Collision functions can be initialized using
the building blocks contained in `col/constitutive.jl`, `col/forcing.jl`,
and `col/mrt_matrices.jl`, and the constructors contained in `col/modcol.jl`.
Or, paramters can be bound to collision functions defined in `col/stdcol.jl`.
Note: collision functions must have the interface:
`(sim::AbstractSim, bounds::Matrix{Int64})`

### boundaries
    # boundaries
    sbounds:
      value: "[1 ni 1 nj;]'"
      expr: true

    cbounds:
      value: "[1 ni 1 nj;]'"
      expr: true

`sbounds` specifies boundaries on which streaming takes place. `cbounds`
specifices boundaries on which collisions take place. The data is structured:
`[i_min_1 i_min_2 ... i_min_n; i_max_1 i_max_2 ... i_max_n; j_min_1 j_min_2 ...
j_min_n; j_max_1 j_max_2 ... j_max_n;]`

### boundary conditions
    # boundary conditions
    bcs:
      - north_bounce_back!
      - south_bounce_back!
      - periodic_east_to_west!

A list of boundary conditions to be applied to the lattice. Boundary condition
functions can be found in `inc/boundary.jl`

### callback functions
    # callback functions
    callbacks:
      - print_step_callback(100, id)
      - write_jld_file_callback(datadir, convert(Int64, round(nsteps/20)))

A list of callback functions. The interface of these function is
`(sim::AbstractSim, k::Int)` where `k` is the frequency (in time steps) in
which the callback is performed. Useful callback functions can be found in
the `inc/io` directory.


### finally
    # clean-up, backup, write out
    finally:
      - >
        (sim::Sim, k::Int) -> begin
          writedlm(joinpath(datadir, "ux_profile.dsv"), 
            extract_ux_prof_callback(convert(Int64, round(ni/2)))(sim), 
            ",");
        end
      - write_jld_file_callback(datadir)

The interface of these functions is consitent with callback functions, with
the difference that these are only executed at the end of the simulation. If
the simulation is ended on an exception or `Ctrl+c` interrupt, then `finally`
functions should still execute.

### test for term (optional)
    # test for conditions to end simulation
    test_for_term: is_steadystate_x

`test_for_term` specifies functions that, if return true, terminate the
simulation. These are usually functions that check for steady-state conditions.
Examples can be found in `inc/convergence.jl`

### Full Example

    # author:         
    #  name:         Matthew Grasinger
    #  email:        grasingerm at gmail dot com
    #
    # date-created:   2015-06-12
    #
    # description: >
    #   Poiseuille bingham plastic flow for verification. Based on parameters from
    #   Chen et. al. 2014, Simulations of Bingham plastic flows with the multiple-
    #   relaxation-time latice Boltzmann model. Yield stress is 4e-5.

    version: 0.2.1

    preamble: >
      const ni = 128;
      const nj = 64;
      const pgrad = -5.6e-6;
      const F = [-pgrad; 0.0];
      const mu_p = 0.2;
      const tau_y = 4.0e-5;
      const m = 1.0e8;
      const max_iters = 11;
      const tol = 1e-3;
      const nsteps = 20000;
      const id = "poise-chen-e-000004";
      const datadir = joinpath("data","poise","chen","explicit","000004");
      const constit_rel_f = init_constit_mrt_bingham_explicit(mu_p, tau_y, m, 1.0e-9, 1.0);
      const forcing_kf = init_korner_Fk(F);

    datadir:    { value: datadir,   expr: true    }

    # material init
    rho_0:      { value: 1.0,       expr: false   }
    nu:         { value: mu_p,      expr: true    }

    # lattice configuration
    dx:         { value: 1.0,       expr: false   }
    dt:         { value: 1.0,       expr: false   }
    ni:         { value: ni,        expr: true    }
    nj:         { value: nj,        expr: true    }

    # simulation parameters
    simtype:    default
    nsteps:     { value: nsteps,    expr: true    }
    col_f:      init_col_mrt(constit_rel_f, forcing_kf, S_luo)

    # boundaries
    sbounds:
      value: "[1 ni 1 nj;]'"
      expr: true

    cbounds:
      value: "[1 ni 1 nj;]'"
      expr: true

    # boundary conditions
    bcs:
      - north_bounce_back!
      - south_bounce_back!
      - periodic_east_to_west!

    # callback functions
    callbacks:
      - print_step_callback(100, id)
      - write_jld_file_callback(datadir, convert(Int64, round(nsteps/20)))

    # clean-up, backup, write out
    finally:
      - >
        (sim::Sim, k::Int) -> begin
          writedlm(joinpath(datadir, "ux_profile.dsv"), 
            extract_ux_prof_callback(convert(Int64, round(ni/2)))(sim), 
            ",");
        end
      - write_jld_file_callback(datadir)

    # test for conditions to end simulation
    test_for_term: is_steadystate_x

## TODO list
* [ ] implement free surface flow
* [x] break collision functions down into forcing method, constitutive eqn, etc.
* [ ] derive and implement more general velocity and pressure boundary conditions 
* [ ] create alt. stream/collide methods using grids of obstacle flags
* [x] write user documentation
* [ ] parallelize streaming and collision steps, multiscale mapping, and BCs
* [ ] wrap simulation code in a module
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
* [ ] consider "DSL" for mapping particle to macroscale and vice versa
