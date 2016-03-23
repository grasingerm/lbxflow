version: 0.2.6

preamble: >
  @init_plot_env();
  const datadir       =   joinpath("data","dam");
  const nu            =   0.2;
  const constit_rel_f =   init_constit_srt_const(nu);
  const F             =   [0.0; -1.0e-6];
  const forcing_kf    =   init_guo_Fk(F);
  const ni            =   50;
  const nj            =   50;

# data directory
datadir:    {   value: datadir,       expr: true    }

# material properties
rho_0:      {   value: 1.0,           expr: false   }
nu:         {   value: nu,            expr: true    }

# lattice parameters
dx:         {   value: 1.0,           expr: false   }
dt:         {   value: 1.0,           expr: false   }
ni:         {   value:  ni,           expr: true    }
nj:         {   value:  nj,           expr: true    }

# simulation parameters
simtype:    free_surface
col_f:      init_col_srt(constit_rel_f, forcing_kf);
nsteps:     {   value: 250000,        expr: false   }

# boundaries
sbounds:
  value:  "[1 ni 1 nj;]'"
  expr:   true

cbounds:
  value:  "[1 ni 1 nj;]'"
  expr:   true

# boundary conditions
bcs:
  - west_bounce_back!
  - east_bounce_back!
  - south_bounce_back!
  - north_bounce_back!

# free surface conditions
rho_g: 1.0

fill_x: { value: 0.5, expr: false }
fill_y: { value: 1.0, expr: false }

# callback functions
callbacks:
  - plot_mass_contours_callback(100, joinpath(datadir, "mass"))
  - print_step_callback(50, "free-surf")
#  - write_jld_file_callback(datadir, 500)

# clean-up, backup, write out
#finally:
#  - write_jld_file_callback(datadir)
