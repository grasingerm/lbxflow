# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

#! Initialize plotting environment
macro init_gadfly_plot_env()
  return quote
    import Gadfly, DataFrames;
    using Colors;
  end
end

#! Calculate analytical solution for newtonain poiseuille flow assuming BBBC
#!
#! \param   mu        Dynamic viscosity
#! \param   pgrad     Pressure gradient
#! \param   nnodes    Number of nodes in the channel width
#! \return            Velocity profile array
function analytical_poise_newton(mu::Real, pgrad::Real, nnodes::Int)
  const h       =   (nnodes - 1) / 2.0;

  results       =   zeros(nnodes);
  for i = 1:nnodes
    const x       =   i - nnodes/2.0 - 0.5;
    results[i]    =   -1.0 / (2.0 * mu) * pgrad * (h^2 - x^2);
  end

  return results;
end

#! Calculate analytical solution for bingham poiseuille flow assuming BBBC
#!
#! \param   mu        Plastic viscosity
#! \param   tau       Yield stress
#! \param   pgrad     Pressure gradient
#! \param   nnodes    Number of nodes in the channel width
#! \return            Velocity profile array
function analytical_poise_bingham(mu::Real, tau::Real, pgrad::Real, nnodes::Int)
  const y_tau   =   -tau / pgrad;
  const h       =   (nnodes - 1) / 2.0;

  results       =   zeros(nnodes);
  for i = 1:nnodes
    const x       =   i - nnodes/2.0 - 0.5;
    if abs(x) <= y_tau
      results[i]    =   (-1.0 / (2.0 * mu) * pgrad * (h^2 - y_tau^2) 
                         - tau / mu * (h - y_tau));
    else
      results[i]    =   (-1.0 / (2.0 * mu) * pgrad * (h^2 - x^2) 
                         - tau / mu * (h - abs(x)));
    end
  end

  return results;
end

#! Calculate analytical solution for HB poiseuille flow assuming BBBC
#!
#! \param   k         Flow consistency index
#! \param   n         Power law index
#! \param   tau       Yield stress
#! \param   pgrad     Pressure gradient
#! \param   nnodes    Number of nodes in the channel width
#! \return            Velocity profile array
function analytical_poise_hb(k::AbstractFloat, n::Real, tau::Real, pgrad::Real, 
                             nnodes::Int)
  error("'analytical_poise_hb' not yet implemented");

  const y_tau   =   -tau / pgrad;
  const h       =   (nnodes - 1) / 2.0;

  results       =   zeros(nnodes);
  for i = 1:nnodes
    const x       =   i - nnodes/2.0 - 0.5;
    if abs(x) <= y_tau
      results[i]    =   (-1.0 / (2.0 * mu) * pgrad * (h^2 - y_tau^2) 
                         - tau / mu * (h - y_tau));
    else
      results[i]    =   (-1.0 / (2.0 * mu) * pgrad * (h^2 - x^2) 
                         - tau / mu * (h - abs(x)));
    end
  end

  return results;
end

#! Calculate analytical solution for power law poiseuille flow assuming BBBC
#!
#! \param   k         Flow consistency index
#! \param   n         Power law index
#! \param   pgrad     Pressure gradient
#! \param   nnodes    Number of nodes in the channel width
#! \return            Velocity profile array
function analytical_poise_power_law(k::AbstractFloat, n::Real, tau::Real, 
                                    pgrad::Real, nnodes::Int)
  const h       =   (nnodes - 1) / 2.0;
  const n_rat   =   convert(Float64, (n + 1.0) / n);

  results       =   zeros(nnodes);
  for i = 1:nnodes
    const x       =   i - nnodes/2.0 - 0.5;
    const l_n     =   1.0 / n_rat * (-pgrad / k);
    results[i]    =   l_n * (h^n_rat - abs(x)^n_rat);
  end

  return results;
end

#! Calculate Lp relative error
#!
#! \param   analyt    Analytical solution
#! \param   approx    Approximate solution
#! \param   p         P
#! \return            L2 relative error
macro Lp_rel_error(analyt, approx, p)
  return :(norm(analyt - approx, p) / norm(analyt, p));
end

#! Calculate error of lbm approximation with respect to an analytical solution
#!
#! \param   sim                   Simulation data structure
#! \param   analytical_solution   analytical_solution(nnodes)
#! \param   d                     Dimension of flow index
#! \param   idx                   Index at which to take cross-section
#! \param   ps                    Lp Ps
#! \return                        Lp errors
function lbm_error(sim::AbstractSim, analytical_solution::LBXFunction, d::Int,
                   idx::Int; ps::Vector=[2, Inf], plot_errors::Bool=false,
                   show_plot::Bool=false, basename::AbstractString="",
                   plot_size::Tuple{Int, Int}=(6, 6))
  @assert(size(sim.msm.u, 1) >= d, 
          "Dimension of direction is larger than dimension of lattice.");
  @assert(size(sim.msm.u, d+1) >= idx, 
          "Index in which to take cross-section is out of the domain.");
  gd = Gadfly;

  const nnodes  =   (d == 1) ? size(sim.msm.u, 3) : size(sim.msm.u, 2);
  const approx  =   ((d == 1) ? vec(sim.msm.u[d, idx, :]) : 
                                vec(sim.msm.u[d, :, idx]));
  const analyt  =   analytical_solution(nnodes);
  rerrors       =   [];
  for p in ps
    push!(rerrors, @Lp_rel_error(analyt, approx, p));
  end

  if plot_errors
    for (model_error, etype) in zip((analyt-approx, (analyt-approx)/maximum(analyt)),
                                    ("Absolute", "Relative"))
      eplt = gd.plot(x=linspace(-0.5, 0.5, nnodes), y=model_error, gd.Geom.line,
                     gd.Guide.XLabel("x"), gd.Guide.YLabel("$etype error"));
      gd.draw(gd.PDF("$(basename)_$(etype)_error.pdf", plot_size[1]gd.inch, 
              plot_size[2]gd.inch), eplt); 
      if show_plot
        eplt 
      end
    end
  end

  return rerrors;
end

#! Plot lbm approximation with an analytical solution
#!
#! \param   sim                   Simulation data structure
#! \param   analytical_solution   analytical_solution(nnodes)
#! \param   d                     Dimension of flow index
#! \param   idx                   Index at which to take cross-section
#! \param   show_plot             Switch for shwoing plot
#! \param   fname                 Filename to save plot to
#! \param   plot_size             Tuple of dimensions for plot (inches)
function plot_lbm_vs_analyt(sim::AbstractSim, analytical_solution::LBXFunction, 
                            d::Int, idx::Int; show_plot::Bool=false, 
                            fname::AbstractString="", 
                            plot_size::Tuple{Int, Int}=(6, 6))
  @assert(size(sim.msm.u, 1) >= d, 
          "Dimension of direction is larger than dimension of lattice.");
  @assert(size(sim.msm.u, d+1) >= idx, 
          "Index in which to take cross-section is out of the domain.");
  gd = Gadfly;

  const nnodes  =   (d == 1) ? size(sim.msm.u, 3) : size(sim.msm.u, 2);
  const approx  =   ((d == 1) ? vec(sim.msm.u[d, idx, :]) : 
                                vec(sim.msm.u[d, :, idx]));
  const analyt  =   analytical_solution(nnodes);

  const x       =   linspace(-0.5, 0.5, nnodes);
  const df1     =   DataFrames.DataFrame(x=x, y=approx, label="LBM"); 
  const df2     =   DataFrames.DataFrame(x=x, y=analyt, label="Analytical");
  const df      =   DataFrames.vcat(df1, df2);
  
  #uplt = gd.plot(df, x="x", y="y", color="label", gd.Geom.line,
  #               gd.Guide.XLabel("x (lat)"), gd.Guide.YLabel("u (lat/sec)"));
  uplt = gd.plot(gd.layer(x=x, y=approx, gd.Geom.point, gd.Theme(default_color=parse(Colors.Colorant, "blue"))),
                 gd.layer(x=x, y=analyt, gd.Geom.line, gd.Theme(default_color=parse(Colors.Colorant, "red"))),
                 gd.Guide.XLabel("x (lat)"), gd.Guide.YLabel("u (lat/sec)"),
                 gd.Guide.manual_color_key("Legend", ["LBM", "Analytical"],
                                           ["blue", "red"]));
  gd.draw(gd.PDF(fname, plot_size[1]gd.inch, plot_size[2]gd.inch), 
          uplt);

  if show_plot
    uplt
  end
end

#! Initialize function for reporting LBM error
function report_lbm_error(f_analyt::LBXFunction, d::Int, idx::Int, 
                          datadir::AbstractString)

  (sim::AbstractSim, k::Int) -> begin
    const errors            =   lbm_error(sim, f_analyt, d, idx; ps=[2, Inf], 
                                          plot_errors=true, 
                                          basename=joinpath(datadir, 
                                                            "error-plots"));
    
    info("Relative L2 error   = $(errors[1])");
    info("Relative LInf error = $(errors[2])");
    h = open(joinpath(datadir, "rerrors.txt"), "w");
    write(h, "Relative L2 error   = $(errors[1])\n");
    write(h, "Relative LInf error = $(errors[2])\n");
    close(h);
  end
end

#! Initialize function for plotting lbm approx against analytical solution
function init_plot_lbm_vs_analyt(f_analyt::LBXFunction, d::Int, idx::Int,
                                 datadir::AbstractString)

  return ((sim::AbstractSim, k::Int) -> begin;
    plot_lbm_vs_analyt(sim, f_analyt, d, 
                       idx; fname=joinpath(datadir, "velocity_profiles.pdf"));
  end);
end
