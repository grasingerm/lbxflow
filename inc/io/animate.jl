# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

#! Create callback for pausing the simulation
function pause_sim_callback(step::Int)
  return (sim::AbstractSim, k::Int) -> begin
    if k % step == 0
      println("Press ENTER to continue...");
      readline(STDIN);
    end
  end
end

#! Create callback for reporting step
function print_step_callback(step::Int)
  return (sim::AbstractSim, k::Int) -> begin
    if k % step == 0
      println("step $k");
    end
  end
end

#! Create callback for reporting step
function print_step_callback(step::Int, name::AbstractString)
  return (sim::AbstractSim, k::Int) -> begin
    if k % step == 0
      println(name * ":\tstep $k");
   end
  end
end

#! Extract velocity profile cut parallel to y-axis
function extract_prof_callback(i::Int)

  return (sim::AbstractSim) -> begin
    const nj = size(sim.msm.u, 3);
    x = Array(Float64, (nj, 3));

    for j=1:nj
      x[j,:] = [j, sim.msm.u[1,i,j], sim.msm.u[2,i,j]];
    end

    return x;
  end;

end

#! Extract ux profile cut parallel to y-axis
function extract_ux_prof_callback(i::Int)

  return (sim::AbstractSim) -> begin
    const ni, nj = size(sim.msm.u, 2), size(sim.msm.u, 3);
    const u = vec(sim.msm.u[1,i,:]);

    x = linspace(-0.5, 0.5, nj);
    y = u;

    return [x y];
  end;

end

#! Extract u_bar profile cut parallel to y-axis
function extract_ubar_prof_callback(i::Int)

  return (sim::AbstractSim) -> begin
    const ni, nj = size(sim.msm.u, 2), size(sim.msm.u, 3);
    const u = vec(sim.msm.u[1,i,:]);

    x = linspace(-0.5, 0.5, nj);
    y = u / maximum(u);

    return [x y];
  end;

end

#! Initialize plotting environment
macro init_plot_env()
  return quote
    import PyPlot;
    PyPlot.ion();
  end
end

#############################################################################
###################### Base plots ###########################################
#############################################################################

#! Plot x-component of velocity profile cut parallel to y-axis
function plot_ux_profile_callback(i::Int, iters_per_frame::Int,
                                  pause::AbstractFloat = 0.025)

  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
      const nj = size(sim.msm.u, 3);

      x = linspace(-0.5, 0.5, nj);
      y = vec(sim.msm.u[1,i,:]);

      PyPlot.clf();
      PyPlot.plot(x,y);
      PyPlot.xlabel("x / width");
      PyPlot.ylabel("ux (lat / sec)");
      PyPlot.draw();
      PyPlot.pause(0.001);
      sleep(pause);
    end
  end

end

#! Plot nondimensional x-component of velocity profile cut parallel to y-axis
function plot_ubar_profile_callback(i::Int, iters_per_frame::Int,
                                    pause::AbstractFloat = 0.025)

  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
      const nj = size(sim.msm.u, 3);
      const u = vec(sim.msm.u[1,i,:]);

      x = linspace(-0.5, 0.5, nj);
      y = u / maximum(u);

      PyPlot.clf();
      PyPlot.plot(x,y);
      PyPlot.xlabel("x / width");
      PyPlot.ylabel("ux / u_max");
      PyPlot.draw();
      PyPlot.pause(0.001);
      sleep(pause);
    end
  end

end

#! Plot x-component of velocity profile cut parallel to y-axis
function plot_umag_contour_callback(iters_per_frame::Int,
                                    pause::AbstractFloat = 0.025)

  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
      PyPlot.clf();
      cs = PyPlot.contourf(transpose(u_mag(sim.msm)));
      PyPlot.colorbar(cs);
      PyPlot.draw();
      PyPlot.pause(0.001);
      sleep(pause);
    end
  end

end

#! Plot velocity vectors for the domain
function plot_uvecs_callback(iters_per_frame::Int,
                             pause::AbstractFloat = 0.025)

  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
      PyPlot.clf();
      PyPlot.quiver(transpose(sim.msm.u[1,:,:]), transpose(sim.msm.u[2,:,:]));
      PyPlot.draw();
      PyPlot.pause(0.001);
      sleep(pause);
    end
  end

end

#! Plot streamlines for the domain
function plot_streamlines_callback(iters_per_frame::Int,
                                   pause::AbstractFloat = 0.025)

  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
      const ni, nj = size(sim.msm.rho);
      x = collect(linspace(0.0, 1.0, ni));
      y = collect(linspace(0.0, 1.0, nj));

      PyPlot.clf();
      PyPlot.streamplot(x, y, transpose(reshape(sim.msm.u[1,:,:], (ni, nj))), 
                        transpose(reshape(sim.msm.u[2,:,:], (ni, nj))));
      PyPlot.ylim(0.0, 1.0);
      PyPlot.xlim(0.0, 1.0);
      PyPlot.draw();
      PyPlot.pause(0.001);
      sleep(pause);
    end
  end

end

#! Plot pressure contours for the domain
function plot_pressure_contours_callback(iters_per_frame::Int,
                                         pause::AbstractFloat = 0.025)

  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
      PyPlot.clf();
      PyPlot.contour(transpose(pmap(rho -> rho*sim.lat.cssq, sim.msm.rho)));
      PyPlot.draw();
      PyPlot.pause(0.001);
      sleep(pause);
    end
  end

end

#! Plot mass matrix for the domain
function plot_mass_contours_callback(iters_per_frame::Int,
                                     pause::AbstractFloat = 0.025)

  return (sim::FreeSurfSim, k::Int) -> begin
    if k % iters_per_frame == 0
      PyPlot.clf();
      PyPlot.contour(transpose(sim.tracker.M));
      PyPlot.draw();
      PyPlot.pause(0.001);
      sleep(pause);
    end
  end

end

#! Plot mass matrix for the domain
function plot_strain_rate_mrt_contours_callback(iters_per_frame::Int,
                                                pause::AbstractFloat = 0.025)
  const M = @DEFAULT_MRT_M();
  const iM = @DEFAULT_MRT_IM();
  return (sim::AbstractSim, k_iter::Int) -> begin
    if k_iter % iters_per_frame == 0
      sr = Array(Float64, ni, nj);
      for j=1:nj, i=1:ni
        Sij = S_luo(@nu(sim.msm.omega[i,j], sim.lat.cssq, sim.lat.dt),
                    sim.msm.rho[i,j], sim.lat.cssq, sim.lat.dt);
        feq = Array(Float64, sim.lat.n);
        for k=1:sim.lat.n
          feq[k] = feq_incomp(sim.lat, sim.msm.rho[i,j], sim.msm.u[:,i,j], k);
        end
        D = strain_rate_tensor(sim.lat, sim.msm.rho[i,j],
                               sim.lat.f[:,i,j] - feq, M, iM, Sij);
        sr[i,j] = @strain_rate(D);
      end
      PyPlot.clf();
      cs = PyPlot.contourf(transpose(sr));
      PyPlot.colorbar(cs);
      PyPlot.draw();
      PyPlot.pause(0.001);
      sleep(pause);
    end
  end

end

#############################################################################
###################### Plots saved to file ##################################
#############################################################################

#! Plot x-component of velocity profile cut parallel to y-axis
function plot_ux_profile_callback(i::Int, iters_per_frame::Int, 
                                  fname::AbstractString,
                                  pause::AbstractFloat = 0.025)

  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
      const nj = size(sim.msm.u, 3);

      x = linspace(-0.5, 0.5, nj);
      y = vec(sim.msm.u[1,i,:]);

      PyPlot.clf();
      PyPlot.plot(x,y);
      PyPlot.xlabel("x / width");
      PyPlot.ylabel("ux (lat / sec)");
      PyPlot.savefig(@sprintf("%s_step-%09d.png", fname, k));
      PyPlot.draw();
      PyPlot.pause(0.001);
      sleep(pause);
    end
  end

end

#! Plot nondimensional x-component of velocity profile cut parallel to y-axis
function plot_ubar_profile_callback(i::Int, iters_per_frame::Int, 
                                    fname::AbstractString,
                                    pause::AbstractFloat = 0.025)

  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
      const nj = size(sim.msm.u, 3);
      const u = vec(sim.msm.u[1,i,:]);

      x = linspace(-0.5, 0.5, nj);
      y = u / maximum(u);

      PyPlot.clf();
      PyPlot.plot(x,y);
      PyPlot.xlabel("x / width");
      PyPlot.ylabel("ux / u_max");
      PyPlot.savefig(@sprintf("%s_step-%09d.png", fname, k));
      PyPlot.draw();
      PyPlot.pause(0.001);
      sleep(pause);
    end
  end

end

#! Plot x-component of velocity profile cut parallel to y-axis
function plot_umag_contour_callback(iters_per_frame::Int, fname::AbstractString,
                                    pause::AbstractFloat = 0.025)

  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
      PyPlot.clf();
      cs = PyPlot.contourf(transpose(u_mag(sim.msm)));
      PyPlot.colorbar(cs);
      PyPlot.savefig(@sprintf("%s_step-%09d.png", fname, k));
      PyPlot.draw();
      PyPlot.pause(0.001);
      sleep(pause);
    end
  end

end

#! Plot velocity vectors for the domain
function plot_uvecs_callback(iters_per_frame::Int, fname::AbstractString,
                             pause::AbstractFloat = 0.025)

  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
      PyPlot.clf();
      PyPlot.quiver(transpose(sim.msm.u[1,:,:]), transpose(sim.msm.u[2,:,:]));
      PyPlot.savefig(@sprintf("%s_step-%09d.png", fname, k));
      PyPlot.draw();
      PyPlot.pause(0.001);
      sleep(pause);
    end
  end

end

#! Plot streamlines for the domain
function plot_streamlines_callback(iters_per_frame::Int, fname::AbstractString,
                                   pause::AbstractFloat = 0.025)

  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
      const ni, nj = size(sim.msm.rho);
      x = collect(linspace(0.0, 1.0, ni));
      y = collect(linspace(0.0, 1.0, nj));

      PyPlot.clf();
      PyPlot.streamplot(x, y, transpose(reshape(sim.msm.u[1,:,:], (ni, nj))), 
                        transpose(reshape(sim.msm.u[2,:,:], (ni, nj))));
      PyPlot.ylim(0.0, 1.0);
      PyPlot.xlim(0.0, 1.0);
      PyPlot.savefig(@sprintf("%s_step-%09d.png", fname, k));
      PyPlot.draw();
      PyPlot.pause(0.001);
      sleep(pause);
    end
  end

end

#! Plot pressure contours for the domain
function plot_pressure_contours_callback(iters_per_frame::Int, 
                                         fname::AbstractString,
                                         pause::AbstractFloat = 0.025)

  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
      PyPlot.clf();
      PyPlot.contour(transpose(pmap(rho -> rho*sim.lat.cssq, sim.msm.rho)));
      PyPlot.savefig(@sprintf("%s_step-%09d.png", fname, k));
      PyPlot.draw();
      PyPlot.pause(0.001);
      sleep(pause);
    end
  end

end

#! Plot mass matrix for the domain
function plot_mass_contours_callback(iters_per_frame::Int, 
                                     fname::AbstractString,
                                     pause::AbstractFloat = 0.025)

  return (sim::FreeSurfSim, k::Int) -> begin
    if k % iters_per_frame == 0
      PyPlot.clf();
      PyPlot.contour(transpose(sim.tracker.M));
      PyPlot.savefig(@sprintf("%s_step-%09d.png", fname, k));
      PyPlot.draw();
      PyPlot.pause(0.001);
      sleep(pause);
    end
  end

end

#! Plot mass matrix for the domain
function plot_strain_rate_mrt_contours_callback(iters_per_frame::Int,
                                                fname::AbstractString,
                                                pause::AbstractFloat = 0.025)
  const M = @DEFAULT_MRT_M();
  const iM = @DEFAULT_MRT_IM();
  return (sim::AbstractSim, k_iter::Int) -> begin
    if k_iter % iters_per_frame == 0
      sr = Array(Float64, ni, nj);
      for j=1:nj, i=1:ni
        Sij = S_luo(@nu(sim.msm.omega[i,j], sim.lat.cssq, sim.lat.dt),
                  sim.msm.rho[i,j], sim.lat.cssq, sim.lat.dt);
        feq = Array(Float64, sim.lat.n);
        for k=1:sim.lat.n
          feq[k] = feq_incomp(sim.lat, sim.msm.rho[i,j], sim.msm.u[:,i,j], k);
        end
        D = strain_rate_tensor(sim.lat, sim.msm.rho[i,j],
                               sim.lat.f[:,i,j] - feq, M, iM, Sij);
        sr[i,j] = @strain_rate(D);
      end
      PyPlot.clf();
      cs = PyPlot.contourf(transpose(sr));
      PyPlot.colorbar(cs);
      PyPlot.savefig(@sprintf("%s_step-%09d.png", fname, k));
      PyPlot.draw();
      PyPlot.pause(0.001);
      sleep(pause);
    end
  end

end

#! Plot mass matrix for the domain
function plot_is_yielded_mrt_contours_callback(iters_per_frame::Int,
                                               fname::AbstractString,
                                               gamma_min::AbstractFloat,
                                               pause::AbstractFloat = 0.025)
  const M = @DEFAULT_MRT_M();
  const iM = @DEFAULT_MRT_IM();
  return (sim::AbstractSim, k_iter::Int) -> begin
    if k_iter % iters_per_frame == 0
      sr = Array(Float64, ni, nj);
      for j=1:nj, i=1:ni
        Sij = S_luo(@nu(sim.msm.omega[i,j], sim.lat.cssq, sim.lat.dt),
                  sim.msm.rho[i,j], sim.lat.cssq, sim.lat.dt);
        feq = Array(Float64, sim.lat.n);
        for k=1:sim.lat.n
          feq[k] = feq_incomp(sim.lat, sim.msm.rho[i,j], sim.msm.u[:,i,j], k);
        end
        D = strain_rate_tensor(sim.lat, sim.msm.rho[i,j],
                               sim.lat.f[:,i,j] - feq, M, iM, Sij);
        sr[i,j] = (@strain_rate(D) > gamma_min) ? 1.0 : 0.0;
      end
      PyPlot.clf();
      PyPlot.contourf(transpose(sr), levels=[0.0, 1.0]);
      PyPlot.savefig(@sprintf("%s_step-%09d.png", fname, k));
      PyPlot.draw();
      PyPlot.pause(0.001);
      sleep(pause);
    end
  end

end

#############################################################################
###################### Plots  with timestep printed #########################
#############################################################################

#! Plot x-component of velocity profile cut parallel to y-axis
function plot_ux_profile_callback(i::Int, iters_per_frame::Int,
                                  xy::Tuple{Number,Number},
                                  pause::AbstractFloat = 0.025)

  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
      const nj = size(sim.msm.u, 3);

      x = linspace(-0.5, 0.5, nj);
      y = vec(sim.msm.u[1,i,:]);

      PyPlot.clf();
      PyPlot.plot(x,y);
      PyPlot.xlabel("x / width");
      PyPlot.ylabel("ux (lat / sec)");
      PyPlot.text(xy[1], xy[2], "step: $k");
      PyPlot.draw();
      PyPlot.pause(0.001);
      sleep(pause);
    end
  end

end

#! Plot nondimensional x-component of velocity profile cut parallel to y-axis
function plot_ubar_profile_callback(i::Int, iters_per_frame::Int,
                                    xy::Tuple{Number,Number},
                                    pause::AbstractFloat = 0.025)

  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
      const nj = size(sim.msm.u, 3);
      const u = vec(sim.msm.u[1,i,:]);

      x = linspace(-0.5, 0.5, nj);
      y = u / maximum(u);

      PyPlot.clf();
      PyPlot.plot(x,y);
      PyPlot.xlabel("x / width");
      PyPlot.ylabel("ux / u_max");
      PyPlot.text(xy[1], xy[2], "step: $k");
      PyPlot.draw();
      PyPlot.pause(0.001);
      sleep(pause);
    end
  end

end

#! Plot x-component of velocity profile cut parallel to y-axis
function plot_umag_contour_callback(iters_per_frame::Int, 
                                    xy::Tuple{Number,Number},
                                    pause::AbstractFloat = 0.025)

  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
      PyPlot.clf();
      cs = PyPlot.contourf(transpose(u_mag(sim.msm)));
      PyPlot.colorbar(cs);
      PyPlot.text(xy[1], xy[2], "step: $k");
      PyPlot.draw();
      PyPlot.pause(0.001);
      sleep(pause);
    end
  end

end

#! Plot velocity vectors for the domain
function plot_uvecs_callback(iters_per_frame::Int, xy::Tuple{Number,Number},
                             pause::AbstractFloat = 0.025)

  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
      PyPlot.clf();
      PyPlot.quiver(transpose(sim.msm.u[1,:,:]), transpose(sim.msm.u[2,:,:]));
      PyPlot.text(xy[1], xy[2], "step: $k");
      PyPlot.draw();
      PyPlot.pause(0.001);
      sleep(pause);
    end
  end

end

#! Plot streamlines for the domain
function plot_streamlines_callback(iters_per_frame::Int, 
                                   xy::Tuple{Number,Number},
                                   pause::AbstractFloat = 0.025)

  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
      const ni, nj = size(sim.msm.rho);
      x = collect(linspace(0.0, 1.0, ni));
      y = collect(linspace(0.0, 1.0, nj));

      PyPlot.clf();
      PyPlot.streamplot(x, y, transpose(reshape(sim.msm.u[1,:,:], (ni, nj))), 
                        transpose(reshape(sim.msm.u[2,:,:], (ni, nj))));
      PyPlot.ylim(0.0, 1.0);
      PyPlot.xlim(0.0, 1.0);
      PyPlot.text(xy[1], xy[2], "step: $k");
      PyPlot.draw();
      PyPlot.pause(0.001);
      sleep(pause);
    end
  end

end

#! Plot pressure contours for the domain
function plot_pressure_contours_callback(iters_per_frame::Int,
                                         xy::Tuple{Number, Number},
                                         pause::AbstractFloat = 0.025)

  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
      PyPlot.clf();
      PyPlot.contour(transpose(pmap(rho -> rho*sim.lat.cssq, sim.msm.rho)));
      PyPlot.text(xy[1], xy[2], "step: $k");
      PyPlot.draw();
      PyPlot.pause(0.001);
      sleep(pause);
    end
  end

end

#! Plot mass matrix for the domain
function plot_mass_contours_callback(iters_per_frame::Int,
                                     xy::Tuple{Number, Number},
                                     pause::AbstractFloat = 0.025)

  return (sim::FreeSurfSim, k::Int) -> begin
    if k % iters_per_frame == 0
      PyPlot.clf();
      PyPlot.contour(transpose(sim.tracker.M));
      PyPlot.text(xy[1], xy[2], "step: $k");
      PyPlot.draw();
      PyPlot.pause(0.001);
      sleep(pause);
    end
  end

end

#! Plot mass matrix for the domain
function plot_strain_rate_mrt_contours_callback(iters_per_frame::Int,
                                                xy::Tuple{Number, Number},
                                                pause::AbstractFloat = 0.025)
  const M = @DEFAULT_MRT_M();
  const iM = @DEFAULT_MRT_IM();
  return (sim::AbstractSim, k_iter::Int) -> begin
    if k_iter % iters_per_frame == 0
      sr = Array(Float64, ni, nj);
      for j=1:nj, i=1:ni
        Sij = S_luo(@nu(sim.msm.omega[i,j], sim.lat.cssq, sim.lat.dt),
                  sim.msm.rho[i,j], sim.lat.cssq, sim.lat.dt);
        feq = Array(Float64, sim.lat.n);
        for k=1:sim.lat.n
          feq[k] = feq_incomp(sim.lat, sim.msm.rho[i,j], sim.msm.u[:,i,j], k);
        end
        D = strain_rate_tensor(sim.lat, sim.msm.rho[i,j],
                               sim.lat.f[:,i,j] - feq, M, iM, Sij);
        sr[i,j] = @strain_rate(D);
      end
      PyPlot.clf();
      cs = PyPlot.contourf(transpose(sr));
      PyPlot.colorbar(cs);
      PyPlot.text(xy[1], xy[2], "step: $k_iter");
      PyPlot.draw();
      PyPlot.pause(0.001);
      sleep(pause);
    end
  end

end

#############################################################################
############ Plots saved to file with timestep printed ######################
#############################################################################

#! Plot x-component of velocity profile cut parallel to y-axis
function plot_ux_profile_callback(i::Int, iters_per_frame::Int,
                                  xy::Tuple{Number, Number},
                                  fname::AbstractString,
                                  pause::AbstractFloat = 0.025)

  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
      const nj = size(sim.msm.u, 3);

      x = linspace(-0.5, 0.5, nj);
      y = vec(sim.msm.u[1,i,:]);

      PyPlot.clf();
      PyPlot.plot(x,y);
      PyPlot.xlabel("x / width");
      PyPlot.ylabel("ux (lat / sec)");
      PyPlot.text(xy[1], xy[2], "step: $k");
      PyPlot.savefig(@sprintf("%s_step-%09d.png", fname, k));
      PyPlot.draw();
      PyPlot.pause(0.001);
      sleep(pause);
    end
  end

end

#! Plot nondimensional x-component of velocity profile cut parallel to y-axis
function plot_ubar_profile_callback(i::Int, iters_per_frame::Int,
                                    xy::Tuple{Number, Number},
                                    fname::AbstractString,
                                    pause::AbstractFloat = 0.025)

  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
      const nj = size(sim.msm.u, 3);
      const u = vec(sim.msm.u[1,i,:]);

      x = linspace(-0.5, 0.5, nj);
      y = u / maximum(u);

      PyPlot.clf();
      PyPlot.plot(x,y);
      PyPlot.xlabel("x / width");
      PyPlot.ylabel("ux / u_max");
      PyPlot.text(xy[1], xy[2], "step: $k");
      PyPlot.savefig(@sprintf("%s_step-%09d.png", fname, k));
      PyPlot.draw();
      PyPlot.pause(0.001);
      sleep(pause);
    end
  end

end

#! Plot x-component of velocity profile cut parallel to y-axis
function plot_umag_contour_callback(iters_per_frame::Int,
                                    xy::Tuple{Number, Number},
                                    fname::AbstractString, 
                                    pause::AbstractFloat = 0.025)

  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
      PyPlot.clf();
      cs = PyPlot.contourf(transpose(u_mag(sim.msm)));
      PyPlot.colorbar(cs);
      PyPlot.text(xy[1], xy[2], "step: $k");
      PyPlot.savefig(@sprintf("%s_step-%09d.png", fname, k));
      PyPlot.draw();
      PyPlot.pause(0.001);
      sleep(pause);
    end
  end

end

#! Plot velocity vectors for the domain
function plot_uvecs_callback(iters_per_frame::Int,
                             xy::Tuple{Number, Number},
                             fname::AbstractString, 
                             pause::AbstractFloat = 0.025)

  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
      PyPlot.clf();
      PyPlot.quiver(transpose(sim.msm.u[1,:,:]), transpose(sim.msm.u[2,:,:]));
      PyPlot.text(xy[1], xy[2], "step: $k");
      PyPlot.savefig(@sprintf("%s_step-%09d.png", fname, k));
      PyPlot.draw();
      PyPlot.pause(0.001);
      sleep(pause);
    end
  end

end

#! Plot streamlines for the domain
function plot_streamlines_callback(iters_per_frame::Int,
                                   xy::Tuple{Number, Number},
                                   fname::AbstractString,
                                   pause::AbstractFloat = 0.025)

  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
      const ni, nj = size(sim.msm.rho);
      x = collect(linspace(0.0, 1.0, ni));
      y = collect(linspace(0.0, 1.0, nj));

      PyPlot.clf();
      PyPlot.streamplot(x, y, transpose(reshape(sim.msm.u[1,:,:], (ni, nj))), 
                        transpose(reshape(sim.msm.u[2,:,:], (ni, nj))));
      PyPlot.ylim(0.0, 1.0);
      PyPlot.xlim(0.0, 1.0);
      PyPlot.text(xy[1], xy[2], "step: $k");
      PyPlot.savefig(@sprintf("%s_step-%09d.png", fname, k));
      PyPlot.draw();
      PyPlot.pause(0.001);
      sleep(pause);
    end
  end

end

#! Plot pressure contours for the domain
function plot_pressure_contours_callback(iters_per_frame::Int,
                                         xy::Tuple{Number, Number},
                                         fname::AbstractString,
                                         pause::AbstractFloat = 0.025)

  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
      PyPlot.clf();
      PyPlot.contour(transpose(pmap(rho -> rho*sim.lat.cssq, sim.msm.rho)));
      PyPlot.text(xy[1], xy[2], "step: $k");
      PyPlot.savefig(@sprintf("%s_step-%09d.png", fname, k));
      PyPlot.draw();
      PyPlot.pause(0.001);
      sleep(pause);
    end
  end

end

#! Plot mass matrix for the domain
function plot_mass_contours_callback(iters_per_frame::Int,
                                     xy::Tuple{Number, Number},
                                     fname::AbstractString,
                                     pause::AbstractFloat = 0.025)

  return (sim::FreeSurfSim, k::Int) -> begin
    if k % iters_per_frame == 0
      PyPlot.clf();
      PyPlot.contour(transpose(sim.tracker.M));
      PyPlot.text(xy[1], xy[2], "step: $k");
      PyPlot.savefig(@sprintf("%s_step-%09d.png", fname, k));
      PyPlot.draw();
      PyPlot.pause(0.001);
      sleep(pause);
    end
  end

end

#! Plot mass matrix for the domain
function plot_strain_rate_mrt_contours_callback(iters_per_frame::Int,
                                                xy::Tuple{Number, Number},
                                                fname::AbstractString,
                                                pause::AbstractFloat = 0.025)
  const M = @DEFAULT_MRT_M();
  const iM = @DEFAULT_MRT_IM();
  return (sim::AbstractSim, k_iter::Int) -> begin
    if k_iter % iters_per_frame == 0
      sr = Array(Float64, ni, nj);
      for j=1:nj, i=1:ni
        Sij = S_luo(@nu(sim.msm.omega[i,j], sim.lat.cssq, sim.lat.dt),
                  sim.msm.rho[i,j], sim.lat.cssq, sim.lat.dt);
        feq = Array(Float64, sim.lat.n);
        for k=1:sim.lat.n
          feq[k] = feq_incomp(sim.lat, sim.msm.rho[i,j], sim.msm.u[:,i,j], k);
        end
        D = strain_rate_tensor(sim.lat, sim.msm.rho[i,j],
                               sim.lat.f[:,i,j] - feq, M, iM, Sij);
        sr[i,j] = @strain_rate(D);
      end
      PyPlot.clf();
      cs = PyPlot.contourf(transpose(sr));
      PyPlot.colorbar(cs);
      PyPlot.text(xy[1], xy[2], "step: $k_iter");
      PyPlot.savefig(@sprintf("%s_step-%09d.png", fname, k));
      PyPlot.draw();
      PyPlot.pause(0.001);
      sleep(pause);
    end
  end

end

#############################################################################
################## Contours with other options ##############################
#############################################################################

#! Plot x-component of velocity profile cut parallel to y-axis
function plot_umag_contour_callback(iters_per_frame::Int, fname::AbstractString,
                                    levs::Vector,
                                    pause::AbstractFloat = 0.025)

  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
      PyPlot.clf();
      cs = PyPlot.contourf(transpose(u_mag(sim.msm)), levels=levs);
      PyPlot.colorbar(cs);
      PyPlot.savefig(@sprintf("%s_step-%09d.png", fname, k));
      PyPlot.draw();
      PyPlot.pause(0.001);
      sleep(pause);
    end
  end

end

#! Plot pressure contours for the domain
function plot_pressure_contours_callback(iters_per_frame::Int, 
                                         fname::AbstractString,
                                         levs::Vector,
                                         pause::AbstractFloat = 0.025)


  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
      PyPlot.clf();
      cs = PyPlot.contour(transpose(pmap(rho -> rho*sim.lat.cssq, sim.msm.rho)),
                          levels=levs);
      PyPlot.colorbar(cs);
      PyPlot.savefig(@sprintf("%s_step-%09d.png", fname, k));
      PyPlot.draw();
      PyPlot.pause(0.001);
      sleep(pause);
    end
  end

end

#! Plot mass matrix for the domain
function plot_mass_contours_callback(iters_per_frame::Int, 
                                     fname::AbstractString,
                                     levs::Vector,
                                     pause::AbstractFloat = 0.025)

  return (sim::FreeSurfSim, k::Int) -> begin
    if k % iters_per_frame == 0
      PyPlot.clf();
      cs = PyPlot.contour(transpose(sim.tracker.M), levels=levs);
      PyPlot.colorbar(cs);
      PyPlot.savefig(@sprintf("%s_step-%09d.png", fname, k));
      PyPlot.draw();
      PyPlot.pause(0.001);
      sleep(pause);
    end
  end

end

#! Plot mass matrix for the domain
function plot_strain_rate_mrt_contours_callback(iters_per_frame::Int,
                                                fname::AbstractString,
                                                levs::Vector,
                                                pause::AbstractFloat = 0.025)
  const M = @DEFAULT_MRT_M();
  const iM = @DEFAULT_MRT_IM();
  return (sim::AbstractSim, k_iter::Int) -> begin
    if k_iter % iters_per_frame == 0
      sr = Array(Float64, ni, nj);
      for j=1:nj, i=1:ni
        Sij = S_luo(@nu(sim.msm.omega[i,j], sim.lat.cssq, sim.lat.dt),
                  sim.msm.rho[i,j], sim.lat.cssq, sim.lat.dt);
        feq = Array(Float64, sim.lat.n);
        for k=1:sim.lat.n
          feq[k] = feq_incomp(sim.lat, sim.msm.rho[i,j], sim.msm.u[:,i,j], k);
        end
        D = strain_rate_tensor(sim.lat, sim.msm.rho[i,j],
                               sim.lat.f[:,i,j] - feq, M, iM, Sij);
        sr[i,j] = @strain_rate(D);
      end
      PyPlot.clf();
      cs = PyPlot.contourf(transpose(sr), levels=levs);
      PyPlot.colorbar(cs);
      PyPlot.savefig(@sprintf("%s_step-%09d.png", fname, k));
      PyPlot.draw();
      PyPlot.pause(0.001);
      sleep(pause);
    end
  end

end

#! Plot x-component of velocity profile cut parallel to y-axis
function plot_umag_contour_callback(iters_per_frame::Int, fname::AbstractString,
                                    levs::Vector, rects::Vector,
                                    pause::AbstractFloat = 0.025)

  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
      PyPlot.clf();
      cs = PyPlot.contourf(transpose(u_mag(sim.msm)), levels=levs);
      PyPlot.colorbar(cs);
      for rect in rects
        PyPlot.axhspan(rect[1], rect[2], xmin=rect[3], xmax=rect[4], 
                       facecolor="black");
      end
      PyPlot.savefig(@sprintf("%s_step-%09d.png", fname, k));
      PyPlot.draw();
      PyPlot.pause(0.001);
      sleep(pause);
    end
  end

end

#! Plot pressure contours for the domain
function plot_pressure_contours_callback(iters_per_frame::Int, 
                                         fname::AbstractString,
                                         levs::Vector, rects::Vector,
                                         pause::AbstractFloat = 0.025)


  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
      PyPlot.clf();
      cs = PyPlot.contour(transpose(pmap(rho -> rho*sim.lat.cssq, sim.msm.rho)),
                          levels=levs);
      PyPlot.colorbar(cs);
      for rect in rects
        PyPlot.axhspan(rect[1], rect[2], xmin=rect[3], xmax=rect[4],
                       facecolor="black");
      end
      PyPlot.savefig(@sprintf("%s_step-%09d.png", fname, k));
      PyPlot.draw();
      PyPlot.pause(0.001);
      sleep(pause);
    end
  end

end

#! Plot mass matrix for the domain
function plot_mass_contours_callback(iters_per_frame::Int, 
                                     fname::AbstractString,
                                     levs::Vector, rects::Vector,
                                     pause::AbstractFloat = 0.025)

  return (sim::FreeSurfSim, k::Int) -> begin
    if k % iters_per_frame == 0
      PyPlot.clf();
      cs = PyPlot.contour(transpose(sim.tracker.M), levels=levs);
      PyPlot.colorbar(cs);
      for rect in rects
        PyPlot.axhspan(rect[1], rect[2], xmin=rect[3], xmax=rect[4],
                       facecolor="black");
      end
      PyPlot.savefig(@sprintf("%s_step-%09d.png", fname, k));
      PyPlot.draw();
      PyPlot.pause(0.001);
      sleep(pause);
    end
  end

end

#! Plot mass matrix for the domain
function plot_strain_rate_mrt_contours_callback(iters_per_frame::Int,
                                                fname::AbstractString,
                                                levs::Vector, rects::Vector,
                                                pause::AbstractFloat = 0.025)
  const M = @DEFAULT_MRT_M();
  const iM = @DEFAULT_MRT_IM();
  return (sim::AbstractSim, k_iter::Int) -> begin
    if k_iter % iters_per_frame == 0
      sr = Array(Float64, ni, nj);
      for j=1:nj, i=1:ni
        Sij = S_luo(@nu(sim.msm.omega[i,j], sim.lat.cssq, sim.lat.dt),
                    sim.msm.rho[i,j], sim.lat.cssq, sim.lat.dt);
        feq = Array(Float64, sim.lat.n);
        for k=1:sim.lat.n
          feq[k] = feq_incomp(sim.lat, sim.msm.rho[i,j], sim.msm.u[:,i,j], k);
        end
        D = strain_rate_tensor(sim.lat, sim.msm.rho[i,j],
                               sim.lat.f[:,i,j] - feq, M, iM, Sij);
        sr[i,j] = @strain_rate(D);
      end
      PyPlot.clf();
      cs = PyPlot.contourf(transpose(sr), levels=levs);
      PyPlot.colorbar(cs);
      for rect in rects
        PyPlot.axhspan(rect[1], rect[2], xmin=rect[3], xmax=rect[4],
                       facecolor="black");
      end
      PyPlot.savefig(@sprintf("%s_step-%09d.png", fname, k));
      PyPlot.draw();
      PyPlot.pause(0.001);
      sleep(pause);
    end
  end

end
