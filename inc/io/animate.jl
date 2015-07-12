const __animate_root__ = dirname(@__FILE__);
require(abspath(joinpath(__animate_root__, "..", "col", "mrt_matrices.jl")));

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
function print_step_callback(step::Int, name::String)
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

#! Plot x-component of velocity profile cut parallel to y-axis
function plot_ux_profile_callback(i::Int, iters_per_frame::Int,
                                  pause::FloatingPoint = 0.025)

  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
      const nj = size(sim.msm.u, 3);

      x = linspace(-0.5, 0.5, nj);
      y = vec(sim.msm.u[1,i,:]);

      clf();
      plot(x,y);
      xlabel("x / width");
      ylabel("ux (lat / sec)");
      sleep(pause);
    end
  end

end

#! Plot nondimensional x-component of velocity profile cut parallel to y-axis
function plot_ubar_profile_callback(i::Int, iters_per_frame::Int,
                                    pause::FloatingPoint = 0.025)

  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
      const nj = size(sim.msm.u, 3);
      const u = vec(sim.msm.u[1,i,:]);

      x = linspace(-0.5, 0.5, nj);
      y = u / maximum(u);

      clf();
      plot(x,y);
      xlabel("x / width");
      ylabel("ux / u_max");
      sleep(pause);
    end
  end

end

#! Plot x-component of velocity profile cut parallel to y-axis
function plot_umag_contour_callback(iters_per_frame::Int,
                                    pause::FloatingPoint = 0.025)

  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
      clf();
      cs = contourf(transpose(u_mag(sim.msm)));
      colorbar(cs);
      sleep(pause);
    end
  end

end

#! Plot velocity vectors for the domain
function plot_uvecs_callback(iters_per_frame::Int,
                             pause::FloatingPoint = 0.025)

  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
      clf();
      quiver(transpose(sim.msm.u[1,:,:]), transpose(sim.msm.u[2,:,:]));
      sleep(pause);
    end
  end

end

#! Plot streamlines for the domain
function plot_streamlines_callback(iters_per_frame::Int,
                                   pause::FloatingPoint = 0.025)

  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
      const ni, nj = size(sim.msm.rho);
      x = linspace(0.0, 1.0, ni);
      y = linspace(0.0, 1.0, nj);

      clf();
      streamplot(x, y, reshape(transpose(sim.msm.u[1,:,:]), (nj, ni)), 
                 reshape(transpose(sim.msm.u[2,:,:]), (nj, ni)));
      ylim(0.0, 1.0);
      xlim(0.0, 1.0);
      sleep(pause);
    end
  end

end

#! Plot pressure contours for the domain
function plot_pressure_contours_callback(iters_per_frame::Int,
                                         pause::FloatingPoint = 0.025)

  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
      clf();
      contour(transpose(pmap(rho -> rho*sim.lat.cssq, sim.msm.rho)));
      sleep(pause);
    end
  end

end

#! Plot mass matrix for the domain
function plot_mass_contours_callback(iters_per_frame::Int,
                                     pause::FloatingPoint = 0.025)

  return (sim::FreeSurfSim, k::Int) -> begin
    if k % iters_per_frame == 0
      clf();
      contour(transpose(sim.tracker.M));
      sleep(pause);
    end
  end

end

#! Plot mass matrix for the domain
function plot_strain_rate_mrt_contours_callback(iters_per_frame::Int,
                                                pause::FloatingPoint = 0.025)
  const M = @DEFAULT_MRT_M();
  const iM = @DEFAULT_MRT_IM();
  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
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
      clf();
      cs = contourf(transpose(sr));
      colorbar(cs);
      sleep(pause);
    end
  end

end

#! Plot x-component of velocity profile cut parallel to y-axis
function plot_ux_profile_callback(i::Int, iters_per_frame::Int, fname::String,
                                  pause::FloatingPoint = 0.025)

  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
      const nj = size(sim.msm.u, 3);

      x = linspace(-0.5, 0.5, nj);
      y = vec(sim.msm.u[1,i,:]);

      clf();
      plot(x,y);
      xlabel("x / width");
      ylabel("ux (lat / sec)");
      savefig(fname*"_step-$k.png");
      sleep(pause);
    end
  end

end

#! Plot nondimensional x-component of velocity profile cut parallel to y-axis
function plot_ubar_profile_callback(i::Int, iters_per_frame::Int, fname::String,
                                    pause::FloatingPoint = 0.025)

  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
      const nj = size(sim.msm.u, 3);
      const u = vec(sim.msm.u[1,i,:]);

      x = linspace(-0.5, 0.5, nj);
      y = u / maximum(u);

      clf();
      plot(x,y);
      xlabel("x / width");
      ylabel("ux / u_max");
      savefig(fname*"_step-$k.png");
      sleep(pause);
    end
  end

end

#! Plot x-component of velocity profile cut parallel to y-axis
function plot_umag_contour_callback(iters_per_frame::Int, fname::String,
                                    pause::FloatingPoint = 0.025)

  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
      clf();
      cs = contourf(transpose(u_mag(sim.msm)));
      colorbar(cs);
      savefig(fname*"_step-$k.png");
      sleep(pause);
    end
  end

end

#! Plot velocity vectors for the domain
function plot_uvecs_callback(iters_per_frame::Int, fname::String,
                             pause::FloatingPoint = 0.025)

  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
      clf();
      quiver(transpose(sim.msm.u[1,:,:]), transpose(sim.msm.u[2,:,:]));
      savefig(fname*"_step-$k.png");
      sleep(pause);
    end
  end

end

#! Plot streamlines for the domain
function plot_streamlines_callback(iters_per_frame::Int, fname::String,
                                   pause::FloatingPoint = 0.025)

  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
      const ni, nj = size(sim.msm.rho);
      x = linspace(0.0, 1.0, ni);
      y = linspace(0.0, 1.0, nj);

      clf();
      streamplot(x, y, reshape(transpose(sim.msm.u[1,:,:]), (nj, ni)), 
                 reshape(transpose(sim.msm.u[2,:,:]), (nj, ni)));
      ylim(0.0, 1.0);
      xlim(0.0, 1.0);
      savefig(fname*"_step-$k.png");
      sleep(pause);
    end
  end

end

#! Plot pressure contours for the domain
function plot_pressure_contours_callback(iters_per_frame::Int, fname::String,
                                         pause::FloatingPoint = 0.025)

  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
      clf();
      contour(transpose(pmap(rho -> rho*sim.lat.cssq, sim.msm.rho)));
      savefig(fname*"_step-$k.png");
      sleep(pause);
    end
  end

end

#! Plot mass matrix for the domain
function plot_mass_contours_callback(iters_per_frame::Int, fname::String,
                                     pause::FloatingPoint = 0.025)

  return (sim::FreeSurfSim, k::Int) -> begin
    if k % iters_per_frame == 0
      clf();
      contour(transpose(sim.tracker.M));
      savefig(fname*"_step-$k.png");
      sleep(pause);
    end
  end

end

#! Plot mass matrix for the domain
function plot_strain_rate_mrt_contours_callback(iters_per_frame::Int,
                                                fname::String,
                                                pause::FloatingPoint = 0.025)
  const M = @DEFAULT_MRT_M();
  const iM = @DEFAULT_MRT_IM();
  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
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
      clf();
      cs = contourf(transpose(sr));
      colorbar(cs);
      savefig(fname*"_step-$k.png");
      sleep(pause);
    end
  end

end

#! Plot mass matrix for the domain
function plot_is_yielded_mrt_contours_callback(iters_per_frame::Int,
                                               fname::String,
                                               pause::FloatingPoint = 0.025)
  const M = @DEFAULT_MRT_M();
  const iM = @DEFAULT_MRT_IM();
  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
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
        sr[i,j] = (@strain_rate(D) > 1e-9) ? 1.0 : 0.0;
      end
      clf();
      contourf(transpose(sr), levels=[0.0, 1.0]);
      savefig(fname*"_step-$k.png");
      sleep(pause);
    end
  end

end


#! Plot x-component of velocity profile cut parallel to y-axis
function plot_ux_profile_callback(i::Int, iters_per_frame::Int,
                                  xy::(Number,Number),
                                  pause::FloatingPoint = 0.025)

  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
      const nj = size(sim.msm.u, 3);

      x = linspace(-0.5, 0.5, nj);
      y = vec(sim.msm.u[1,i,:]);

      clf();
      plot(x,y);
      xlabel("x / width");
      ylabel("ux (lat / sec)");
      text(xy[1], xy[2], "step: $k");
      sleep(pause);
    end
  end

end

#! Plot nondimensional x-component of velocity profile cut parallel to y-axis
function plot_ubar_profile_callback(i::Int, iters_per_frame::Int,
                                    xy::(Number,Number),
                                    pause::FloatingPoint = 0.025)

  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
      const nj = size(sim.msm.u, 3);
      const u = vec(sim.msm.u[1,i,:]);

      x = linspace(-0.5, 0.5, nj);
      y = u / maximum(u);

      clf();
      plot(x,y);
      xlabel("x / width");
      ylabel("ux / u_max");
      text(xy[1], xy[2], "step: $k");
      sleep(pause);
    end
  end

end

#! Plot x-component of velocity profile cut parallel to y-axis
function plot_umag_contour_callback(iters_per_frame::Int, xy::(Number,Number),
                                    pause::FloatingPoint = 0.025)

  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
      clf();
      cs = contourf(transpose(u_mag(sim.msm)));
      colorbar(cs);
      text(xy[1], xy[2], "step: $k");
      sleep(pause);
    end
  end

end

#! Plot velocity vectors for the domain
function plot_uvecs_callback(iters_per_frame::Int, xy::(Number,Number),
                             pause::FloatingPoint = 0.025)

  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
      clf();
      quiver(transpose(sim.msm.u[1,:,:]), transpose(sim.msm.u[2,:,:]));
      text(xy[1], xy[2], "step: $k");
      sleep(pause);
    end
  end

end

#! Plot streamlines for the domain
function plot_streamlines_callback(iters_per_frame::Int, xy::(Number,Number),
                                   pause::FloatingPoint = 0.025)

  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
      const ni, nj = size(sim.msm.rho);
      x = linspace(0.0, 1.0, ni);
      y = linspace(0.0, 1.0, nj);

      clf();
      streamplot(x, y, reshape(transpose(sim.msm.u[1,:,:]), (nj, ni)), 
                 reshape(transpose(sim.msm.u[2,:,:]), (nj, ni)));
      ylim(0.0, 1.0);
      xlim(0.0, 1.0);
      text(xy[1], xy[2], "step: $k");
      sleep(pause);
    end
  end

end

#! Plot pressure contours for the domain
function plot_pressure_contours_callback(iters_per_frame::Int,
                                         xy::(Number, Number),
                                         pause::FloatingPoint = 0.025)

  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
      clf();
      contour(transpose(pmap(rho -> rho*sim.lat.cssq, sim.msm.rho)));
      text(xy[1], xy[2], "step: $k");
      sleep(pause);
    end
  end

end

#! Plot mass matrix for the domain
function plot_mass_contours_callback(iters_per_frame::Int,
                                     xy::(Number, Number),
                                     pause::FloatingPoint = 0.025)

  return (sim::FreeSurfSim, k::Int) -> begin
    if k % iters_per_frame == 0
      clf();
      contour(transpose(sim.tracker.M));
      text(xy[1], xy[2], "step: $k");
      sleep(pause);
    end
  end

end

#! Plot mass matrix for the domain
function plot_strain_rate_mrt_contours_callback(iters_per_frame::Int,
                                                xy::(Number, Number),
                                                pause::FloatingPoint = 0.025)
  const M = @DEFAULT_MRT_M();
  const iM = @DEFAULT_MRT_IM();
  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
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
      clf();
      cs = contourf(transpose(sr));
      colorbar(cs);
      text(xy[1], xy[2], "step: $k");
      sleep(pause);
    end
  end

end

#! Plot x-component of velocity profile cut parallel to y-axis
function plot_ux_profile_callback(i::Int, iters_per_frame::Int,
                                  xy::(Number,Number), fname::String,
                                  pause::FloatingPoint = 0.025)

  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
      const nj = size(sim.msm.u, 3);

      x = linspace(-0.5, 0.5, nj);
      y = vec(sim.msm.u[1,i,:]);

      clf();
      plot(x,y);
      xlabel("x / width");
      ylabel("ux (lat / sec)");
      text(xy[1], xy[2], "step: $k");
      savefig(fname*"_step-$k.png");
      sleep(pause);
    end
  end

end

#! Plot nondimensional x-component of velocity profile cut parallel to y-axis
function plot_ubar_profile_callback(i::Int, iters_per_frame::Int,
                                    xy::(Number,Number), fname::String,
                                    pause::FloatingPoint = 0.025)

  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
      const nj = size(sim.msm.u, 3);
      const u = vec(sim.msm.u[1,i,:]);

      x = linspace(-0.5, 0.5, nj);
      y = u / maximum(u);

      clf();
      plot(x,y);
      xlabel("x / width");
      ylabel("ux / u_max");
      text(xy[1], xy[2], "step: $k");
      savefig(fname*"_step-$k.png");
      sleep(pause);
    end
  end

end

#! Plot x-component of velocity profile cut parallel to y-axis
function plot_umag_contour_callback(iters_per_frame::Int, xy::(Number,Number),
                                    fname::String, pause::FloatingPoint = 0.025)

  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
      clf();
      cs = contourf(transpose(u_mag(sim.msm)));
      colorbar(cs);
      text(xy[1], xy[2], "step: $k");
      savefig(fname*"_step-$k.png");
      sleep(pause);
    end
  end

end

#! Plot velocity vectors for the domain
function plot_uvecs_callback(iters_per_frame::Int, xy::(Number,Number),
                             fname::String, pause::FloatingPoint = 0.025)

  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
      clf();
      quiver(transpose(sim.msm.u[1,:,:]), transpose(sim.msm.u[2,:,:]));
      text(xy[1], xy[2], "step: $k");
      savefig(fname*"_step-$k.png");
      sleep(pause);
    end
  end

end

#! Plot streamlines for the domain
function plot_streamlines_callback(iters_per_frame::Int, xy::(Number,Number),
                                   fname::String, pause::FloatingPoint = 0.025)

  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
      const ni, nj = size(sim.msm.rho);
      x = linspace(0.0, 1.0, ni);
      y = linspace(0.0, 1.0, nj);

      clf();
      streamplot(x, y, reshape(transpose(sim.msm.u[1,:,:]), (nj, ni)), 
                 reshape(transpose(sim.msm.u[2,:,:]), (nj, ni)));
      ylim(0.0, 1.0);
      xlim(0.0, 1.0);
      text(xy[1], xy[2], "step: $k");
      savefig(fname*"_step-$k.png");
      sleep(pause);
    end
  end

end

#! Plot pressure contours for the domain
function plot_pressure_contours_callback(iters_per_frame::Int,
                                         xy::(Number, Number),
                                         fname::String,
                                         pause::FloatingPoint = 0.025)

  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
      clf();
      contour(transpose(pmap(rho -> rho*sim.lat.cssq, sim.msm.rho)));
      text(xy[1], xy[2], "step: $k");
      savefig(fname*"_step-$k.png");
      sleep(pause);
    end
  end

end

#! Plot mass matrix for the domain
function plot_mass_contours_callback(iters_per_frame::Int,
                                     xy::(Number, Number),
                                     fname::String,
                                     pause::FloatingPoint = 0.025)

  return (sim::FreeSurfSim, k::Int) -> begin
    if k % iters_per_frame == 0
      clf();
      contour(transpose(sim.tracker.M));
      text(xy[1], xy[2], "step: $k");
      savefig(fname*"_step-$k.png");
      sleep(pause);
    end
  end

end

#! Plot mass matrix for the domain
function plot_strain_rate_mrt_contours_callback(iters_per_frame::Int,
                                                xy::(Number, Number),
                                                fname::String,
                                                pause::FloatingPoint = 0.025)
  const M = @DEFAULT_MRT_M();
  const iM = @DEFAULT_MRT_IM();
  return (sim::AbstractSim, k::Int) -> begin
    if k % iters_per_frame == 0
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
      clf();
      cs = contourf(transpose(sr));
      colorbar(cs);
      text(xy[1], xy[2], "step: $k");
      savefig(fname*"_step-$k.png");
      sleep(pause);
    end
  end

end
