const __lbxio_root__ = dirname(@__FILE__);
require(abspath(joinpath(__lbxio_root__, "multiscale.jl")));

#! Create a callback function for writing to a delimited file
#!
#! \param pre Prefix of filename
#! \param stepout Number of steps in between writing
#! \param A Function for extracting a 2D array from the multiscale map
#! \param delim Delimiter to separate values with
function write_datafile_callback (pre::String, stepout::Int, A::Function,
  dir=".", delim=',')

  return (msm::MultiscaleMap, k::Int) -> begin
    if k % stepout == 0
      writedlm(joinpath(dir, pre*"_step-$k.dsv"), A(msm), delim);
    end
  end;

end

#! Extract velocity profile cut parallel to y-axis
function extract_prof_callback(i::Int)

  return (msm::MultiscaleMap) -> begin
    const nj = size(msm.u)[2];
    x = Array(Float64, (nj, 3));

    for j=1:nj
      x[j,:] = [j, msm.u[i,j,1], msm.u[i,j,2]];
    end

    return x;
  end;

end

#! Extract ux profile cut parallel to y-axis
function extract_ux_prof_callback(i::Int)

  return (msm::MultiscaleMap) -> begin
    const ni, nj = size(msm.u);
    const u = vec(msm.u[i,:,1]);

    x = linspace(-0.5, 0.5, nj);
    y = u;

    return [x y];
  end;

end

#! Extract u_bar profile cut parallel to y-axis
function extract_ubar_prof_callback(i::Int)

  return (msm::MultiscaleMap) -> begin
    const ni, nj = size(msm.u);
    const u = vec(msm.u[i,:,1]);

    x = linspace(-0.5, 0.5, nj);
    y = u / maximum(u);

    return [x y];
  end;

end

#! Create callback for reporting step
function print_step_callback(step::Int)
  return (msm::MultiscaleMap, k::Int) -> begin
    if k % step == 0
      println("step $k");
    end
  end
end

#! Plot x-component of velocity profile cut parallel to y-axis
function plot_ux_profile_callback(i::Int, iters_per_frame::Int,
  pause::FloatingPoint = 0.1)

  return (msm::MultiscaleMap, k::Int) -> begin
    if k % iters_per_frame == 0
      const nj = size(msm.u)[2];

      x = linspace(-0.5, 0.5, nj);
      y = vec(msm.u[i,:,1]);

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
  pause::FloatingPoint = 0.1)

  return (msm::MultiscaleMap, k::Int) -> begin
    if k % iters_per_frame == 0
      const nj = size(msm.u)[2];
      const u = vec(msm.u[i,:,1]);

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
  pause::FloatingPoint = 0.1)

  return (msm::MultiscaleMap, k::Int) -> begin
    if k % iters_per_frame == 0
      clf();
      contour(transpose(u_mag(msm)));
      sleep(pause);
    end
  end

end

#! Plot velocity vectors for the domain
function plot_uvecs_callback(iters_per_frame::Int, 
  pause::FloatingPoint = 0.1)

  return (msm::MultiscaleMap, k::Int) -> begin
    if k % iters_per_frame == 0
      clf();
      quiver(msm.u[:,:,1], msm.u[:,:,2]);
      sleep(pause);
    end
  end

end

#! Plot x-component of velocity profile cut parallel to y-axis
function plot_ux_profile_callback(i::Int, iters_per_frame::Int, fname::String,
  pause::FloatingPoint = 0.1)

  return (msm::MultiscaleMap, k::Int) -> begin
    if k % iters_per_frame == 0
      const nj = size(msm.u)[2];

      x = linspace(-0.5, 0.5, nj);
      y = vec(msm.u[i,:,1]);

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
  pause::FloatingPoint = 0.1)

  return (msm::MultiscaleMap, k::Int) -> begin
    if k % iters_per_frame == 0
      const nj = size(msm.u)[2];
      const u = vec(msm.u[i,:,1]);

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
  pause::FloatingPoint = 0.1)

  return (msm::MultiscaleMap, k::Int) -> begin
    if k % iters_per_frame == 0
      clf();
      contour(transpose(u_mag(msm)));
      savefig(fname*"_step-$k.png");
      sleep(pause);
    end
  end

end

#! Plot velocity vectors for the domain
function plot_uvecs_callback(iters_per_frame::Int, fname::String, 
  pause::FloatingPoint = 0.1)

  return (msm::MultiscaleMap, k::Int) -> begin
    if k % iters_per_frame == 0
      clf();
      quiver(msm.u[:,:,1], msm.u[:,:,2]);
      savefig(fname*"_step-$k.png");
      sleep(pause);
    end
  end

end

#=
#! Animate x-component of velocity profile cut parallel to y-axis
function animate_ux_profile_callback(i::Int, iters_per_frame::Int,
  pause::FloatingPoint = 0.1)

  return (msm::MultiscaleMap, k::Int) -> begin
    if k % iters_per_frame == 0
      const nj = size(msm.u)[2];

      x = linspace(-0.5, 0.5, nj);
      y = vec(msm.u[i,:,1]);

      plot(x,y);
      xlabel("x / width");
      ylabel("ux (lat / sec)");
      sleep(pause);
    end
  end

end

#! Animate nondimensional x-component of velocity profile cut parallel to y-axis
function animate_ubar_profile_callback(i::Int, iters_per_frame::Int,
  pause::FloatingPoint = 0.1)

  return (msm::MultiscaleMap, k::Int) -> begin
    if k % iters_per_frame == 0
      const nj = size(msm.u)[2];

      x = linspace(-0.5, 0.5, nj);
      y = vec(msm.u[i,:,1]) / max(msm.u[i,:,1]);

      plot(x,y);
      xlabel("x / width");
      ylabel("ux / u_max");
      sleep(pause);
    end
  end

end

#! Animate x-component of velocity profile cut parallel to y-axis
function animate_umag_contour_callback(iters_per_frame::Int,
  pause::FloatingPoint = 0.1)

  return (msm::MultiscaleMap, k::Int) -> begin
    if k % iters_per_frame == 0
      contour(u_mag(msm));
      sleep(pause);
    end
  end

end
=#
