# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

using Printf;
import PyPlot;

#! Create callback for pausing the simulation
function pause_sim_callback(stepout::Real)
  return (sim::AbstractSim, k::Real) -> begin
    if k % stepout < sim.Δt
      println("Press ENTER to continue...");
      readline(stdin);
    end
  end
end

#! Create callback for reporting step
function print_step_callback(stepout::Real)
  return (sim::AbstractSim, k::Real) -> begin
    if k % stepout < sim.Δt
      println("step $k");
    end
  end
end

#! Create callback for reporting step
function print_step_callback(stepout::Real, name::AbstractString)
  return (sim::AbstractSim, k::Real) -> begin
    if k % stepout < sim.Δt
      println(name * ":\tstep $k");
    end
  end
end

#! Call garbage collector callback
function gc_callback(stepout::Real)
  return (sim::AbstractSim, k::Real) -> begin
    if k % stepout < sim.Δt
      gc();
    end
  end
end

#! Initialize plotting environment
macro init_plot_env()
  println("init_plot_env() is deprecated");
  #=return quote
    import PyPlot;
  end=#
end

#! Change default figure size
macro change_default_figsize(w, h)
  return :(PyPlot.rc("figure", figsize=($w, $h)));
end

#! Create a pyplot callback function for visualizing the simulation
#!
#! \param     stepout             Time per frame
#! \param     accessor            Accessor function variable of interest
#! \param     showfig             Show figure flag
#! \param     fname               Filename, if empty do not save figure
#! \param     title               Figure title
#! \param     xlabel              x-axis label
#! \param     ylabel              y-axis label
#! \param     xlim                x-axis limits, if false use default limits
#! \param     ylim                y-axis limits, if false use default limits
#! \param     xlogscale           Flag for converting x-axis to logscale
#! \param     ylogscale           Flag for converting y-axis to logscale
#! \param     grid                Flag for adding grid lines
#! \param     xticks              x-axis tick marks, if false use default
#! \param     yticks              y-axis tick marks, if false use default
#! \param     rects               List of rectangles
#! \param     rcolor              Rectangle color
#! \return                        Anonymous callback function for visualization
function pyplot_callback(stepout::Real, accessor::LBXFunction; 
                         showfig::Bool=true, fname::AbstractString="",
                         title::AbstractString="",
                         xlabel::AbstractString="", ylabel::AbstractString="", 
                         xlim=false, ylim=false, xlogscale=false, 
                         ylogscale=false, grid::Bool=false, xticks=false, 
                         yticks=false, rects=false, 
                         rcolor::AbstractString="black")

  return (sim::AbstractSim, k::Real) -> begin
    if k % stepout < sim.Δt
      if showfig == true
        PyPlot.ion();
      else
        PyPlot.ioff();
      end;

      x, y = accessor(sim);

      PyPlot.clf();
      fig = PyPlot.plot(x, y);

      if title != ""
        PyPlot.title(title);
      end
      if xlabel != ""
        PyPlot.xlabel(xlabel);
      end;
      if ylabel != ""
        PyPlot.ylabel(ylabel);
      end;

      if xlim != false
        PyPlot.xlim(xlim);
      end;
      if ylim != false
        PyPlot.ylim(ylim);
      end;

      if xlogscale != false
        PyPlot.xscale("log");
      end;
      if ylogscale != false
        PyPlot.yscale("log");
      end;
      
      if grid != false
        PyPlot.grid();
      end;

      if xticks != false
        PyPlot.xticks(xticks);
      end;
      if yticks != false
        PyPlot.yticks(yticks);
      end;

      if rects != false
        for rect in rects
          PyPlot.axhspan(rect[1], rect[2], xmin=rect[3], xmax=rect[4], 
                         facecolor=rcolor);
        end
      end;

      PyPlot.draw();

      if fname != ""
        PyPlot.savefig(@sprintf("%s_step-%09d.png", fname, convert(Int, round(k))));
      end

      if showfig == true
        PyPlot.pause(0.00001);
      else
        fig = 0;
      end;
    end
  end

end

#! Create a pycontour callback function for visualizing the simulation
#!
#! \param     stepout             Time per frame
#! \param     accessor            Accessor function variable of interest
#! \param     showfig             Show figure flag
#! \param     filled              Flag to use filled contours
#! \param     colorbar            Flag for showing colorbar
#! \param     levels              Colorbar levels, if false use default
#! \param     fname               Filename, if empty do not save figure
#! \param     title               Figure title
#! \param     xlabel              x-axis label
#! \param     ylabel              y-axis label
#! \param     xlim                x-axis limits, if false use default limits
#! \param     ylim                y-axis limits, if false use default limits
#! \param     grid                Flag for adding grid lines
#! \param     xticks              x-axis tick marks, if false use default
#! \param     yticks              y-axis tick marks, if false use default
#! \param     rects               List of rectangles
#! \param     rcolor              Rectangle color
#! \return                        Anonymous callback function for visualization
function pycontour_callback(stepout::Real, accessor::LBXFunction; 
                            showfig::Bool=true, filled=true, colorbar=false, 
                            levels=false, fname::AbstractString="",
                            title::AbstractString="", xlabel::AbstractString="", 
                            ylabel::AbstractString="", xlim=false, ylim=false, 
                            grid::Bool=false, xticks=false, yticks=false,
                            rects=false, rcolor::AbstractString="black")

  return (sim::AbstractSim, k::Real) -> begin
    if k % stepout < sim.Δt
      if showfig == true
        PyPlot.ion();
      else
        PyPlot.ioff();
      end;

      mat = accessor(sim);

      PyPlot.clf();
      if filled
       if colorbar
         if levels != false
             cs = PyPlot.contourf(mat, levels=levels);
             PyPlot.colorbar(cs);
         else
             cs = PyPlot.contourf(mat);
             PyPlot.colorbar(cs);
         end
       else
         cs = PyPlot.contourf(mat)
       end
     else
       if colorbar
         if levels != false
             cs = PyPlot.contour(mat, levels=levels);
             PyPlot.colorbar(cs);
         else
             cs = PyPlot.contour(mat);
             PyPlot.colorbar(cs);
         end
       else
         cs = PyPlot.contour(mat)
       end
     end
      
      if title != ""
        PyPlot.title(title);
      end
      if xlabel != ""
        PyPlot.xlabel(xlabel);
      end;
      if ylabel != ""
        PyPlot.ylabel(ylabel);
      end;

      if xlim != false
        PyPlot.xlim(xlim);
      end;
      if ylim != false
        PyPlot.ylim(ylim);
      end

      if grid != false
        PyPlot.grid();
      end;

      if xticks != false
        PyPlot.xticks(xticks);
      end
      if yticks != false
        PyPlot.yticks(yticks);
      end;

      if rects != false
        for rect in rects
          PyPlot.axhspan(rect[1], rect[2], xmin=rect[3], xmax=rect[4], 
                         facecolor=rcolor);
        end
      end

      PyPlot.draw();

      if fname != ""
        PyPlot.savefig(@sprintf("%s_step-%09d.png", fname, convert(Int, round(k))));
      end

      if showfig == true
        PyPlot.pause(0.00001);
      else
        cs = 0;
      end

    end
  end
end

#! Create a pyquiver callback function for visualizing the simulation
#!
#! \param     stepout             Time per frame
#! \param     accessor            Accessor function variable of interest
#! \param     showfig             Show figure flag
#! \param     fname               Filename, if empty do not save figure
#! \param     title               Figure title
#! \param     xlabel              x-axis label
#! \param     ylabel              y-axis label
#! \param     xlim                x-axis limits, if false use default limits
#! \param     ylim                y-axis limits, if false use default limits
#! \param     grid                Flag for adding grid lines
#! \param     xticks              x-axis tick marks, if false use default
#! \param     yticks              y-axis tick marks, if false use default
#! \param     rects               List of rectangles
#! \param     rcolor              Rectangle color
#! \return                        Anonymous callback function for visualization
function pyquiver_callback(stepout::Real, accessor::LBXFunction; 
                           showfig::Bool=true, fname::AbstractString="",
                           title::AbstractString="",
                           xlabel::AbstractString="", ylabel::AbstractString="", 
                           xlim=false, ylim=false, grid::Bool=false, 
                           xticks=false, yticks=false, rects=false, 
                           rcolor::AbstractString="black")

  return (sim::AbstractSim, k::Real) -> begin
    if k % stepout < sim.Δt
      if showfig == true
             PyPlot.ion();
           else
             PyPlot.ioff()
           end;

      u, v = accessor(sim);

      PyPlot.clf();
      fig = PyPlot.quiver(u, v);
      
      if title != ""
             PyPlot.title(title);
           end
      if xlabel != ""
             PyPlot.xlabel(xlabel);
           end;
      if ylabel != ""
             PyPlot.ylabel(ylabel);
           end;

      if xlim != false
             PyPlot.xlim(xlim);
           end;
      if ylim != false
             PyPlot.ylim(ylim);
           end;

      if grid != false
             PyPlot.grid();
           end;

      if xticks != false
             PyPlot.xticks(xticks);
           end;
      if yticks != false
             PyPlot.yticks(yticks);
           end;

      if rects != false
              for rect in rects
                PyPlot.axhspan(rect[1], rect[2], xmin=rect[3], xmax=rect[4], 
                               facecolor=rcolor);
              end
          end;

      PyPlot.draw();

      if fname != ""
        PyPlot.savefig(@sprintf("%s_step-%09d.png", fname, convert(Int, round(k))));
      end

      if showfig == true
               PyPlot.pause(0.00001);
           else
             fig = 0;
       end
     end
   end
end

#! Create a pyplot callback function for visualizing the simulation
#!
#! \param     stepout             Time per frame
#! \param     accessor            Accessor function variable of interest
#! \param     showfig             Show figure flag
#! \param     fname               Filename, if empty do not save figure
#! \param     title               Figure title
#! \param     xlabel              x-axis label
#! \param     ylabel              y-axis label
#! \param     xlim                x-axis limits, if false use default limits
#! \param     ylim                y-axis limits, if false use default limits
#! \param     grid                Flag for adding grid lines
#! \param     xticks              x-axis tick marks, if false use default
#! \param     yticks              y-axis tick marks, if false use default
#! \param     rects               List of rectangles
#! \param     rcolor              Rectangle color
#! \return                        Anonymous callback function for visualization
function pystream_callback(stepout::Real, accessor::LBXFunction; 
                           showfig::Bool=true, fname::AbstractString="",
                           title::AbstractString="",
                           xlabel::AbstractString="", ylabel::AbstractString="", 
                           xlim=false, ylim=false, grid::Bool=false, 
                           xticks=false, yticks=false, rects=false, 
                           rcolor::AbstractString="black")

  return (sim::AbstractSim, k::Real) -> begin
    if k % stepout < sim.Δt
      if showfig == true
             PyPlot.ion();
           else
             PyPlot.ioff();
           end;

      x, y, u, v = accessor(sim);

      PyPlot.clf();
      fig = PyPlot.streamplot(x, y, u, v);
      
      if title != ""
             PyPlot.title(title);
           end
      if xlabel != ""
             PyPlot.xlabel(xlabel);
           end;
      if ylabel != ""
             PyPlot.ylabel(ylabel);
           end;

      if xlim != false
             PyPlot.xlim(xlim);
           end;
      if ylim != false
             PyPlot.ylim(ylim);
           end;

      if grid != false
             PyPlot.grid();
           end;

      if xticks != false
             PyPlot.xticks(xticks);
           end;
      if yticks != false
             PyPlot.yticks(yticks);
           end;

      if rects != false
              for rect in rects
                PyPlot.axhspan(rect[1], rect[2], xmin=rect[3], xmax=rect[4], 
                               facecolor=rcolor);
              end
          end;

      PyPlot.draw();

      if fname != ""
        PyPlot.savefig(@sprintf("%s_step-%09d.png", fname, convert(Int, round(k))));
      end

      if showfig == true
               PyPlot.pause(0.00001);
           else
             fig = 0;
       end
     end
   end

end
