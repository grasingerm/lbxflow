include("../inc/LBXFlow.jl");
using LBXFlow
using Gadfly
using DataFrames

τs  =   map(x -> 2^x-1, 0:4);
exs =   vcat([1], map(x -> (1-2.0^(-x)), 4:-1:1));
ν   =   1.0;
∇p  =   -3.0;
n   =   25;
xs  =   linspace(-0.5, 0.5, n);
df1 =   DataFrame(x = Float64[], y = Float64[], Series = AbstractString[]);
df2 =   DataFrame(x = Float64[], y = Float64[], Series = AbstractString[]);
df3 =   DataFrame(x = Float64[], y = Float64[], Series = AbstractString[]);
df4 =   DataFrame(x = Float64[], y = Float64[], Series = AbstractString[]);

for (τ, ex) in zip(τs, exs) 
  us  =   LBXFlow.analytical_poise_bingham(ν, τ, ∇p, n);
  df1 =   vcat(df1, DataFrame(x = xs, y = us, 
                              Series = fill("τ<sub>y</sub> = $(τ)", n)));
  df2 =   vcat(df2, DataFrame(x = xs, y = us / maximum(us), 
                              Series = fill("τ<sub>y</sub> = $(τ)", n)));
  us  =   LBXFlow.analytical_poise_power_law(ν, ex, ∇p / 10, n);
  df3 =   vcat(df3, DataFrame(x = xs, y = us, 
                              Series = fill("n = $(ex)", n)));
  df4 =   vcat(df4, DataFrame(x = xs, y = us / maximum(us), 
                              Series = fill("n = $(ex)", n)));
end

draw(PNG("bingham.png", 5.5inch, 3inch),
       plot(df1, x=:x, y=:y, color=:Series,
       Geom.line,
       Guide.XLabel("y / H"),
       Guide.YLabel("u (lat / sec)")
       ));

draw(PNG("bingham_normalized.png", 5.5inch, 3inch),
       plot(df2, x=:x, y=:y, color=:Series,
       Geom.line,
       Guide.XLabel("y / H"),
       Guide.YLabel("u (lat / sec)")
       ));

draw(PNG("power-law.png", 5.5inch, 3inch),
       plot(df3, x=:x, y=:y, color=:Series,
       Geom.line,
       Guide.XLabel("y / H"),
       Guide.YLabel("u (lat / sec)")
       ));

draw(PNG("power-law_normalized.png", 5.5inch, 3inch),
       plot(df4, x=:x, y=:y, color=:Series,
       Geom.line,
       Guide.XLabel("y / H"),
       Guide.YLabel("u (lat / sec)")
       ));
