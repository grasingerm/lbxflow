using Base.Test
include(joinpath("..", "inc","multiscale.jl"));

println("Testing @nu and @omega macros...");
for i=1:100000
  lat1 = LatticeD2Q9(1.0, 1.0, 10, 10);
  lat2 = LatticeD2Q4(1.0, 1.0, 10, 10);
  rand_cssq = rand(0.1:0.1:0.8);
  rand_dt = rand(0.1:0.1:2.0);

  omega = rand(1.0e-3:1.0e-3:2.0);

  nu = @nu(omega, lat1.cssq, lat1.dt);
  @test_approx_eq omega @omega(nu, lat1.cssq, lat1.dt)
  
  nu = @nu(omega, lat2.cssq, lat2.dt);
  @test_approx_eq omega @omega(nu, lat2.cssq, lat2.dt)

  nu = @nu(omega, rand_cssq, rand_dt);
  @test_approx_eq omega @omega(nu, rand_cssq, rand_dt)
  
  nu = rand(1.0e-3:1.0e-3:10.0);

  omega = @omega(nu, lat1.cssq, lat1.dt);
  @test_approx_eq nu @nu(omega, lat1.cssq, lat1.dt)
  
  omega = @omega(nu, lat2.cssq, lat2.dt);
  @test_approx_eq nu @nu(omega, lat2.cssq, lat2.dt)

  omega = @omega(nu, rand_cssq, rand_dt);
  @test_approx_eq nu @nu(omega, rand_cssq, rand_dt)

  for j=1:rand(1:1000)
    nu = @nu(omega, rand_cssq, rand_dt);
    omega = @omega(nu, rand_cssq, rand_dt);
  end
  @test_approx_eq nu @nu(omega, rand_cssq, rand_dt)
  @test_approx_eq omega @omega(nu, rand_cssq, rand_dt)
end
println("Tests passed.");
