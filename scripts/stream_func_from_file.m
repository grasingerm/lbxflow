function phi, psi = stream_func_from_file(u_file, v_file)
  u = transpose(csvread(u_file));
  v = transpose(csvread(v_file));

  phi, psi = flowfun(u, v);
end
