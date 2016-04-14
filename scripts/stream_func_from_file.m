function [phi, psi] = stream_func_from_file(u_file, v_file)
  u = csvread(u_file);
  v = csvread(v_file);

  [phi, psi] = flowfun(u, v);
end
