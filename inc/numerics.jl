#! Compute the double dot product of two matrices
function ddot(A::Array{Float64,2}, B::Array{Float64,2})
  sum = 0;
  for (a,b) in zip(A,B)
    sum += a*b;
  end
  return sum;
end
