#! Compute the double dot product of two matrices
function ddot(A::Array{Number,2}, B::Array{Number,2})
  sum = 0;
  for (a,b) in zip(A,B)
    sum += a*b;
  end
  return sum;
end
