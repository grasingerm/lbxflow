function  f = cumsimp(y)

% F = CUMSIMP(Y)    Simpson-rule column-wise cumulative summation.
%       Numerical approximation of a function F(x) such that 
%       Y(X) = dF/dX.  Each column of the input matrix Y represents
%       the value of the integrand  Y(X)  at equally spaced points
%       X = 0,1,...size(Y,1).
%       The output is a matrix  F of the same size as Y.
%       The first row of F is equal to zero and each following row
%       is the approximation of the integral of each column of matrix
%       Y up to the givem row.
%       CUMSIMP assumes continuity of each column of the function Y(X)
%       and uses Simpson rule summation.
%       Similar to the command F = CUMSUM(Y), exept for zero first
%       row and more accurate summation (under the assumption of
%       continuous integrand Y(X)).
% 
%    See also CUMSUM, SUM, TRAPZ, QUAD

%  Kirill K. Pankratov, March 7, 1994.

 % 3-points interpolation coefficients to midpoints.
 % Second-order polynomial (parabolic) interpolation coefficients
 % from  Xbasis = [0 1 2]  to  Xint = [.5 1.5]
c1 = 3/8; c2 = 6/8; c3 = -1/8;

 % Determine the size of the input and make column if vector
ist = 0;         % If to be transposed
lv = size(y,1);
if lv==1, ist = 1; y = y(:); lv = length(y); end
f = zeros(size(y));

 % If only 2 elements in columns - simple sum divided by 2
if lv==2
  f(2,:) = (y(1,:)+y(2))/2;
  if ist, f = f'; end   % Transpose output if necessary
  return
end

 % If more than two elements in columns - Simpson summation
num = 1:lv-2;
   % Interpolate values of Y to all midpoints
f(num+1,:) = c1*y(num,:)+c2*y(num+1,:)+c3*y(num+2,:);
f(num+2,:) = f(num+2,:)+c3*y(num,:)+c2*y(num+1,:)+c1*y(num+2,:);
f(2,:) = f(2,:)*2; f(lv,:) = f(lv,:)*2;
   % Now Simpson (1,4,1) rule
f(2:lv,:) = 2*f(2:lv,:)+y(1:lv-1,:)+y(2:lv,:);
f = cumsum(f)/6;  % Cumulative sum, 6 - denom. from the Simpson rule

if ist, f = f'; end     % Transpose output if necessary