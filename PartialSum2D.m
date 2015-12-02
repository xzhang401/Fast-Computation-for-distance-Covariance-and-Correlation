function gamma = PartialSum2D(x,y,c)
% A subroutine for FaDCor; compute a type of partial sum
% Inputs: x, y, c
% Output: gamma({c})
% See paper for details 

n = length(x); temp=[1:n]';
[vx, Ix0] = sort(x); Ix(Ix0) = temp; %Ix = order stat
x = x(Ix0); y = y(Ix0); c = c(Ix0); % so x is at increasing order 
[vy, Iy0] = sort(y); Iy(Iy0) = temp; y = Iy';      % y is a perm of {1,...,n}
sy = cumsum(c(Iy0)) - c(Iy0); 
sx = cumsum(c) - c; 
cdot = sum(c);

gamma1 = DyadUpdate_c(y,c);

gamma = cdot - c -2*sy(Iy) -2*sx + 4*gamma1; 
gamma = gamma(Ix); 