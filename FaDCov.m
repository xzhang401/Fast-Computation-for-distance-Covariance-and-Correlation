function dCov = FaDCov(x,y)
% Fast computing for the distance covariance
% Inputs: x,y   observations of two univariate random variables
% Output: dCov  distance covariance

n = length(x); temp=[1:n]';
[vx, Ix0] = sort(x); Ix(Ix0) = temp; Ix=Ix'; 
[vy, Iy0] = sort(y); Iy(Iy0) = temp; Iy=Iy'; 
sx = cumsum(vx); 
sy = cumsum(vy); 
alphax = Ix-1; 
alphay = Iy-1; 
betax = sx(Ix) - vx(Ix); 
betay = sy(Iy) - vy(Iy); 
xdot = sum(x); 
ydot = sum(y); 

aidot = xdot + (2*alphax-n).*x - 2*betax; 
bidot = ydot + (2*alphay-n).*y - 2*betay; 
Sab = sum(aidot.*bidot); 

adotdot = 2*sum(alphax.*x) -2*sum(betax); 
bdotdot = 2*sum(alphay.*y) -2*sum(betay); 

gamma_1  = PartialSum2D(x,y, ones(n,1));
gamma_x  = PartialSum2D(x,y, x);
gamma_y  = PartialSum2D(x,y, y);
gamma_xy = PartialSum2D(x,y, x.*y);

aijbij = sum(x.*y.*gamma_1 + gamma_xy - x.*gamma_y - y.*gamma_x); 
dCov = aijbij/n/(n-3) - 2*Sab/n/(n-2)/(n-3) + adotdot*bdotdot/n/(n-1)/(n-2)/(n-3); 
