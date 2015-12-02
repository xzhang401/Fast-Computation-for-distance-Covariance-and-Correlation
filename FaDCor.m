function dCor = FaDCor(x,y)
% Fast computing for the distance covariance
% Inputs: x,y   observations of two univariate random variables
% Output: dCor  distance coefficient

dcovXY = FaDCov(x,y); 
dcovX  = FaDCov(x,x); 
dcovY  = FaDCov(y,y);
if abs(dcovX * dcovY)<1e-10
    dCor = 0; 
else 
    dCor = sign(dcovXY)*sqrt(abs(dcovXY)/sqrt(dcovX*dcovY));
end