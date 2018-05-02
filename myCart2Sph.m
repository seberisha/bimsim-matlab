
function [theta,phi,r] = myCart2Sph(x,y,z)

r = sqrt(x.*x + y.*y + z.*z);
theta = atan2(y,x);
phi=zeros(size(r));
phi = acos(y(r~=0)./r(r~=0));