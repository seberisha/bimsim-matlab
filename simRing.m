function r = simRing(ring, uc, vc,  N, R, radius, D_E, D_x, D_y)

%
% Input:
%   - N = # of sample in the ring
%   - R - # of radial sample rings used for integration
%   - radius - length of the radius of samples
%   - D_E - field at the detector for a point source
%   - D_x, D_y - cartesian coordinates for the image at the detector
%
% Output:
%   - r - 1xR vector representing the contribution of the point source
%           at each ring


% store the detector image as a function of rho in spherical coordinates
% rho = 0 is the sphere center
r = zeros(1,R);

rho = linspace(0, radius, R);

% 1st ps at center (0,0)
% scale the value in the center by the number of ps samples since
% this is not accounted below
r(1) = interp2(D_x, D_y, D_E, uc, vc)*N(ring);

for i=2:R
    %N = ceil(2*pi*(rho(i)));
    theta = 2*pi/N(i):2*pi/N(i):2*pi;
    
    %convert to cartesian coordinates  
    [x,y] = pol2cart(theta,rho(i)); 
    x = x + uc;
    y = vc - y;
    temp = interp2(D_x, D_y, D_E, x, y);
    temp(isnan(temp)) = 0;
    r(i) = sum(temp);
end


