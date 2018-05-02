function samples = monteCarlo(s, N, k, NAin, NAout)
%   Description:    Create Monte-Carlo samples of a cassegrain objective by 
%                   performing uniform sampling of a sphere and projecting 
%                   these samples onto an inscribed sphere.
%
%   Input:  -N - number of samples
%           -k - incident light direction in cartesian coordinates, e.g. [0 0 1]
%           -NAin - internal obscuration NA, e.g. 0.2
%           -NAout - outer cassegrain NA, e.g. 0.63
%   Output: -samples - a 3xN matrix, it's columns are vector samples
%   

rng(s)
if (size(k,2)==1)
    k=k';
end
cos_angle = k*[0; 0; 1];
R = eye(3);
if(cos_angle ~= 1.0)
    
    axis = cross([0 0 1], k);
    axis = axis./(norm(axis));
    angle = acos(cos_angle);
    
    q = [cos(angle/2) sin(angle/2)*axis];
    R = quat2rotm(q);
end

%find the phi values associated with the cassegrain ring
inPhi = asin(NAin);
outPhi = asin(NAout);

%calculate the z-values associated with these angles
inZ = cos(inPhi);
outZ = cos(outPhi);

rangeZ = inZ - outZ;

%draw a distribution of random phi, z values
samples = zeros(3, N);
for i=1:N
    z = rand()* rangeZ + outZ;
    theta = rand() * 2 *pi;
    
    %calculate theta
    phi = acos(z);
    
    %this is phi if we use Matlab's sph2cart
    %phi = asin(z);
    %compute and store cartesian coordinates
    %[x, y, z] = sph2cart(theta, phi, 1);
    r=1;
    x = r * cos(theta) * sin(phi); y = (r* sin(theta) * sin(phi));
    z = (r * cos(phi));
    cart=[x; y; z];
    samples(:,i) = R *cart;
end
