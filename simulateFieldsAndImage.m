

%% setup
% single plane wave case
clear
addpath(genpath('~/source/stim-matlab/'))

params.a = 3;                     %radius of the sphere

params.fov = 50;
%specify padding
padding = 1;

%specify the size of the field plane in wavelength units (microns)
params.gridSize = params.fov*(2*padding + 1);

params.n=1.4;

%wavelength
params.lambda = 2.5;


params.res = 256;

%specify the spatial resolution of the field plane
params.simRes = params.res*(2*padding + 1);


%create a parameter structure for the simulation

params.pf = [0 0 0];
params.ps = [10 10 0];
params.orderEf=200;

theta=0; phi=0;
kdir = [phi theta 1]
params.E0 = 1;
params.NA_in = 0;
params.NA_out = 1;

params.samples=200;

params.numOrd = computeN_l(params.a, params.lambda);
%params.numOrd = 50;

params.wavNum = 2*pi/params.lambda;        %wavenumber

%direction of the incident light
[x,y,z] = sph2cart(theta,phi,1);

params.kVec=[x y z]*params.wavNum;

% get r and rVecs
%create a grid of points representing pixel positions in the field plane
gridPoints = linspace(-params.gridSize,params.gridSize,params.simRes);
[x,y] = meshgrid(gridPoints, gridPoints);
z = zeros(params.simRes,params.simRes);   %field plane z = 0

%convert the field plane pixel positions to cartesian coordinates
[theta, phi, r] = cart2sph(x,y,z);


%store the point positions in rVecs
params.rVecs = zeros(params.simRes*params.simRes, 3);
params.rVecs(:,1) = x(:); params.rVecs(:,2) = y(:); params.rVecs(:,3) = z(:);
params.rVecs = bsxfun(@minus, params.rVecs,params.ps);

%norm of the position vectors
normPMinPs = sqrt(sum(params.rVecs.^2,2));
params.r=reshape(normPMinPs,params.simRes, params.simRes);                         %r value at each pixel position


params
%% compute and display Ef
display('Compting E_f...')
E_f = newComputeEf(params);
display('Done computing E_f.')

brewer = brewermap(1000);
figure, imagesc((abs((E_f)))),title('abs(E_f)'), colorbar, axis image
colormap(brewer)

%% compute and display E_s and E_i and E_t
display('Computing E_s, E_i, E_t...')
[E_s, E_i] = computeEsEi(params);
brewer = brewermap(1000);

figure, imagesc((abs((E_s)))),title('E_s'), colorbar, axis image
colormap(brewer)

figure, imagesc((abs((E_i)))),title('E_s'), colorbar, axis image
colormap(brewer)

E_t = zeros(params.simRes, params.simRes);
E_t(params.r<params.a) = E_i(params.r<params.a);
E_t(params.r>=params.a) = E_f(params.r>=params.a) + E_s(params.r>=params.a);
figure, imagesc((abs((E_t)))),title('E_t'), colorbar, axis image
colormap(brewer)
display('Done.')
