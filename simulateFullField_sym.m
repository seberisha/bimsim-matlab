% Description:      -Script to compute the forward model for simulating extended source 
%                       Mie scattering in spheres.
%                   -Computations are performed on CPU only.
% Input:            See the Setup section for the parameters that need to be chosen.
% Output:           Absorbance images of a sphere at one or more wavelengths.
% References:       BIM-Sim: Interactive Simulation of Broadband Imaging Using Mie Theory
% Authors:          S. Berisha
% Last modified:    11/24/15

%% Setup
%clear
addpath(genpath('~/source/stim-matlab/'))
%create a parameter structure for the simulation

load pmma.mat
numWav = size(material,1);

params.a = 6.5;                     %radius of the sphere
params.ps = [0 0 0];
params.pf=([0 0 0]);
params.res = 32;

params.fov = round(params.a)*4; %field of view in micrometers


params.objectiveMin = .2;   %inner objective NA
params.objectiveMax=.62;    %outer objective NA

params.samples=10;    %number of samples used for the Monte Carlo simulation of the scattered and internal fields
params.orderEf=10; %number of orders used for computing the incident/focused field
params.numPS=1; %number of point sources 
params.a = 6.5; %radius of the sphere in micrometeres
padding = 1;    %use padding = 1, do not padd =0; padding is used for performing more accurate computations for
                %measuring the field at the detector.
params.gridSize = round(params.fov/2)*(2*padding + 1);   %specify the size of the field plane in wavelength units (micrometers)

%%%%%
params.res = 128;   %spatial resolution - for symmetry this is a vector, half the size of the original x dimension of the resolution
%%%%%

params.simRes = params.res*(2*padding + 1); %compute the spatial resolution of the field plane
params.ps = [0 0 0];    %vector representing center of sphere
params.E0 = 1;  %amplitude of the field
params.NA_in = .2; %inner and outer NA of the condenser
params.NA_out = 0.62;
params.wavNum = 2*pi/params.lambda; %wavenumber
params.numOrd = computeN_l(params.a, params.lambda);    %the maximum order required for convergence - this is used for computing the scattered
                                                        %and internal
                                                        %fields
% compute alpha1 and alpha2 from NA_in and NA_out, respectively
params.alpha1 = asin(params.NA_in); params.alpha2 = asin(params.NA_out);
%create a vector of orders [0 1 ... Nl]
ordVec = (0:params.numOrd)';
%calculate the prefix term (2l + 1)*i^l
twolp1 = 2.*ordVec+1;
il = 1i.^ordVec;
twolp1_il = twolp1.*il;
%compute the arguments needed to evaluate spherical bessel functions,
%hankel functions, and their derivatives
ka=params.wavNum*params.a;
kna = params.wavNum*params.n*params.a;
%evaluate the spherical bessel functions of the first kind at ka
jl_ka = sphbesselj(params.numOrd,ka,'multiple');
%evaluate the derivate of the spherical bessel functions of the first kind at kna
jl_kna_p = derivSphBes(params.numOrd, kna);
%evaluate the spherical bessel functions of the first kind at kna
jl_kna = sphbesselj(params.numOrd,kna,'multiple');
%evaluate the derivative of the spherical bessel functions of the first kind at ka
jl_ka_p = derivSphBes(params.numOrd, ka);
%compute the numerator for B coefficients
numB = jl_ka.*jl_kna_p.*params.n - jl_kna.*jl_ka_p;
%evaluate the derivative of the hankel functions of the first kind at ka
hl_ka_p = derivSphHan(params.numOrd, ka);
%evaluate the hankel functions of the first kind at ka
hl_ka = shank1(params.numOrd, ka, 'multiple');
%compute the denominator for coefficients A and B
denAB = jl_kna.*hl_ka_p - hl_ka.*jl_kna_p*params.n;
%compute B
params.B = twolp1_il.*(numB./denAB);
%calculate the numerator for the scattering coefficients A
numA = jl_ka.*hl_ka_p - jl_ka_p.*hl_ka;
%calculate the scattering coefficients A
params.A = twolp1_il.*(numA./denAB);
%compute the amplitude that makes it through the condenser
params.subA = 2 * pi * params.E0 * ( (1 - cos(params.alpha2)) - (1 - cos(params.alpha1)) );


% D = params.gridSize*2;  %used to compute delta d needed for computing the indices for the band pass filter
% % Set up range of variables.
% u = 0:(params.simRes-1);
% v = 0:(params.simRes-1);
% % Compute the indices for use in meshgrid
% idx = find(u > params.simRes/2);
% u(idx) = u(idx) - params.simRes;
% idy = find(v > params.simRes/2);
% v(idy) = v(idy) - params.simRes;
% % Compute the meshgrid arrays
% [v, u] = meshgrid(v, u);
% df = 1/D;
% u=u.*df; 
% v=v.*df;
% 
% params.fmag = sqrt(u.*u + v.*v);    %compute the magnitude of the frequencies



params.objMinParam = params.objectiveMin/params.lambda; %min cut off of frequencies
params.objMaxParam = params.objectiveMax/params.lambda; %max cut off of frequencies
params.BPF = fftshift((params.fmag < params.objMinParam | params.fmag > params.objMaxParam));   %compute the band pass filter
cropSize = padding*params.res;  %amount of cropping needed after padding

if padding==0
    %no cropping if no padding
    params.startIdx=1;
    params.endIdx=params.simRes;
else
    %crop the middle of the padded image
    params.startIdx = round((params.simRes  - cropSize)/2);
    params.endIdx = params.startIdx + cropSize-1;
end

%%%%
Et_d = zeros(params.res, params.res);   %total field at the detector - vector in case of symmetry
%%%%

Ef_d = Et_d;    %incident/focused field at the detector
theta=1.5708; phi=0; %direction of the incident light: theta (0 - 2pi), phi (0-pi), r=1
[x,y,z] = sph2cart(theta,phi,1);    %direction of the incident light in cartesian coordinates
params.kVec=[x y z]*params.wavNum;  %incident light vector scaled by the wavenumber

%%%%
%create a grid of points representing pixel positions in the field plane
gridPoints = linspace(-params.gridSize,params.gridSize,params.simRes); %in case of symmetry this goes from 0 to gridSize, instead of -gridSize to gridSize
%%%%%

%generate a 'params.samples' number of k vectors in the direction of
%incident light using Monte Carlo simulation
params.k_j = monteCarlo(params.samples, params.kVec, params.NA_in, params.NA_out);  

%%%%
x = gridPoints'; % field slice in the x plane -- in case of symmetry choose x to be between 0 and gridSize
y = zeros(params.simRes,1);   %field plane y = 0 --vector if symmetry
z = zeros(params.simRes,1);   %field plane z = 0 --vector if symmetry

params.rVecs = zeros(params.simRes, 3);
%%%%

params.rVecs(:,1) = x(:); params.rVecs(:,2) = y(:); params.rVecs(:,3) = z(:); %r value at each pixel position
params.psVecs = bsxfun(@minus, params.rVecs, params.ps);
normPMinPs = sqrt(sum(params.psVecs.^2,2));

%%%%
params.psMask=normPMinPs;   %used to mask out the values inside the sphere area -- vector if symmetry 
%%%%

params.pf=[0 0 0];
params.displaySubplots=0; 
alpha1 = asin(params.NA_in); alpha2 = asin(params.NA_out);
params.Pl_cosalpha1 = myLegendre(params.orderEf+1,cos(alpha1));
params.Pl_cosalpha2 = myLegendre(params.orderEf+1,cos(alpha2));

%create a grid of points representing pixel positions in the 2d field plane
%gridPoints = linspace(-params.gridSize, params.gridSize,128);

[x_2d,z_2d] = meshgrid(gridPoints, gridPoints); % 2d field slice in the x z plane
y_2d = zeros(params.simRes,params.simRes);   % 2d field plane y = 0
[~,~,r_2d] = cart2sph(x_2d,y_2d,z_2d);
halfIdx = params.simRes/2 + 1;
[~,~,r_1d] = cart2sph(gridPoints, zeros(size(gridPoints)), zeros(size(gridPoints)));

params

%%
tic 
display('simulating the full field....')

for p = 1:params.numPS
   
    params.rVecs(:,1) = x(:); params.rVecs(:,2) = y(:); params.rVecs(:,3) = z(:); 
    params.rVecs = bsxfun(@minus, params.rVecs,params.pf);
    %norm of the position vectors with respect to the focal point 
    normPMinPf = sqrt(sum(params.rVecs.^2,2));
    
    %%%%
    params.r= normPMinPf; %r value at each pixel position --vector if symmetry
    %%%%
    
    %compute and display Ef
    display('Compting E_f...')
    E_f = computeEf_sym(params);
        
    display('Done computing E_f.')

    params.rVecs(:,1) = x(:); params.rVecs(:,2) = y(:); params.rVecs(:,3) = z(:);
    params.rVecs = bsxfun(@minus, params.rVecs,params.ps);
    normPMinPs = sqrt(sum(params.rVecs.^2,2));
    
    %%%%
    params.r=normPMinPs; %r value at each pixel position with respect to the center of the sphere --vector if symmetry
    %%%%
    
    %compute and display E_s and E_i and E_t
    display('Computing E_s, E_i, E_t...')
    [E_s, E_i] = computeScatteredFields_sym(params);
    %%%%
    E_t = zeros(params.simRes, 1); %vector if symmetry
    %%%%
    
    E_t(params.psMask<params.a) = E_i(params.psMask<params.a);
    E_t(params.psMask>=params.a) = E_f(params.psMask>=params.a) + E_s(params.psMask>=params.a);
    
    
    
    %replace NaNs with zeros
    %note: If the input to fft contains any NaN elements, then all the output elements will be NaN
    E_t(isnan(E_t))=0;
    
    Et_interp = interp1(params.r(params.simRes/2+1:end)', E_t(params.simRes/2+1:end), r_2d);
    Ef_interp = interp1(params.r(params.simRes/2+1:end)', E_f(params.simRes/2+1:end), r_2d);
    
    Et_interp(isnan(Et_interp))=0;
        
    % this is not correct if the slice goes through the sphere
    E_d = fft2(Et_interp);
    E_d = fftshift((E_d));
    E_d(params.BPF) = 0;
    iftE_d = ifft2(ifftshift(E_d));
    %first crop the filtered near-field image of the source and scattered fields
    cropEd=iftE_d(params.startIdx:params.endIdx, params.startIdx:params.endIdx);
    cropEd(isnan(cropEd))=0;     cropEd(isinf(cropEd))=0;
    cropEf=Ef_interp(params.startIdx:params.endIdx, params.startIdx:params.endIdx);
    cropEf(isinf(cropEf))=0;     cropEf(isinf(cropEf))=0;
    %integrate and resample
    Et_d = Et_d + (abs(cropEd)).^2;
    Ef_d = Ef_d + (abs(cropEf)).^2;
    %calculate absorbance
    A = -log10(Et_d./Ef_d);
    %figure, plot(A)
    
    %Ai = interp1(r_1d(65:end)', A(65:end), r_2d);
    %figure, imagesc(Ai),axis image, colormap(brewer)
    
    params.pf = randn([1 3]);
    params.pf(2)=0;
    params.pf
end

%perform interpolation


toc

%save output 
coutput.E_f = E_f;
coutput.E_s = E_s;
coutput.E_i = E_i;
coutput.E_t = E_t;
coutput.Et_d = Et_d;
coutput.Ef_d = Ef_d;
coutput.A = A;

display('Done.')
