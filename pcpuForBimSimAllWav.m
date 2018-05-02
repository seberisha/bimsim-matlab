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

load pmma.mat
numWav = size(material,1);

params.a = 6.5;                     %radius of the sphere
params.ps = [0 0 0];
params.pf=([0 0 0]);
params.res = 32;

params.fov = round(params.a)*4; %field of view in micrometers

params.samples=10;
params.orderEf=10;
params.numPS=1;

params.NA_in = .2; %inner and outer NA of the condenser
params.NA_out = 0.62;

theta=1.5708; phi=0; %direction of the incident light: theta (0 - 2pi), phi (0-pi), r=1


params.objectiveMin = .2;   %inner objective NA
params.objectiveMax=.62;    %outer objective NA


params.numPS=1; %number of point sources
padding = 1;    %use padding = 1, do not padd =0; padding is used for performing more accurate computations for
%measuring the field at the detector.
params.gridSize = round(params.fov/2)*(2*padding + 1);   %specify the size of the field plane in wavelength units (micrometers)
params.simRes = params.res*(2*padding + 1); %compute the spatial resolution of the field plane
params.E0 = 1;  %amplitude of the field
%and internal
%fields
% compute alpha1 and alpha2 from NA_in and NA_out, respectively
params.alpha1 = asin(params.NA_in); params.alpha2 = asin(params.NA_out);


%compute the amplitude that makes it through the condenser
params.subA = 2 * pi * params.E0 * ( (1 - cos(params.alpha2)) - (1 - cos(params.alpha1)) );
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

Et_d = zeros(params.res, params.res);   %total field at the detector
Ef_d = Et_d;    %incident/focused field at the detector
[x,y,z] = sph2cart(theta,phi,1);    %direction of the incident light in cartesian coordinates
lightDirection = [x, y, z];
%create a grid of points representing pixel positions in the field plane
gridPoints = linspace(-params.gridSize,params.gridSize,params.simRes);

[z,x] = meshgrid(gridPoints, gridPoints); % field slice in the x z plane
y = ones(params.simRes,params.simRes).*params.a;   %field plane y = 0


params.rVecs = zeros(params.simRes*params.simRes, 3);
params.rVecs(:,1) = x(:); params.rVecs(:,2) = y(:); params.rVecs(:,3) = z(:); %r value at each pixel position
params.psVecs = bsxfun(@minus, params.rVecs, params.ps);
normPMinPs = sqrt(sum(params.psVecs.^2,2));
params.psMask=reshape(normPMinPs,params.simRes, params.simRes);




params.displaySubplots=0;
alpha1 = asin(params.NA_in); alpha2 = asin(params.NA_out);
params.Pl_cosalpha1 = myLegendre(params.orderEf+1,cos(alpha1));
params.Pl_cosalpha2 = myLegendre(params.orderEf+1,cos(alpha2));


params

%%
tic
display('gpu forward BimSim for all wavelengths ....')

h = waitbar(0, 'Per wavelength computation...');

A = (zeros(params.res,params.res,numWav));
absSpec = (zeros(numWav,1));
        E_t = zeros(params.simRes, params.simRes);


for i=1:numWav
    
    
    tic
    params.n = material(i,2) + 1i*material(i,3);
    %params.n
    
    params.lambda = material(i,1);
    
    params.wavNum = 2*pi/params.lambda;        %wavenumber
    
    params.numOrd = computeN_l(params.a, params.lambda);
    
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
    
    BPF = pbandpass(1/(params.gridSize*2), params.simRes, params.NA_in, params.NA_out, params.lambda);
    
    
    params.kVec = lightDirection*params.wavNum;  %incident light vector scaled by the wavenumber
    
    %generate a 'params.samples' number of k vectors in the direction of
    %incident light using Monte Carlo simulation
    params.k_j = monteCarlo(params.samples, params.kVec, params.NA_in, params.NA_out);
    
    for p = 1:params.numPS
        
        params.rVecs(:,1) = x(:); params.rVecs(:,2) = y(:); params.rVecs(:,3) = z(:);
        params.rVecs = bsxfun(@minus, params.rVecs,params.pf);
        %norm of the position vectors with respect to the focal point
        normPMinPf = sqrt(sum(params.rVecs.^2,2));
        params.r=reshape(normPMinPf,params.simRes, params.simRes); %r value at each pixel position
        %compute and display Ef
        % display('Compting E_f...')
        E_f = newComputeEf(params);
        % display('Done computing E_f.')
        
        params.rVecs(:,1) = x(:); params.rVecs(:,2) = y(:); params.rVecs(:,3) = z(:);
        params.rVecs = bsxfun(@minus, params.rVecs,params.ps);
        normPMinPs = sqrt(sum(params.rVecs.^2,2));
        params.r=reshape(normPMinPs,params.simRes, params.simRes); %r value at each pixel position with respect to the center of the sphere
        %compute and display E_s and E_i and E_t
        %  display('Computing E_s, E_i, E_t...')
        [E_s, E_i] = computeScatteredFields(params);
        E_t(params.psMask<params.a) = E_i(params.psMask<params.a);
        E_t(params.psMask>=params.a) = E_f(params.psMask>=params.a) + E_s(params.psMask>=params.a);
        
        E_t(isnan(E_t)) = 0;    %near field = out_n in cuda bimsim
        % this is not correct if the slice goes through the sphere
        
        E_d = fft2(E_t); 
        E_d(BPF) = 0;
        iftE_d = ifft2((E_d));
        %first crop the filtered near-field image of the source and scattered fields
        cropEd=iftE_d(params.startIdx:params.endIdx,params.startIdx:params.endIdx);
        cropEd(isnan(cropEd))=0;     cropEd(isinf(cropEd))=0;   %far field  = out_f in bimsim
        cropEf=E_f(params.startIdx:params.endIdx, params.startIdx:params.endIdx);
        cropEf(isinf(cropEf))=0;     cropEf(isinf(cropEf))=0;
        %integrate and resample
        Et_d = Et_d + (abs(cropEd)).^2; %intensity image = out_i in cuda bimsim
        Ef_d = Ef_d + (abs(cropEf)).^2; %incident field image = out_inc in cuda bimsim
        %calculate absorbance
        
        
        %         params.pf = randn([1 3]);
        %         params.pf(2)=0;
        %         params.pf
        E_t = arrayfun(@(x) 0, E_t);
    end
    A(:,:,i) = -log10(Et_d./Ef_d);
    absSpec(i) = -log10(sum(Et_d(:))/sum(Ef_d(:)));
    
    Et_d = arrayfun(@(x) 0, Et_d);
    Ef_d = arrayfun(@(x) 0, Ef_d);
    waitbar(i/numWav, h, sprintf('wavelength %f',material(i)));
    
    toc
    i
end

%save output
coutput.E_f = E_f;
coutput.E_s = E_s;
coutput.E_i = E_i;
coutput.E_t = E_t;
coutput.Et_d = Et_d;
coutput.Ef_d = Ef_d;
coutput.A = A;

%show results in subplots
%showBimSim(1, E_f, E_s, E_i, E_t, Et_d, Ef_d, A)

display('Done.')
