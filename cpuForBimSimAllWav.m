% Description:      -Script to compute the forward model for simulating extended source
%                       Mie scattering in spheres.
%                   -Computations are performed on CPU only.
% Input:            See the Setup section for the parameters that need to be chosen.
% Output:           Absorbance images of a sphere at one or more wavelengths.
% References:       BIM-Sim: Interactive Simulation of Broadband Imaging Using Mie Theory
% Authors:          S. Berisha
% Last modified:    11/24/15

%% Setup
clear
addpath(genpath('~/source/stim-matlab/'))
brewer = brewermap(1000);

s=rng;

%load the file containing the wavelength and refractive index of the
%material
load pmma.mat
numWav = size(material,1);
wavenumbers = material(:,1);
wavenumbers = 1e4./wavenumbers;

params.a = 6.5;                     %radius of the sphere
params.ps = [0 0 0];
params.pf=([0 0 0]);
params.res = 32;
params.samples=10;
params.orderEf=10;
params.numPS=10;
params.objectiveMin = .2;   %inner objective NA
params.objectiveMax=.62;    %outer objective NA
params.NA_in = .2;
params.NA_out = 0.62;


params.fov = round(params.a)*4; %field of view in micrometers


%specify padding
padding = 1;    %use padding = 1, do not padd =0; padding is used for performing more accurate computations for
%specify the size of the field plane in wavelength units (microns)
params.gridSize = round(params.fov/2)*(2*padding + 1);

%specify the spatial resolution of the field plane
params.simRes = params.res*(2*padding + 1);

params.E0 = 1;


% compute alpha1 and alpha2 from NA_in and NA_out, respectively
params.alpha1 = asin(params.NA_in); params.alpha2 = asin(params.NA_out);


%compute the amplitude that makes it through the condenser
params.subA = 2 * pi * params.E0 * ( (1 - cos(params.alpha2)) - (1 - cos(params.alpha1)) );


%generate grid points
gridPoints = (2*params.gridSize)*(0:params.simRes-1)/params.simRes - params.gridSize;

df = 1/(params.gridSize*2);

[iv, iu] = meshgrid(0:params.simRes-1, 0:params.simRes-1);
%[iv, iu] = meshgrid(gridPoints, gridPoints);
u=zeros(size(iu)); v = u;
idx = find(iu <= params.simRes/2);
u(idx) = iu(idx);
idx = find(iu > params.simRes/2);
u(idx) = (iu(idx) - params.simRes+1);
u=u.*df;
idx = find(iv <= params.simRes/2);
v(idx) = iv(idx);
idx = find(iv > params.simRes/2);
v(idx) = iv(idx) - params.simRes + 1;
v=v.*df;
fmag = sqrt(u.*u + v.*v);    %compute the magnitude of the frequencies

%amount of cropping needed after padding
cropSize = padding*params.res;


if padding==0
    params.startIdx=1;
    params.endIdx=params.simRes;
else
    %params.startIdx = round((params.simRes  - cropSize)/2);
    
    params.startIdx = round(params.simRes/2) - floor(cropSize/2);
    params.endIdx = params.startIdx + cropSize-1;
end


%this gives the kVec direction from the y axis
theta=1.5708;
phi=0;
%direction of the incident light
[x,y,z] = sph2cart(theta,phi,1);
lightDirection  = [x y z];

[x,z] = meshgrid(gridPoints, gridPoints); % field slice in the x z plane
y = ones(params.simRes,params.simRes)*(params.a);   %field plane y = 0

params.rVecs = zeros(params.simRes*params.simRes, 3);
params.rVecs(:,1) = x(:); params.rVecs(:,2) = y(:); params.rVecs(:,3) = z(:); %r value at each pixel position
params.psVecs = bsxfun(@minus, params.rVecs, params.ps);
normPMinPs = sqrt(sum(params.psVecs.^2,2));
params.psMask=reshape(normPMinPs,params.simRes, params.simRes);
alpha1 = asin(params.NA_in); alpha2 = asin(params.NA_out);
params.Pl_cosalpha1 = myLegendre(params.orderEf+1,cos(alpha1));
params.Pl_cosalpha2 = myLegendre(params.orderEf+1,cos(alpha2));


E_t = zeros(params.simRes, params.simRes);

Et_d = zeros(params.res, params.res);   %total field at the detector
Ef_d = Et_d;    %incident/focused field at the detector

A = (zeros(params.res,params.res,numWav));
absSpec = (zeros(numWav,1));

origRVecs(:,1) = x(:); origRVecs(:,2) = y(:); origRVecs(:,3) = z(:);
rVecs_ps = bsxfun(@minus, origRVecs,params.ps);
params.r_ps=reshape(sqrt(sum(rVecs_ps.^2,2)),params.simRes, params.simRes); %r value at each pixel position with respect to the center of the sphere


params.normPMinPs = bsxfun(@rdivide, rVecs_ps,  sqrt(sum(rVecs_ps.^2,2)));


allEtd = zeros(params.res, params.res,numWav);
allEfd = allEtd;

D_Et = zeros(params.res, params.res);
D_Ef = D_Et;





params

%%
%tic
display('cpu forward BimSim for all wavelengths ....')

h = waitbar(0, 'Per wavelength computation...');
%% get point sources at a particular radius
samples = 32;
fraction = 3;
[pf_x, pf_z] = pointSources(samples, fraction, params.numPS);
xlim([-params.gridSize params.gridSize])
ylim([-params.gridSize params.gridSize])

for i=1:numWav
    
    
    if i==1
        tic
    end
    
    params.n = material(i,2) + 1i*material(i,3);
    
    params.lambda = material(i,1);
    
    params.wavNum = 2*pi/params.lambda;        %wavenumber
    
    params.numOrd = computeN_l(params.a, params.lambda);
    
    %create a vector of orders [0 1 ... Nl]
    ordVec = ((0:params.numOrd)');
    
    %The scattered field for a single incident plane-wave k produced by
    %a sphere with radius r positioned at point pf
    
    %calculate the prefix term (2l + 1)*i^l
    twolp1 = 2.*ordVec+1;
    il = 1i.^ordVec;
    params.twolp1_il = twolp1.*il;
    
    
    
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
    hl_ka_p = (derivSphHan(params.numOrd, ka));
    %evaluate the hankel functions of the first kind at ka
    hl_ka = (shank1(params.numOrd, ka, 'multiple'));
    %compute the denominator for coefficients A and B
    denAB = jl_kna.*hl_ka_p - hl_ka.*jl_kna_p*params.n;
    
    
    %compute B
    params.B = params.twolp1_il.*(numB./denAB);
    params.B = reshape(params.B,[1 1 params.numOrd+1]);
    %calculate the numerator for the scattering coefficients A
    numA = jl_ka.*hl_ka_p - jl_ka_p.*hl_ka;
    %calculate the scattering coefficients A
    params.A = params.twolp1_il.*(numA./denAB);
    params.A = reshape(params.A,[1 1 params.numOrd+1]);
    
    
    objMinParam = params.NA_in/params.lambda; %min cut off of frequencies
    objMaxParam = params.NA_out/params.lambda; %max cut off of frequencies
    BPF = (fmag < objMinParam | fmag > objMaxParam);   %compute the band pass filter
    
    params.kVec=(lightDirection*params.wavNum);
    params.normKvec = params.kVec./params.wavNum;
    
    
    
    params.k_j = (monteCarlo(s,params.samples, params.kVec, params.NA_in, params.NA_out));
    
    %% compute E_f and E_t at pf and ps [0 0 0]
    params.rVecs = bsxfun(@minus, origRVecs,params.pf);
    %norm of the position vectors with respect to the focal point
    normPMinPf = sqrt(sum(params.rVecs.^2,2));
    
    params.r=reshape(normPMinPf,params.simRes, params.simRes); %r value at each pixel position
    
    
    %compute and display Ef
    % display('Compting E_f...')
    E_f = newComputeEf(params);
    % display('Done computing E_f.')
    
    %compute and display E_s and E_i and E_t
    %  display('Computing E_s, E_i, E_t...')
    [E_s, E_i] = computeScatteredFields(params);
    
    
    
    E_t(params.psMask<params.a) = E_i(params.psMask<params.a);
    E_t(params.psMask>=params.a) = E_f(params.psMask>=params.a) + E_s(params.psMask>=params.a);
    
    E_t(isnan(E_t)) = 0;    E_t(isinf(E_t))=0;%near field = out_n in cuda bimsim
    
    
    Et_d = (fft2(E_t));
    fftEf = fft2(E_f);
    
    Et_d(BPF) = 0;
    fftEf(BPF) = 0;
    
    iftEt_d = ifft2((Et_d));
    iftEf = ifft2(fftEf);
    
    %first crop the filtered near-field image of the source and scattered fields
    cropEt_d=iftEt_d(params.startIdx:params.endIdx,params.startIdx:params.endIdx);
    cropEt_d(isnan(cropEt_d))=0;     cropEt_d(isinf(cropEt_d))=0;
    cropEf_d= iftEf(params.startIdx:params.endIdx, params.startIdx:params.endIdx);
    cropEf_d(isnan(cropEf_d))=0;     cropEf_d(isinf(cropEf_d))=0;
    %integrate and resample
    D_Et = D_Et + (abs(cropEt_d)).^2;   %this is out_i in cuda bimsim -- the measured intesity
    D_Ef = D_Ef + (abs(cropEf_d)).^2;   %this out_inc in cuda bimsim -- the measured incident field
    E_t = arrayfun(@(x) 0, E_t);
    %%
    
    for p = 1:params.numPS
        params.pf(1) = pf_x(p); params.pf(3) = pf_z(p);
        
        params.rVecs = bsxfun(@minus, origRVecs,params.pf);
        %norm of the position vectors with respect to the focal point
        normPMinPf = sqrt(sum(params.rVecs.^2,2));
        
        params.r=reshape(normPMinPf,params.simRes, params.simRes); %r value at each pixel position
        
        
        %compute and display Ef
        % display('Compting E_f...')
        E_f = newComputeEf(params);
        % display('Done computing E_f.')
        
        %compute and display E_s and E_i and E_t
        %  display('Computing E_s, E_i, E_t...')
        [E_s, E_i] = computeScatteredFields(params);
        
        
        
        E_t(params.psMask<params.a) = E_i(params.psMask<params.a);
        E_t(params.psMask>=params.a) = E_f(params.psMask>=params.a) + E_s(params.psMask>=params.a);
        
        E_t(isnan(E_t)) = 0;    E_t(isinf(E_t))=0;%near field = out_n in cuda bimsim
        
        
        Et_d = (fft2(E_t));
        fftEf = fft2(E_f);
        
        Et_d(BPF) = 0;
        fftEf(BPF) = 0;
        
        iftEt_d = ifft2((Et_d));
        iftEf = ifft2(fftEf);
        
        %first crop the filtered near-field image of the source and scattered fields
        cropEt_d=iftEt_d(params.startIdx:params.endIdx,params.startIdx:params.endIdx);
        cropEt_d(isnan(cropEt_d))=0;     cropEt_d(isinf(cropEt_d))=0;
        cropEf_d= iftEf(params.startIdx:params.endIdx, params.startIdx:params.endIdx);
        cropEf_d(isnan(cropEf_d))=0;     cropEf_d(isinf(cropEf_d))=0;
        %integrate and resample
        D_Et = D_Et + (abs(cropEt_d)).^2;   %this is out_i in cuda bimsim -- the measured intesity
        D_Ef = D_Ef + (abs(cropEf_d)).^2;   %this out_inc in cuda bimsim -- the measured incident field
        
        E_t = arrayfun(@(x) 0, E_t);
    end
    
    
    %calculate absorbance
    A(:,:,i) = -log10((D_Et)./(D_Ef));
    absSpec(i) = -log10(sum(D_Et(:))/sum(D_Ef(:)));
    absSpec(i)
    allEtd(:,:,i) = D_Et;
    allEfd(:,:,i) = D_Ef;
    
    D_Et = arrayfun(@(x) 0, D_Et);
    D_Ef = arrayfun(@(x) 0, D_Ef);
    
    if i==1
        toc
    end
    
    waitbar(i/numWav, h, sprintf('wavelength %f',material(i)));
    
end
close(h)

%save output
coutput.E_f = E_f;
coutput.E_s = E_s;
coutput.E_i = E_i;
coutput.E_t = E_t;
coutput.Et_d = D_Et;
coutput.Ef_d = D_Ef;
coutput.A = A;
coutput.allEtd = allEtd;
coutput.allEfd = allEfd;

%%
%%
figure;ax=axes;

plot(wavenumbers(1:i),absSpec(1:i))
set(ax, 'Xdir','reverse')
%%

%show results in subplots
%showBimSim(1, E_f, E_s, E_i, E_t, Et_d, Ef_d, A)

display('Done.')
