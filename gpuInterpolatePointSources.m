%% setup
clear

addpath(genpath('~/source/stim-matlab/'))
brewer = brewermap(1000);
s=rng;

load pmma.mat
numWav = size(material,1);
%numWav=1;

wavenumbers = material(:,1);
wavenumbers = 1e4./wavenumbers;

params.a = 6.5;                     %radius of the sphere
params.ps = [0 0 0];
params.pf=([0 0 0]);
params.res = 33;

params.fov = round(params.a)*4; %field of view in micrometers

params.samples=100;
params.orderEf=100;
params.numPS=1;

%specify padding
padding = 1;
%specify the size of the field plane in wavelength units (microns)
params.gridSize = round(params.fov/2)*(2*padding + 1);

%specify the spatial resolution of the field plane
params.simRes = params.res*(2*padding + 1);
%create a parameter structure for the simulation
params.gpu_ps = gpuArray(params.ps);

params.E0 = 1;
params.NA_in = .2;
params.NA_out = 0.62;


% compute alpha1 and alpha2 from NA_in and NA_out, respectively
params.alpha1 = asin(params.NA_in); params.alpha2 = asin(params.NA_out);



%calculate the prefix term (2l + 1)*i^l
ordVecEf=gpuArray((0:params.orderEf)');
params.il = 1i.^ordVecEf;
params.il = reshape(params.il,[1 1 params.orderEf+1]);


%compute the amplitude that makes it through the condenser
params.subA = 2 * pi * params.E0 * ( (1 - cos(params.alpha2)) - (1 - cos(params.alpha1)) );

% D = params.gridSize*2;
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
% params.objectiveMin = .2;
% params.objectiveMax=.62;
% params.fmag = sqrt(u.*u + v.*v);

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

BPF = false(params.simRes, 'gpuArray');

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
% matlab coordinates
% get r and rVecs
%create a grid of points representing pixel positions in the field plane
%gridPoints = linspace(-params.gridSize,params.gridSize,params.simRes);



[x,z] = meshgrid(gridPoints, gridPoints); % field slice in the x z plane
y = ones(params.simRes,params.simRes)*ceil(params.a);   %field plane y = 0

params.rVecs = zeros(params.simRes*params.simRes, 3);
params.rVecs(:,1) = x(:); params.rVecs(:,2) = y(:); params.rVecs(:,3) = z(:); %r value at each pixel position
params.psVecs = bsxfun(@minus, params.rVecs, params.ps);
normPMinPs = sqrt(sum(params.psVecs.^2,2));
params.psMask=reshape(normPMinPs,params.simRes, params.simRes);
params.gpu_pf = gpuArray(params.pf);
params.displaySubplots=0;
alpha1 = asin(params.NA_in); alpha2 = asin(params.NA_out);
params.P  = ones(params.orderEf+1,1,'gpuArray');
params.Pl_cosalpha1 = gpuLegendre(params.orderEf+1,cos(alpha1),params.P);
params.Pl_cosalpha2 = gpuLegendre(params.orderEf+1,cos(alpha2),params.P);

origRVecs(:,1) = x(:); origRVecs(:,2) = y(:); origRVecs(:,3) = z(:);
rVecs_ps = bsxfun(@minus, origRVecs,params.ps);
params.normPMinPs = bsxfun(@rdivide, rVecs_ps,  sqrt(sum(rVecs_ps.^2,2)));
params.r_ps=reshape(sqrt(sum(rVecs_ps.^2,2)),params.simRes, params.simRes); %r value at each pixel position with respect to the center of the sphere
params.gpu_r_ps = gpuArray(params.r_ps);

params.gpu_Es = zeros(params.simRes,params.simRes,'gpuArray');
params.gpu_E_i = params.gpu_Es;
E_t = zeros(params.simRes,params.simRes,'gpuArray');

params.P  = ones(params.simRes,params.simRes,'gpuArray');
D_Et = zeros(params.res, params.res,'gpuArray');
D_Ef = D_Et;
allEtd = zeros(params.res, params.res,numWav,'gpuArray');
allEfd = allEtd;

A = gpuArray(zeros(params.res,params.res,numWav));
absSpec = gpuArray(zeros(numWav,1));


params

%%
display('gpu forward BimSim for all wavelengths ....')

h = waitbar(0, 'Per wavelength computation...');



for i=1:numWav
    
    if i==1
        tic
    end
    params.n = material(i,2) + 1i*material(i,3);
    
    params.lambda = material(i,1);
    
    params.wavNum = 2*pi/params.lambda;        %wavenumber
    
    params.numOrd = computeN_l(params.a, params.lambda);
    
    %create a vector of orders [0 1 ... Nl]
    ordVec = gpuArray((0:params.numOrd)');
    
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
    jl_ka = gpuArray(sphbesselj(params.numOrd,ka,'multiple'));
    %evaluate the derivate of the spherical bessel functions of the first kind at kna
    jl_kna_p = gpuArray(derivSphBes(params.numOrd, kna));
    %evaluate the spherical bessel functions of the first kind at kna
    jl_kna = gpuArray(sphbesselj(params.numOrd,kna,'multiple'));
    %evaluate the derivative of the spherical bessel functions of the first kind at ka
    jl_ka_p = gpuArray(derivSphBes(params.numOrd, ka));
    
    %compute the numerator for B coefficients
    numB = jl_ka.*jl_kna_p.*params.n - jl_kna.*jl_ka_p;
    %evaluate the derivative of the hankel functions of the first kind at ka
    hl_ka_p = gpuArray(derivSphHan(params.numOrd, ka));
    %evaluate the hankel functions of the first kind at ka
    hl_ka = gpuArray(shank1(params.numOrd, ka, 'multiple'));
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
    
    %     params.objMinParam = params.objectiveMin/params.lambda;
    %     params.objMaxParam = params.objectiveMax/params.lambda;
    %  BPF = gpuArray(bandpass(fmag, params.NA_in, params.NA_out, params.lambda));
    
    
    %params.BPF = fftshift((params.fmag < params.objMinParam | params.fmag > params.objMaxParam));
    
    params.kVec=(lightDirection*params.wavNum);
    params.normKvec = params.kVec./params.wavNum;
    
    params.k_j = gpuArray(monteCarlo(s,params.samples, params.kVec, params.NA_in, params.NA_out));
    
    %params.normKvec = params.kVec./params.wavNum;
    
    
    for p = 1:params.numPS
        
        params.rVecs = bsxfun(@minus, origRVecs,params.pf);
        %norm of the position vectors with respect to the focal point
        params.normPMinPf = bsxfun(@rdivide, params.rVecs,  sqrt(sum(params.rVecs.^2,2)));
        
        params.r=reshape(sqrt(sum(params.rVecs.^2,2)),params.simRes, params.simRes); %r value at each pixel position
        %compute and display Ef
        %display('Computing E_f...')
        params.P  = arrayfun(@(x) 1, params.P);
        
        %tic 
        E_f = gpuComputeEf(params);
        
        %toc 
        %display('Done computing E_f.')
        %compute and display E_s and E_i and E_t
        %display('Computing E_s, E_i, E_t...')
        %arrayfun(@(x) 0, params.gpu_Es);
        %arrayfun(@(x) 0, params.gpu_E_i);
        params.P = arrayfun(@(x) 1, params.P);
        %tic 
        [E_s, E_i] = gpuComputeScatteredFields(params);
        %toc 
        %display('Done computing E_s E_i.')
        E_t(params.psMask<params.a) = E_i(params.psMask<params.a);
        E_t(params.psMask>=params.a) = E_f(params.psMask>=params.a) + E_s(params.psMask>=params.a);
        
        
        E_t(isnan(E_t))=0; E_t(isinf(E_t))=0;
        % this is not correct if the slice goes through the sphere
        %Et_d = fftshift(fft2(E_t));
        Et_d = (fft2(E_t));
        
        
        fftEf = fft2(E_f);
        %Et_d(params.BPF) = 0;
        Et_d(BPF) = 0;
        fftEf(BPF) = 0;
        %iftEt_d = ifft2(ifftshift(Et_d));
        
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
        
        
%         params.pf = randn([1 3])*params.fov/5;
%         params.pf(2)=params.a;
        %     params.pf
        E_t = arrayfun(@(x) 0, E_t);
    end
    
    
    waitbar(i/numWav, h, sprintf('wavelength %f',material(i)));
    %calculate absorbance
    A(:,:,i) = -log10((D_Et)./(D_Ef));
    absSpec(i) = -log10(sum(D_Et(:))/sum(D_Ef(:)));
    allEtd(:,:,i) = D_Et;
    allEfd(:,:,i) = D_Ef;
    
    D_Et = arrayfun(@(x) 0, D_Et);
    D_Ef = arrayfun(@(x) 0, D_Ef);
    
    if i==1
        toc
    end
end
close(h)

%%
figure;ax=axes;
plot(wavenumbers(1:i),absSpec(1:i))
set(ax, 'Xdir','reverse')

%%
%save output
goutput.E_f = E_f;
goutput.E_s = E_s;
goutput.E_i = E_i;
goutput.E_t = E_t;
goutput.Et_d = D_Et;
goutput.Ef_d = D_Ef;
goutput.A = A;
goutput.absSpec = absSpec;
goutput.allEtd = allEtd;
goutput.allEfd = allEfd;

%show results in subplots
%showBimSim(1, E_f, E_s, E_i, E_t, D_Et, D_Ef, A)
% brewer = brewermap(1000);
% figure,
% imagesc((A)), title(sprintf('A after p = %i',p)),axis image, colormap(brewer), colorbar
display('Done.')
