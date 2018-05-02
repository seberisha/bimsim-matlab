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

%%
load pmma.mat
numWav = size(material,1);
%numWav=1;

wavenumbers = material(:,1);
wavenumbers = 1e4./wavenumbers;

%%

params.a = 6.5;                     %radius of the sphere
params.ps = [0 0 0];
params.pf=([0 0 0]);
params.res = 33;

%sphere

params.fov = ceil(params.a)*2; %field of view in micrometers


numPS=1;
% load fpv100_6.5.mat
% numPS = size(fpv,1);

% load 529PointSources_r32.mat
% numPS = size(focalPoints,1);
%numPS=1;
%create a parameter structure for the simulation

params.samples=100;    %number of samples used for the Monte Carlo simulation of the scattered and internal fields

params.orderEf=100; %number of orders used for computing the incident/focused field

params.objectiveMin = .2;   %inner objective NA
params.objectiveMax=.62;    %outer objective NA


theta=1.5708; phi=0; %direction of the incident light: theta (0 - 2pi), phi (0-pi), r=1



%params.numPS=1; %number of point sources
padding = 1;    %use padding = 1, do not padd =0; padding is used for performing more accurate computations for
%measuring the field at the detector.
params.gridSize = round(params.fov/2)*(2*padding + 1);   %specify the size of the field plane in wavelength units (micrometers)

params.simRes = params.res*(2*padding + 1); %compute the spatial resolution of the field plane
params.E0 = 1;  %amplitude of the field
params.NA_in = .2; %inner and outer NA of the condenser
params.NA_out = 0.62;
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

[x,y,z] = sph2cart(theta,phi,1);    %direction of the incident light in cartesian coordinates
lightDirection = [x y z];
%%%%
%create a grid of points representing pixel positions in the field plane
gridPoints = linspace(-params.gridSize,params.gridSize,params.simRes); %in case of symmetry this goes from 0 to gridSize, instead of -gridSize to gridSize
%%%%%

%generate grid points
gridPoints = (2*params.gridSize)*(0:params.simRes-1)/params.simRes - params.gridSize;


%%%%
x = gridPoints'; % field slice in the x plane -- in case of symmetry choose x to be between 0 and gridSize
y = ones(params.simRes,1)*ceil(params.a);   %field plane y = 0 --vector if symmetry
z = zeros(params.simRes,1);   %field plane z = 0 --vector if symmetry
params.rVecs = zeros(params.simRes, 3);
%%%%

params.rVecs(:,1) = x(:); params.rVecs(:,2) = y(:); params.rVecs(:,3) = z(:); %r value at each pixel position
params.psVecs = bsxfun(@minus, params.rVecs, params.ps);
params.normPMinPs = sqrt(sum(params.psVecs.^2,2));
params.normPsVecs = bsxfun(@rdivide, params.psVecs,  params.normPMinPs);

%%%%
params.psMask = params.normPMinPs;   %used to mask out the values inside the sphere area -- vector if symmetry
%%%%


params.displaySubplots=0;
alpha1 = asin(params.NA_in); alpha2 = asin(params.NA_out);
params.Pl_cosalpha1 = myLegendre(params.orderEf+1,cos(alpha1));
params.Pl_cosalpha2 = myLegendre(params.orderEf+1,cos(alpha2));

%create a grid of points representing pixel positions in the 2d field plane
%gridPoints = linspace(-params.gridSize, params.gridSize,128);

[x_2d,z_2d] = meshgrid(gridPoints, gridPoints); % 2d field slice in the x z plane
y_2d = ones(params.simRes,params.simRes)*ceil(params.a);   % 2d field plane y = 0

rVecs_2d = zeros(params.simRes*params.simRes, 3);
rVecs_2d(:,1) = x_2d(:); rVecs_2d(:,2) = y_2d(:); rVecs_2d(:,3) = z_2d(:); %r value at each pixel position
psVecs_2d = bsxfun(@minus, rVecs_2d, params.ps);
%normPMinPs = sqrt(sum(params.psVecs.^2,2));
normPMinPs_2d = sqrt(sum(psVecs_2d.^2,2));
r_2d=reshape(normPMinPs_2d,params.simRes, params.simRes);


%[~,~,r_2d] = cart2sph(x_2d,y_2d,z_2d);
% halfIdx = params.simRes/2 + 1;
% [~,~,r_1d] = cart2sph(gridPoints, zeros(size(gridPoints)), zeros(size(gridPoints)));



A = (zeros(params.res,params.res,numWav));
absSpec = (zeros(numWav,1));


%%%%
E_t = zeros(params.simRes, 1); %vector if symmetry
%%%%

%%%%
Et_d = zeros(params.simRes, 1);   %total field at the detector - vector in case of symmetry
%%%%
Ef_d = Et_d;

brewer = brewermap(1000);

% D = params.gridSize*2;  %used to compute delta d needed for computing the indices for the band pass filter
% % Set up range of variables.
% u = 0:(params.simRes-1);
% %v = 0:(params.simRes-1);
% % Compute the indices for use in meshgrid
% idx = find(u > params.simRes/2);
% u(idx) = u(idx) - params.simRes;
% %idy = find(v > params.simRes/2);
% %v(idy) = v(idy) - params.simRes;
% % Compute the meshgrid arrays
% %[v, u] = meshgrid(v, u);
% df = 1/D;
% u=u.*df;
% %v=v.*df;

df = 1/(params.gridSize*2);

iu = 0:params.simRes-1;
u=zeros(size(iu));
idx = find(iu <= params.simRes/2);
u(idx) = iu(idx);
idx = find(iu > params.simRes/2);
u(idx) = (iu(idx) - params.simRes+1);
u=u.*df;

fmag = sqrt(u.*u);    %compute the magnitude of the frequencies

h = waitbar(0, 'Per wavelength computation...');
%figure;ax=axes;
%figure, ax = axes;drawnow

allEtd = zeros(params.res,params.res,numWav);
allEfd = allEtd;

%%
for i=1:numWav
    params.pf=([0 0 0]);
    
    if i==1
        tic
    end
    params.n = material(i,2) + 1i*material(i,3);
    params.lambda = material(i,1);
    
    
    params.wavNum = 2*pi/params.lambda; %wavenumber
    params.numOrd = computeN_l(params.a, params.lambda);    %the maximum order required for convergence - this is used for computing the scattered
    %and internal
    %fields
    
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
    
    %     params.fmag = sqrt(u.*u);    %compute the magnitude of the frequencies
    %     params.objMinParam = params.objectiveMin/params.lambda; %min cut off of frequencies
    %     params.objMaxParam = params.objectiveMax/params.lambda; %max cut off of frequencies
    %     params.BPF = fftshift((params.fmag < params.objMinParam | params.fmag > params.objMaxParam));   %compute the band pass filter
    
    objMinParam = params.NA_in/params.lambda; %min cut off of frequencies
    objMaxParam = params.NA_out/params.lambda; %max cut off of frequencies
    BPF = (fmag < objMinParam | fmag > objMaxParam);   %compute the band pass filter
    
    params.kVec=lightDirection*params.wavNum;  %incident light vector scaled by the wavenumber
    
    %generate a 'params.samples' number of k vectors in the direction of
    %incident light using Monte Carlo simulation
    params.k_j = monteCarlo(params.samples, params.kVec, params.NA_in, params.NA_out);
    
    
    %%
    %tic
    % display('starting computations for point sources....')
    %numPS=100;
    params.pf= params.ps;
    for p = 1:numPS
        %params.pf = fpv(p,:);
        %params.pf
        %         params.pf = focalPoints(p,:);
        %         params.pf(3) = params.pf(2);
        %         params.pf(2) = params.a/2;
        %params.pf
        %         params.pf(2) = 0;
        
        params.rVecs(:,1) = x(:); params.rVecs(:,2) = y(:); params.rVecs(:,3) = z(:);
        params.rVecs = bsxfun(@minus, params.rVecs,params.pf);
        %norm of the position vectors with respect to the focal point
        normPMinPf = sqrt(sum(params.rVecs.^2,2));
        
        %%%%
        params.r= normPMinPf; %r value at each pixel position --vector if symmetry
        %%%%
        
        %compute and display Ef
        %display('Compting E_f...')
        E_f = computeEf_sym(params);
        
        %display('Done computing E_f.')
        
        
        
        %compute and display E_s and E_i and E_t
        % display('Computing E_s, E_i, E_t...')
        [E_s, E_i] = computeScatteredFields_sym(params);
        
        E_t(params.psMask<params.a) = E_i(params.psMask<params.a);
        E_t(params.psMask>=params.a) = E_f(params.psMask>=params.a) + E_s(params.psMask>=params.a);
        
        %replace NaNs with zeros
        %note: If the input to fft contains any NaN elements, then all the output elements will be NaN
        E_t(isnan(E_t))=0;
        
       % this is not correct if the slice goes through the sphere
        E_d = fft(E_t);
        
        %         N = length(E_d);
        %         k = ifftshift(-floor(N/2):floor((N-1)/2))'; %//compute the frequency vector
        %         E_d = E_d.*(abs(k).*df<= objMaxParam); %// zero out all frequencies larger than 'm'
        %         E_d = E_d.*(abs(k).*df>= objMinParam);
        
        
        
        E_d(BPF) = 0;
        iftE_d = ifft((E_d));
        
        %integrate and resample
        Et_d = Et_d + (abs(iftE_d)).^2;
        Ef_d = Ef_d + (abs(E_f)).^2;
        %calculate absorbance
        %figure, plot(A)
        
        %Ai = interp1(r_1d(65:end)', A(65:end), r_2d);
        %figure, imagesc(Ai),axis image, colormap(brewer)
        
        
        
        %params.pf
       % params.pf = randn([1 3])*params.fov/5;
        %         params.pf(2)=0;
        %params.pf
        E_t = arrayfun(@(x) 0, E_t);
    end
    %display('time per wavelength:')
    %toc
    
    Et_d(isnan(Et_d))=0;     Et_d(isinf(Et_d))=0;
    Ef_d(isnan(Ef_d))=0;     Ef_d(isinf(Ef_d))=0;
    
    
    Etd_interp = interp1(params.normPMinPs(params.simRes/2+1:end)', Et_d(params.simRes/2+1:end), r_2d);
    Efd_interp = interp1(params.normPMinPs(params.simRes/2+1:end)', Ef_d(params.simRes/2+1:end), r_2d);
    
    %first crop the filtered near-field image of the source and scattered fields
    cropEd=Etd_interp(params.startIdx:params.endIdx,params.startIdx:params.endIdx);
    cropEd(isnan(cropEd))=0;     cropEd(isinf(cropEd))=0;
    cropEf=Efd_interp(params.startIdx:params.endIdx,params.startIdx:params.endIdx);
    cropEf(isnan(cropEf))=0;     cropEf(isinf(cropEf))=0;
    
    %calculate absorbance
    A(:,:,i) = -log10((cropEd)./(cropEf));
    % subplot(1,2,1), imagesc(A(:,:,i)),colormap(brewer),colorbar, axis image, title(sprintf('lambda = %f',params.lambda))
    absSpec(i) = -log10(sum(cropEd(:))/sum(cropEf(:)));
    
    allEtd(:,:,i) = cropEd;
    allEfd(:,:,i) = cropEf;
    
    
    %   subplot(1,2,2),
%     plot(wavenumbers(1:i),absSpec(1:i))
%     set(ax, 'Xdir','reverse')
%     drawnow
    
    Et_d = arrayfun(@(x) 0, Et_d);
    Ef_d = arrayfun(@(x) 0, Ef_d);
    if i==1
        toc
    end
    waitbar(i/numWav, h, sprintf('wavelength %f',material(i)));
    %i
    
end

close(h)

%%
figure;ax=axes;
plot(wavenumbers(1:i),absSpec(1:i))
set(ax, 'Xdir','reverse')

%%
%save output
soutput.E_f = E_f;
soutput.E_s = E_s;
soutput.E_i = E_i;
soutput.E_t = E_t;
soutput.Et_d = Et_d;
soutput.Ef_d = Ef_d;
soutput.A = A;
soutput.absSpec = absSpec;
soutput.allEtd = allEtd;
soutput.allEfd = allEfd;
%save /home/sberisha/source/bimsim_matlab/data/a6.5_fo100_r32_ps529/mc10/cpuSymResults

display('Done.')
