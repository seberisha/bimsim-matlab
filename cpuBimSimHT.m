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
addpath(genpath('~/source/bimsim-matlab/HankelTransform/'));

%%
%load the file containing the wavelength and refractive index of the
%material
load pmma.mat
numWav = size(material,1);
wavenumbers = material(:,1);
wavenumbers = 1e4./wavenumbers;

params.numPS=1;
params.a = 6.5;                     %radius of the sphere
params.ps = [0 0 0];
params.pf=([0 0 0]);
params.res = 33;
params.samples=100;
params.orderEf=100;
params.objectiveMin = .2;   %inner objective NA
params.objectiveMax=.62;    %outer objective NA
params.NA_in = .2;
params.NA_out = 0.62;
params.E0 = 1;



%this gives the kVec direction from the y axis
theta=1.5708;
phi=0;
%direction of the incident light
[x,y,z] = sph2cart(theta,phi,1);
lightDirection  = [x y z];



params.fov = round(params.a)*3; %field of view in micrometers


%specify padding
padding = 1;    %use padding = 1, do not padd =0; padding is used for performing more accurate computations for
%specify the size of the field plane in wavelength units (microns)
params.gridSize = round(params.fov/2)*(2*padding + 1);

%specify the spatial resolution of the field plane
params.simRes = params.res*(2*padding + 1);

% compute alpha1 and alpha2 from NA_in and NA_out, respectively
params.alpha1 = asin(params.NA_in); params.alpha2 = asin(params.NA_out);


%compute the amplitude that makes it through the condenser
params.subA = 2 * pi * params.E0 * ( (1 - cos(params.alpha2)) - (1 - cos(params.alpha1)) );


%generate grid points
gridPoints = (2*params.gridSize)*(0:params.simRes-1)/params.simRes - params.gridSize;




[x_2d,z_2d] = meshgrid(gridPoints, gridPoints); % field slice in the x z plane
y_2d = ones(params.simRes,params.simRes)*(params.a);   %field plane y = 0

origRVecs_2d(:,1) = x_2d(:); origRVecs_2d(:,2) = y_2d(:); origRVecs_2d(:,3) = z_2d(:);
rVecs_ps_2d = bsxfun(@minus, origRVecs_2d,params.ps);
r_ps_2d=reshape(sqrt(sum(rVecs_ps_2d.^2,2)),params.simRes, params.simRes);

middleIdx = round(params.simRes/2);
x = x_2d(middleIdx,:); y = y_2d(middleIdx,:); z = z_2d(middleIdx,:);





df = 1/(params.gridSize*2);
iu = 0:params.simRes-1;

u=zeros(size(iu));
idx = find(iu <= params.simRes/2);
u(idx) = iu(idx);
idx = find(iu > params.simRes/2);
u(idx) = (iu(idx) - params.simRes+1);
u=u.*df;

fmag = sqrt(u.*u);    %compute the magnitude of the frequencies

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


%evaluate the field at the line where z=0 (middle of the plane)
% x = gridPoints;
% z = zeros(size(x));
% y = ones(size(x))*(params.a);   %field plane y = 0

params.rVecs = zeros(params.simRes, 3);
params.rVecs(:,1) = x(:); params.rVecs(:,2) = y(:); params.rVecs(:,3) = z(:); %r value at each pixel position
params.psVecs = bsxfun(@minus, params.rVecs, params.ps);
normPMinPs = sqrt(sum(params.psVecs.^2,2));
params.psMask=normPMinPs;
alpha1 = asin(params.NA_in); alpha2 = asin(params.NA_out);
params.Pl_cosalpha1 = myLegendre(params.orderEf+1,cos(alpha1));
params.Pl_cosalpha2 = myLegendre(params.orderEf+1,cos(alpha2));


E_t = zeros(params.simRes, 1);

Et_d = zeros(params.res, 1);   %total field at the detector
Ef_d = Et_d;    %incident/focused field at the detector



A = (zeros(params.res,params.res,numWav));
absSpec = (zeros(numWav,1));

origRVecs(:,1) = x(:); origRVecs(:,2) = y(:); origRVecs(:,3) = z(:);
rVecs_ps = bsxfun(@minus, origRVecs,params.ps);
params.r_ps= sqrt(sum(rVecs_ps.^2,2)); %r value at each pixel position with respect to the center of the sphere
params.normPMinPs = bsxfun(@rdivide, rVecs_ps,  sqrt(sum(rVecs_ps.^2,2)));


allEtd = zeros(params.res, params.res,numWav);
allEfd = allEtd;

D_Et = zeros(params.res, 1); %size required if we crop for each point source
%D_Et = zeros(params.simRes, 1); %size required if we crop after computing for all point
%sources
D_Ef = D_Et;

%interpolated fields at the detector
D_Eti = zeros(params.res,params.res);
D_Efi = zeros(params.res,params.res);

params

%%
%tic
display('cpu forward BimSim for all wavelengths ....')

h = waitbar(0, 'Per wavelength computation...');
for i=1:10
    mypf(i,:) = randn([1 3]);
end

for i=1:numWav
    
    if i==1
        tic
    end
    
    %tic
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
    %calculate the numerator for the scattering coefficients A
    numA = jl_ka.*hl_ka_p - jl_ka_p.*hl_ka;
    %calculate the scattering coefficients A
    params.A = params.twolp1_il.*(numA./denAB);
    
    
    objMinParam = params.NA_in/params.lambda; %min cut off of frequencies
    objMaxParam = params.NA_out/params.lambda; %max cut off of frequencies
    BPF = (fmag < objMinParam | fmag > objMaxParam);   %compute the band pass filter
    
    params.kVec=(lightDirection*params.wavNum);
    params.normKvec = params.kVec./params.wavNum;
    
    
    
    params.k_j = (monteCarlo(params.samples, params.kVec, params.NA_in, params.NA_out));
    params.pf=[0 0 0];
    for p = 1:params.numPS
        
        params.rVecs = bsxfun(@minus, origRVecs,params.pf);
        %norm of the position vectors with respect to the focal point
        normPMinPf = sqrt(sum(params.rVecs.^2,2));
        
        params.r=normPMinPf; %r value at each pixel position
        
        
        %compute and display Ef
        % display('Compting E_f...')
        E_f = computeEfSym(params);
        % display('Done computing E_f.')
        
        %compute and display E_s and E_i and E_t
        %  display('Computing E_s, E_i, E_t...')
        [E_s, E_i] = computeEsEiSym(params);
        
        
        
        E_t(params.psMask<params.a) = E_i(params.psMask<params.a);
        E_t(params.psMask>=params.a) = E_f(params.psMask>=params.a) + E_s(params.psMask>=params.a);
        
        E_t(isnan(E_t)) = 0;    E_t(isinf(E_t))=0;%near field = out_n in cuda bimsim
        % this is not correct if the slice goes through the sphere
        
        %Et_d = (fft(E_t));
        Et_d = ht(E_t).';
        
        %Et_d((BPF)) = 0;
        
        %iftEt_d = ifft((Et_d));
        iftEt_d = iht(Et_d).';
        
        % crop at each point source
        %first crop the filtered near-field image of the source and scattered fields
        cropEt_d=iftEt_d(params.startIdx:params.endIdx);
        cropEt_d(isnan(cropEt_d))=0;     cropEt_d(isinf(cropEt_d))=0;
        cropEf_d= E_f(params.startIdx:params.endIdx);
        cropEf_d(isnan(cropEf_d))=0;     cropEf_d(isinf(cropEf_d))=0;
        
        %integrate and resample
        D_Et = D_Et + (abs(cropEt_d)).^2;   %this is out_i in cuda bimsim -- the measured intesity
        D_Ef = D_Ef + (abs(cropEf_d)).^2;   %this out_inc in cuda bimsim -- the measured incident field
        
        
        %         %no cropping
        %         %integrate and resample
        %         D_Et = D_Et + (abs(iftEt_d)).^2;   %this is out_i in cuda bimsim -- the measured intesity
        %         D_Ef = D_Ef + (abs(E_f)).^2;   %this out_inc in cuda bimsim -- the measured incident field
        %
        
        % params.pf = mypf(p,:);
        %         params.pf(2)=params.a;
        %     params.pf
        E_t = arrayfun(@(x) 0, E_t);
    end
    
    %interpolate the cropped intensities
    cropped_r_ps = params.r_ps(params.startIdx:params.endIdx);
    cropped_r_ps_2d = r_ps_2d(params.startIdx:params.endIdx,params.startIdx:params.endIdx);
    zero_idx = round(params.res/2)+1;
    
    D_Eti = interp1(cropped_r_ps(zero_idx:end), D_Et(zero_idx:end), cropped_r_ps_2d);
    D_Efi = interp1(cropped_r_ps(zero_idx:end), D_Ef(zero_idx:end), cropped_r_ps_2d);
    
    D_Eti(isnan(D_Eti)) =0; D_Eti(isnan(D_Eti)) =0;
    D_Efi(isnan(D_Efi)) =0; D_Efi(isnan(D_Efi)) =0;
    
    
    %interpolate full intensities and then crop
    %     zero_idx = round(params.simRes/2)+1;
    %     D_Eti = interp1(params.r_ps(zero_idx:end), D_Et(zero_idx:end), r_ps_2d);
    %     D_Efi = interp1(params.r_ps(zero_idx:end), D_Ef(zero_idx:end), r_ps_2d);
    %
    %     cropEt_d=D_Eti(params.startIdx:params.endIdx,params.startIdx:params.endIdx);
    %     cropEt_d(isnan(cropEt_d))=0;     cropEt_d(isinf(cropEt_d))=0;
    %     cropEf_d= D_Efi(params.startIdx:params.endIdx,params.startIdx:params.endIdx);
    %     cropEf_d(isnan(cropEf_d))=0;     cropEf_d(isinf(cropEf_d))=0;
    
    
    %     cropEt_d(isnan(cropEt_d)) =0; cropEt_d(isnan(cropEt_d)) =0;
    %     cropEf_d(isnan(cropEf_d)) =0; cropEf_d(isnan(cropEf_d)) =0;
    
    
    %calculate absorbance
    % version of cropping for each point source
    A(:,:,i) = -log10((D_Eti)./(D_Efi));
    absSpec(i) = -log10(sum(D_Eti(:))/sum(D_Efi(:)));
    allEtd(:,:,i) = D_Eti;
    allEfd(:,:,i) = D_Efi;
    
    %version of cropping after computing for all point sources
    %     A(:,:,i) = -log10((cropEt_d)./(cropEf_d));
    %     absSpec(i) = -log10(sum(cropEt_d(:))/sum(cropEf_d(:)));
    %     allEtd(:,:,i) = cropEt_d;
    %     allEfd(:,:,i) = cropEf_d;
    
    
    D_Et = arrayfun(@(x) 0, D_Et);
    D_Ef = arrayfun(@(x) 0, D_Ef);
    
    if i==1
        toc
    end
    
    waitbar(i/numWav, h, sprintf('wavelength %f',material(i)));
    
end
close(h)

%save output
soutput.E_f = E_f;
soutput.E_s = E_s;
soutput.E_i = E_i;
soutput.E_t = E_t;
soutput.Et_d = D_Et;
soutput.Ef_d = D_Ef;
soutput.A = A;
soutput.allEtd = allEtd;
soutput.allEfd = allEfd;

%%
%%
figure;ax=axes;

plot(wavenumbers(1:i),absSpec(1:i))
set(ax, 'Xdir','reverse')
%%

%show results in subplots
%showBimSim(1, E_f, E_s, E_i, E_t, Et_d, Ef_d, A)

display('Done.')
