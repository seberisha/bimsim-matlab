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

params.numPS=1;
params.a = 6.5;                     %radius of the sphere
params.ps = [0 0 0];
params.pf=([0 0 0]);
params.res = 129;
params.samples=1000;
params.orderEf=100;
params.objectiveMin = .2;   %inner objective NA
params.objectiveMax=.62;    %outer objective NA
params.NA_in = .2;
params.NA_out = 0.62;
params.E0 = 1;


params.fov = round(params.a)*4; %field of view in micrometers
%this gives the kVec direction from the y axis
theta=1.5708;
phi=0;
%direction of the incident light
[x,y,z] = sph2cart(theta,phi,1);
lightDirection  = [x y z];



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



%%%%%

[x_2d,z_2d] = meshgrid(gridPoints, gridPoints); % field slice in the x z plane
y_2d = ones(params.simRes,params.simRes)*(params.a);   %field plane y = 0

origRVecs_2d(:,1) = x_2d(:); origRVecs_2d(:,2) = y_2d(:); origRVecs_2d(:,3) = z_2d(:);
rVecs_ps_2d = bsxfun(@minus, origRVecs_2d,params.ps);
r_ps_2d=reshape(sqrt(sum(rVecs_ps_2d.^2,2)),params.simRes, params.simRes);

middleIdx = round(params.simRes/2)+1;
x = x_2d(middleIdx,middleIdx:end); y = y_2d(middleIdx,middleIdx:end); z = z_2d(middleIdx,middleIdx:end);

df = 1/(params.gridSize*2);

iu = 0:params.simRes-1;
%[iv, iu] = meshgrid(gridPoints, gridPoints);
u=zeros(size(iu));
idx = find(iu <= params.simRes/2);
u(idx) = iu(idx);
idx = find(iu > params.simRes/2);
u(idx) = (iu(idx) - params.simRes+1);
u=u.*df;
fmag = sqrt(u.*u);    %compute the magnitude of the frequencies
% 
% t = exp(-fmag.^2./(2*std(fmag)^2));


% df = 1/(params.gridSize*2);
% 
% [iv, iu] = meshgrid(0:params.simRes-1, 0:params.simRes-1);
% %[iv, iu] = meshgrid(gridPoints, gridPoints);
% u=zeros(size(iu)); v = u;
% idx = find(iu <= params.simRes/2);
% u(idx) = iu(idx);
% idx = find(iu > params.simRes/2);
% u(idx) = (iu(idx) - params.simRes+1);
% u=u.*df;
% idx = find(iv <= params.simRes/2);
% v(idx) = iv(idx);
% idx = find(iv > params.simRes/2);
% v(idx) = iv(idx) - params.simRes + 1;
% v=v.*df;
% 
% fmag = sqrt(u.*u + v.*v);    %compute the magnitude of the frequencies
% 
% fmag = fmag(middleIdx,:);

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
% y = ones(size(x))*ceil(params.a);   %field plane y = 0

params.rVecs = zeros(floor(params.simRes/2), 3);
params.rVecs(:,1) = x(:); params.rVecs(:,2) = y(:); params.rVecs(:,3) = z(:); %r value at each pixel position
params.psVecs = bsxfun(@minus, params.rVecs, params.ps);
normPMinPs = sqrt(sum(params.psVecs.^2,2));
params.psMask=normPMinPs;
alpha1 = asin(params.NA_in); alpha2 = asin(params.NA_out);
params.Pl_cosalpha1 = myLegendre(params.orderEf+1,cos(alpha1));
params.Pl_cosalpha2 = myLegendre(params.orderEf+1,cos(alpha2));


E_t = zeros(floor(params.simRes/2), 1);

Et_d = zeros(params.simRes, params.simRes);   %total field at the detector
Ef_d = Et_d;    %incident/focused field at the detector



A = (zeros(params.res,params.res,numWav));
absSpec = (zeros(numWav,1));

origRVecs(:,1) = x(:); origRVecs(:,2) = y(:); origRVecs(:,3) = z(:);
rVecs_ps = bsxfun(@minus, origRVecs,params.ps);
params.r_ps= sqrt(sum(rVecs_ps.^2,2)); %r value at each pixel position with respect to the center of the sphere


params.normPMinPs = bsxfun(@rdivide, rVecs_ps,  sqrt(sum(rVecs_ps.^2,2)));


allEtd = zeros(params.res, params.res,numWav);
allEfd = allEtd;

D_Et = zeros(floor(params.simRes/2), 1);
D_Ef = D_Et;

%interpolated fields at the detector
D_Eti = zeros(params.res,params.res);
D_Efi = zeros(params.res,params.res);

[x_2d,z_2d] = meshgrid(gridPoints, gridPoints); % field slice in the x z plane
y_2d = ones(params.simRes,params.simRes)*ceil(params.a);   %field plane y = 0

origRVecs_2d(:,1) = x_2d(:); origRVecs_2d(:,2) = y_2d(:); origRVecs_2d(:,3) = z_2d(:);
% rVecs_ps_2d = bsxfun(@minus, origRVecs_2d,params.ps);
% r_ps_2d=reshape(sqrt(sum(rVecs_ps_2d.^2,2)),params.simRes, params.simRes);

params

[H,k,r,I,K,R,h]=dht([],33,99,0);


%%
%tic
display('cpu forward BimSim for all wavelengths ....')

hx = waitbar(0, 'Per wavelength computation...');

for i=1:10
   mypf(i,:) = randn([1 3]); 
end



pOrder = 4; % change to 1 to get Figure 1(a)
R = params.gridSize;
Nr = floor(params.simRes/2);
h = hankel_matrix(pOrder, R, Nr, 1e-13);


% Samples of f1(r) are on the zeros of Bessel function
r = h.J_roots / (2*pi*h.vmax);
% Transformed vector's sample points
v = h.J_roots / (2*pi*h.rmax);

% The transforms
HT  = @(f1) (h.T  * (f1 ./ h.J * h.rmax)) .* h.J / h.vmax;
IHT = @(f2) (h.T' * (f2 ./ h.J * h.vmax)) .* h.J / h.rmax;



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
  %  nBPF = (t < objMinParam | t > objMaxParam);   %compute the band pass filter
    
    BPF = (fmag < objMinParam | fmag > objMaxParam);   %compute the band pass filter
     
    filter_1d = (v < objMinParam | v > objMaxParam);
    
    params.kVec=(lightDirection*params.wavNum);
    params.normKvec = params.kVec./params.wavNum;

    
    
    params.k_j = (monteCarlo(s,params.samples, params.kVec, params.NA_in, params.NA_out));
    params.pf=[0 0 0];
    
%     fs = 1/(2*params.gridSize);
%     idx = 1:params.simRes;
%     dist = ((idx - middleIdx).^2).^.5;
%     
 
    %dist = dist.*fs;
%     fc = linspace(-fs/2,fs/2,params.simRes);
    %fc = fc.*fs;
%     dist = abs(fc);
    %dist = abs(fftshift(u));
    %(1/(sqrt(2*pi)*objMaxParam))
%     filter1 = exp(-dist./(2*objMaxParam^2));
%     filter2 = exp(-dist./(2*objMinParam^2));
%     filter3 = 1 - filter2;
%     filter3 = filter1.*filter3;
    
%     g = mat2gray(psfGauss([1, params.simRes],objMinParam,fs));
%     hg = 1 - mat2gray(psfGauss([1, params.simRes],objMaxParam,fs));
%     bf = mat2gray(g.*hg);
    
    %filter3 = (gaussianbpf_1d(99,objMinParam,objMaxParam,(fmag)));
    
%     
%     tk = (gaussianBPF(params.simRes,objMinParam,objMaxParam,(k)));
%    
%     tt = (I*(tk./R).*K);
%     
%     tt = mat2gray(tt);
%     %filter3 = mat2gray(filter3);
%     t = changeRange(filter3, 0, 1);
    
    gbpf = gaussianBPF(v,objMinParam,objMaxParam);
    
    
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
        %         t = zeros(params.simRes,params.simRes);
        %         t(50,:)=E_t';
        %         tt = imrotate(t,0,'crop');
        %
        %         for ii=1:180
        %             er = imrotate(t,ii,'crop');
        % %             subplot(1,2,1), imagesc(abs(er)),axis image,colormap(brewer),colorbar
        %             tt(tt==0) = tt(tt==0) + er(tt==0);
        % %             subplot(1,2,2),imagesc(abs(tt)),axis image,colormap(brewer),colorbar
        % %             pause(1e-7)
        %         end
        %         Eti = tt;
        
        
        %%
%         fEt = fftshift(fft(E_t));
%         absFC = abs(fc);
%         kernelBP = absFC > objMaxParam & absFC < objMinParam;
%         bpfEt = kernelBP.*fEt;
%         figure, plot(fc, bpfEt)
%         ibpfEt = ifft(ifftshift(bpfEt));
%         figure, plot(gridPoints, real(ibpfEt))
        
        %% use hankel with gaussian filter
        
        Et_d = HT(E_t);
        %Et_d(filter_1d)=0;
        Et_d = Et_d.*gbpf;
        iftEt_d = IHT(Et_d);
        
        Ef_d = HT(E_f);
        Ef_d = Ef_d.*gbpf;
       % Ef_d(filter_1d)=0;
        iftEf_d = IHT(Ef_d);
        
%         Et_d = I*(E_t./R).*K;
%         Ef_d = I*(E_f./R).*K;
        
        %id = linspace(-params.simRes/2,params.simRes/2,params.simRes).*fs;
        %fmag = (fftshift(u));
        %fft of gaussian
%         lp = 1- mat2gray(exp(-2*pi^2*k.^2*objMinParam^2));
%         hp = 1-mat2gray(exp(-2*pi^2*k.^2*objMaxParam^2));
%         bp = 1-mat2gray(lp.*hp);
        
        %fft of circular disk
%         tl = (objMaxParam*besselj(1,pi*objMaxParam*k.*fs)./(k.*fs));
%         th = (1 - objMinParam*besselj(1,pi*objMinParam*(k.*fs))./(k.*fs));
%         bp = mat2gray(tl.*th);
        
%         fmag = fftshift(abs(u));
%         lp = exp( -(fmag) / (2*objMaxParam^2) ) / ( sqrt((2*pi)^2)*objMaxParam^2 );
%         temp = exp( -(fmag) / (2*objMinParam^2) ) / ( sqrt((2*pi)^2)*objMinParam^2 );
%         hp = 1-temp;
        
%         Et_d = Et_d.*gbpf;
%         Ef_d = Ef_d.*gbpf;
%         
%         iftEt_d=idht(Et_d,I,K,R);
%         iftEf_d=idht(Ef_d,I,K,R);
        
        %% use fft with gaussian filter
        %         Et_d = (fft(E_t));
        %         Ef_d = (fft(E_f));
        %
        %         Et_d = fftshift(Et_d).* filter3.';
        %         Ef_d = fftshift(Ef_d).* filter3.';
        %
        %         iftEt_d = ifft(ifftshift(Et_d));
        %         iftEf_d = ifft(ifftshift(Ef_d));
        
        %% use hankel with ideal filter
        
%         Et_d = I*(E_t./R).*K;
%         Ef_d = I*(E_f./R).*K;
%         
%         Et_d(BPF)=0;
%         Ef_d(BPF)=0;
%         
%         iftEt_d=idht(Et_d,I,K,R);
%         iftEf_d=idht(Ef_d,I,K,R);
        
        
        %% use fft with ideal filter
        
%         Et_d = (fft(E_t));
%         Ef_d = (fft(E_f));
%         
%         Et_d(BPF)=0;
%         Ef_d(BPF)=0;
%         
%         iftEt_d = ifft(ifftshift(Et_d));
%         iftEf_d = ifft(ifftshift(Ef_d));
        
        %%
        %integrate and resample
        D_Et = D_Et + (abs(iftEt_d)).^2;   %this is out_i in cuda bimsim -- the measured intesity
        D_Ef = D_Ef + (abs(iftEf_d)).^2;   %this out_inc in cuda bimsim -- the measured incident field
        
        params.pf = mypf(p,:);
      %  params.pf = randn([1 3]);
%         params.pf(2)=params.a;
        %     params.pf
        E_t = arrayfun(@(x) 0, E_t);
    end
    
    D_Eti = interp1(params.r_ps, D_Et, r_ps_2d);
    D_Efi = interp1(params.r_ps, D_Ef, r_ps_2d);
    
    
    %first crop the filtered near-field image of the source and scattered fields
    cropEt_d=D_Eti(params.startIdx:params.endIdx,params.startIdx:params.endIdx);
    cropEt_d(isnan(cropEt_d))=0;     cropEt_d(isinf(cropEt_d))=0;
    cropEf_d= D_Efi(params.startIdx:params.endIdx,params.startIdx:params.endIdx);
    cropEf_d(isnan(cropEf_d))=0;     cropEf_d(isinf(cropEf_d))=0;
    
    %calculate absorbance
    A(:,:,i) = -log10((cropEt_d)./(cropEf_d));
    absSpec(i) = -log10(sum(cropEt_d(:))/sum(cropEf_d(:)));
    allEtd(:,:,i) = cropEt_d;
    allEfd(:,:,i) = cropEf_d;
    
    D_Et = arrayfun(@(x) 0, D_Et);
    D_Ef = arrayfun(@(x) 0, D_Ef);
    
    if i==1
        toc
    end
    
    waitbar(i/numWav, hx, sprintf('wavelength %f',material(i)));
   
end
close(hx)

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
