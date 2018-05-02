%% setup
%clear
addpath(genpath('~/source/stim-matlab/'))

n=1.5 + 1i*0.1;
fov = 5;
%wavenumber in cm^-1
waveNumber = 3000;
%wavelength in micron
lambda = 1e4/waveNumber;
samples=200;
numPS=30 ;
a = 2.5;                     %radius of the sphere

%specify padding
padding = 1;
%specify the size of the field plane in wavelength units (microns)
gridSize = fov*(2*padding + 1);
res = 128;
%specify the spatial resolution of the field plane
simRes = res*(2*padding + 1);
%create a parameter structure for the simulation
ps = [0 0 0];
orderEf=200;
E0 = 1;
NA_in = .62;
NA_out = 0.2;
wavNum = 2*pi/lambda;        %wavenumber

numOrd = computeN_l(a, lambda);


% compute alpha1 and alpha2 from NA_in and NA_out, respectively
alpha1 = asin(NA_in); alpha2 = asin(NA_out);


%create a vector of orders [0 1 ... Nl]
ordVec = (0:numOrd)';

%The scattered field for a single incident plane-wave k produced by
%a sphere with radius r positioned at point pf

%calculate the prefix term (2l + 1)*i^l
twolp1 = 2.*ordVec+1;
il = 1i.^ordVec;
twolp1_il = twolp1.*il;

%compute the arguments needed to evaluate spherical bessel functions,
%hankel functions, and their derivatives
ka=wavNum*a;

kna = wavNum*n*a;

%evaluate the spherical bessel functions of the first kind at ka
jl_ka = sphbesselj(numOrd,ka,'multiple');
%evaluate the derivate of the spherical bessel functions of the first kind at kna
jl_kna_p = derivSphBes(numOrd, kna);
%evaluate the spherical bessel functions of the first kind at kna
jl_kna = sphbesselj(numOrd,kna,'multiple');
%evaluate the derivative of the spherical bessel functions of the first kind at ka
jl_ka_p = derivSphBes(numOrd, ka);

%compute the numerator for B coefficients
numB = jl_ka.*jl_kna_p.*n - jl_kna.*jl_ka_p;
%evaluate the derivative of the hankel functions of the first kind at ka
hl_ka_p = derivSphHan(numOrd, ka);
%evaluate the hankel functions of the first kind at ka
hl_ka = shank1(numOrd, ka, 'multiple');
%compute the denominator for coefficients A and B
denAB = jl_kna.*hl_ka_p - hl_ka.*jl_kna_p*n;

%compute B
B = twolp1_il.*(numB./denAB);

%calculate the numerator for the scattering coefficients A
numA = jl_ka.*hl_ka_p - jl_ka_p.*hl_ka;
%calculate the scattering coefficients A
A = twolp1_il.*(numA./denAB);

%compute the amplitude that makes it through the condenser
subA = 2 * pi * E0 * ( (1 - cos(alpha2)) - (1 - cos(alpha1)) );



D = gridSize*2;
% Set up range of variables.
u = 0:(simRes-1);
v = 0:(simRes-1);
% Compute the indices for use in meshgrid
idx = find(u > simRes/2);
u(idx) = u(idx) - simRes;
idy = find(v > simRes/2);
v(idy) = v(idy) - simRes;
% Compute the meshgrid arrays
[v, u] = meshgrid(v, u);
df = 1/D;
u=u.*df;
v=v.*df;
objectiveMin = .2;
objectiveMax=.62;
fmag = sqrt(u.*u + v.*v);
objMinParam = objectiveMin/lambda;
objMaxParam = objectiveMax/lambda;
BPF = fftshift((fmag < objMinParam | fmag > objMaxParam));
cropSize = padding*res;

if padding==0
    startIdx=1;
    endIdx=simRes;
else
    startIdx = round((simRes  - cropSize)/2);
    endIdx = startIdx + cropSize-1;
end

D_Ed = zeros(res, res);
D_Ef = D_Ed;
theta=1.5708; phi=0; %this gives the kVec direction from the y axis
%direction of the incident light
[x,y,z] = sph2cart(theta,phi,1);
kVec=[x y z]*wavNum;
% matlab coordinates
% get r and rVecs
%create a grid of points representing pixel positions in the field plane
gridPoints = linspace(-gridSize,gridSize,simRes);
k_j = monteCarlo(samples, kVec, NA_in, NA_out);
rVecs = zeros(simRes*simRes, 3);
[x,z] = meshgrid(gridPoints, gridPoints); % field slice in the x z plane
y = zeros(simRes,simRes);   %field plane y = 0
rVecs(:,1) = x(:); rVecs(:,2) = y(:); rVecs(:,3) = z(:);

psVecs = bsxfun(@minus, rVecs, ps);
normPMinPs = sqrt(sum(psVecs.^2,2));
psMask=reshape(normPMinPs,simRes, simRes);
pf=[0 0 0];
displaySubplots=0;
alpha1 = asin(NA_in); alpha2 = asin(NA_out);
Pl_cosalpha1 = squeeze(myLegendre(orderEf+1,cos(alpha1)));
Pl_cosalpha2 = squeeze(myLegendre(orderEf+1,cos(alpha2)));




%%
POOL = parpool('local',6);
tic
display('simulating the full field....')
absCroppedEds = zeros(res,res,numPS-1);
absCroppedEfs = zeros(res,res,numPS-1);


rVecs = bsxfun(@minus, rVecs,pf);
%norm of the position vectors with respect to the focal point
normPMinPf = sqrt(sum(rVecs.^2,2));
r=reshape(normPMinPf,simRes, simRes); %r value at each pixel position
%compute and display Ef
%E_f = parComputeEf(rVecs,kVec,wavNum,simRes,orderEf,r,Pl_cosalpha1,Pl_cosalpha2,E0);



%% focused field
%apply the partial wave expansion to get the focused field

normrVecs = bsxfun(@rdivide, rVecs,  sqrt(sum(rVecs.^2,2)));
normKvec = kVec./wavNum;
cos_theta = (normKvec*normrVecs')';
cos_theta = reshape(cos_theta, simRes, simRes);

%calculate the prefix term (2l + 1)*i^l
ordVecEf=(0:orderEf)';
il = 1i.^ordVecEf;

jl_kr = sphbesselj(orderEf,wavNum.*r,'multiple');
Pl_costheta = myLegendre(orderEf,cos_theta);
jlkr_Pcostheta = jl_kr.*Pl_costheta;

il_jlkr_Pcostheta = (bsxfun(@times, jlkr_Pcostheta, reshape(il,[1 1 orderEf+1])));

ord=0;
il_jlkr_Pcostheta(:,:,ord+1) = il_jlkr_Pcostheta(:,:,ord+1).*(Pl_cosalpha1(ord+2)-Pl_cosalpha2(ord+2)-Pl_cosalpha1(1)+Pl_cosalpha2(1));

ord=1;
il_jlkr_Pcostheta(:,:,ord+1) = il_jlkr_Pcostheta(:,:,ord+1).*(Pl_cosalpha1(ord+2)-Pl_cosalpha2(ord+2)-Pl_cosalpha1(1)+Pl_cosalpha2(1));

slicePl_cosalpha1 = zeros(size(Pl_cosalpha1));
slicePl_cosalpha1(4:end) = Pl_cosalpha1(4:end);

slicePl_cosalpha2 = zeros(size(Pl_cosalpha2));
slicePl_cosalpha2(4:end) = Pl_cosalpha2(4:end);

parfor ord=2:orderEf
    il_jlkr_Pcostheta(:,:,ord+1) = il_jlkr_Pcostheta(:,:,ord+1).*(slicePl_cosalpha1(ord+2)-slicePl_cosalpha2(ord+2)-Pl_cosalpha1(ord)+Pl_cosalpha2(ord));
end

E_f = 2*pi*E0.*sum(il_jlkr_Pcostheta,3);

%%

rVecs(:,1) = x(:); rVecs(:,2) = y(:); rVecs(:,3) = z(:);
rVecs = bsxfun(@minus, rVecs,ps);
normPMinPs = sqrt(sum(rVecs.^2,2));
r=reshape(normPMinPs,simRes, simRes); %r value at each pixel position with respect to the center of the sphere
%compute and display E_s and E_i and E_t
%[E_s, E_i] = parComputeScatteredFields(a,psMask, n,rVecs,wavNum,r,numOrd,simRes,ps,pf,samples,k_j,subA,B,A);




%normalize the position vectors
normRvecs = bsxfun(@rdivide, rVecs,  sqrt(sum(rVecs.^2,2)));

kr = wavNum.*r;
hl_kr = shank1(numOrd, kr, 'multiple');

%The internal field specifying the field inside of a sphere
%for an incident plane wave

%compute the argument for the spherical bessel function k*n*r
knr = wavNum * n .* r;

%compute the spherical bessel function
jl_knr = sphbesselj(numOrd,knr,'multiple');

E_s = zeros(simRes,simRes);
E_i = E_s;
%h = waitbar(0, 'Monte Carlo integration...');


% r = sqrt(normKvec(1)^2 + normKvec(2)^2 + normKvec(3)^2);
% theta =  atan2(normKvec(2), normKvec(1));
% phi = acos(normKvec(2)/r);
% kSpherical = [r theta phi]

%the vector from the focal point to the center of the sphere
c = ps - pf;
display('Starting Monte Carlo simulation...')


parfor i=1:samples
    cos_theta = (k_j(:,i)'*normRvecs')';
    cos_theta = reshape(cos_theta, simRes, simRes);
    %compute the legendre polynomials needed for computing E_s and E_i
    Pl_costheta = myLegendre(numOrd,cos_theta);
    hlkr_Plcostheta = hl_kr.*Pl_costheta;
    
    %multiply by the legendre polynomial
    jlknr_Plcostheta = jl_knr.*Pl_costheta;
    
    %multiply by the scattering coefficients B
    %sum all orders
    phase = exp(1i.*wavNum.*k_j(:,i)'*c');
    E_s = E_s + (1/samples).*subA.*phase.*sum(bsxfun(@times, hlkr_Plcostheta, reshape(B,[1 1 numOrd+1])),3);
    %E_s(r<A) = 0;%E_i(r<a);
    %subplot(1,2,1), imagesc((abs((E_s)))),title('E_s'), colorbar, axis image, colormap(brewer)
    E_i = E_i + (1/samples).*subA.*phase.*sum(bsxfun(@times, jlknr_Plcostheta, reshape(A,[1 1 numOrd+1])),3);
    %E_i(r>A) = 0;
    %subplot(1,2,2), imagesc((abs((E_i)))),title('E_i'), colorbar, axis image, colormap(brewer)
    %pause(0.1)
    %waitbar(i/samples);
end
display('Finished Monte Carlo simulation...')

%close(h)

E_s(psMask<a) = 0;%E_i(r<a);

E_i(psMask>a) = 0;






E_t = zeros(simRes, simRes);
E_t(psMask<a) = E_i(psMask<a);
E_t(psMask>=a) = E_f(psMask>=a) + E_s(psMask>=a);
% this is not correct if the slice goes through the sphere
E_d = fft2(E_t);
E_d = fftshift((E_d));
E_d(BPF) = 0;
iftE_d = ifft2(ifftshift(E_d));
%first crop the filtered near-field image of the source and scattered fields
cropEd=iftE_d(startIdx:endIdx,startIdx:endIdx);
cropEd(isnan(cropEd))=0;     cropEd(isinf(cropEd))=0;
cropEf=E_f(startIdx:endIdx, startIdx:endIdx);
cropEf(isinf(cropEf))=0;     cropEf(isinf(cropEf))=0;
%save the absolute values of the cropped detected total field and the
%detected focused field
D_Ed = (abs(cropEd)).^2;
D_Ef = (abs(cropEf)).^2;

%%




for p = 2:numPS
    pf = randn([1 3]);
    pf(2)=0;
    rVecs = zeros(simRes*simRes, 3);
    [x,z] = meshgrid(gridPoints, gridPoints); % field slice in the x z plane
    y = zeros(simRes,simRes);   %field plane y = 0
    rVecs(:,1) = x(:); rVecs(:,2) = y(:); rVecs(:,3) = z(:);
    rVecs = bsxfun(@minus, rVecs,pf);
    %norm of the position vectors with respect to the focal point
    normPMinPf = sqrt(sum(rVecs.^2,2));
    r=reshape(normPMinPf,simRes, simRes); %r value at each pixel position
    %compute and display Ef
    %display('Compting E_f...')
    % E_f = parComputeEf(rVecs,kVec,wavNum,simRes,orderEf,r,Pl_cosalpha1,Pl_cosalpha2,E0);
    
    %% focused field
    %apply the partial wave expansion to get the focused field
    
    normrVecs = bsxfun(@rdivide, rVecs,  sqrt(sum(rVecs.^2,2)));
    normKvec = kVec./wavNum;
    cos_theta = (normKvec*normrVecs')';
    cos_theta = reshape(cos_theta, simRes, simRes);
    
    %calculate the prefix term (2l + 1)*i^l
    ordVecEf=(0:orderEf)';
    il = 1i.^ordVecEf;
    
    jl_kr = sphbesselj(orderEf,wavNum.*r,'multiple');
    Pl_costheta = myLegendre(orderEf,cos_theta);
    jlkr_Pcostheta = jl_kr.*Pl_costheta;
    
    il_jlkr_Pcostheta = (bsxfun(@times, jlkr_Pcostheta, reshape(il,[1 1 orderEf+1])));
    
    ord=0;
    il_jlkr_Pcostheta(:,:,ord+1) = il_jlkr_Pcostheta(:,:,ord+1).*(Pl_cosalpha1(ord+2)-Pl_cosalpha2(ord+2)-Pl_cosalpha1(1)+Pl_cosalpha2(1));
    
    ord=1;
    il_jlkr_Pcostheta(:,:,ord+1) = il_jlkr_Pcostheta(:,:,ord+1).*(Pl_cosalpha1(ord+2)-Pl_cosalpha2(ord+2)-Pl_cosalpha1(1)+Pl_cosalpha2(1));
    
    slicePl_cosalpha1 = zeros(size(Pl_cosalpha1));
    slicePl_cosalpha1(4:end) = Pl_cosalpha1(4:end);
    
    slicePl_cosalpha2 = zeros(size(Pl_cosalpha2));
    slicePl_cosalpha2(4:end) = Pl_cosalpha2(4:end);
    
    parfor ord=2:orderEf
        il_jlkr_Pcostheta(:,:,ord+1) = il_jlkr_Pcostheta(:,:,ord+1).*(slicePl_cosalpha1(ord+2)-slicePl_cosalpha2(ord+2)-Pl_cosalpha1(ord)+Pl_cosalpha2(ord));
    end
    
    E_f = 2*pi*E0.*sum(il_jlkr_Pcostheta,3);
    
    
    
    %display('Done computing E_f.')
    
    rVecs(:,1) = x(:); rVecs(:,2) = y(:); rVecs(:,3) = z(:);
    rVecs = bsxfun(@minus, rVecs,ps);
    normPMinPs = sqrt(sum(rVecs.^2,2));
    r=reshape(normPMinPs,simRes, simRes); %r value at each pixel position with respect to the center of the sphere
    %compute and display E_s and E_i and E_t
    display('Computing E_s, E_i, E_t...')
    
    %[E_s, E_i] = parComputeScatteredFields(a,psMask, n,rVecs,wavNum,r,numOrd,simRes,ps,pf,samples,k_j,subA,B,A);
    
    
    %normalize the position vectors
    normRvecs = bsxfun(@rdivide, rVecs,  sqrt(sum(rVecs.^2,2)));
    
    kr = wavNum.*r;
    hl_kr = shank1(numOrd, kr, 'multiple');
    
    %The internal field specifying the field inside of a sphere
    %for an incident plane wave
    
    %compute the argument for the spherical bessel function k*n*r
    knr = wavNum * n .* r;
    
    %compute the spherical bessel function
    jl_knr = sphbesselj(numOrd,knr,'multiple');
    
    E_s = zeros(simRes,simRes);
    E_i = E_s;
    %h = waitbar(0, 'Monte Carlo integration...');
    
    
    % r = sqrt(normKvec(1)^2 + normKvec(2)^2 + normKvec(3)^2);
    % theta =  atan2(normKvec(2), normKvec(1));
    % phi = acos(normKvec(2)/r);
    % kSpherical = [r theta phi]
    
    %the vector from the focal point to the center of the sphere
    c = ps - pf;
    display('Starting Monte Carlo simulation...')
    
    
    parfor i=1:samples
        cos_theta = (k_j(:,i)'*normRvecs')';
        cos_theta = reshape(cos_theta, simRes, simRes);
        %compute the legendre polynomials needed for computing E_s and E_i
        Pl_costheta = myLegendre(numOrd,cos_theta);
        hlkr_Plcostheta = hl_kr.*Pl_costheta;
        
        %multiply by the legendre polynomial
        jlknr_Plcostheta = jl_knr.*Pl_costheta;
        
        %multiply by the scattering coefficients B
        %sum all orders
        phase = exp(1i.*wavNum.*k_j(:,i)'*c');
        E_s = E_s + (1/samples).*subA.*phase.*sum(bsxfun(@times, hlkr_Plcostheta, reshape(B,[1 1 numOrd+1])),3);
        %E_s(r<A) = 0;%E_i(r<a);
        %subplot(1,2,1), imagesc((abs((E_s)))),title('E_s'), colorbar, axis image, colormap(brewer)
        E_i = E_i + (1/samples).*subA.*phase.*sum(bsxfun(@times, jlknr_Plcostheta, reshape(A,[1 1 numOrd+1])),3);
        %E_i(r>A) = 0;
        %subplot(1,2,2), imagesc((abs((E_i)))),title('E_i'), colorbar, axis image, colormap(brewer)
        %pause(0.1)
        %waitbar(i/samples);
    end
    display('Finished Monte Carlo simulation...')
    
    %close(h)
    
    E_s(psMask<a) = 0;%E_i(r<a);
    
    E_i(psMask>a) = 0;
    
    E_t = zeros(simRes, simRes);
    E_t(psMask<a) = E_i(psMask<a);
    E_t(psMask>=a) = E_f(psMask>=a) + E_s(psMask>=a);
    % this is not correct if the slice goes through the sphere
    E_d = fft2(E_t);
    E_d = fftshift((E_d));
    E_d(BPF) = 0;
    iftE_d = ifft2(ifftshift(E_d));
    %first crop the filtered near-field image of the source and scattered fields
    cropEd=iftE_d(startIdx:endIdx,startIdx:endIdx);
    cropEd(isnan(cropEd))=0;     cropEd(isinf(cropEd))=0;
    cropEf=E_f(startIdx:endIdx, startIdx:endIdx);
    cropEf(isinf(cropEf))=0;     cropEf(isinf(cropEf))=0;
    %save the absolute values of the cropped detected total field and the
    %detected focused field
    D_Ed = D_Ed + (abs(cropEd)).^2;
    D_Ef = D_Ef + (abs(cropEf)).^2;
    
end
toc
delete(POOL);
%calculate absorbance
A = -log10(D_Ed./D_Ef);

if displaySubplots==0
    brewer = brewermap(1000);
    %     figure(1),
    %     imagesc((abs((E_f)))),title(sprintf('abs(E_f) for p = %i',p)), colorbar, axis image
    %     colormap(brewer)
    %     figure(2),
    %     imagesc((abs((E_s)))),title(sprintf('abs(E_s) for p = %i',p)), colorbar, axis image
    %     colormap(brewer)
    %     figure(3),
    %     imagesc((abs((E_i)))),title(sprintf('abs(E_i) for p = %i',p)), colorbar, axis image
    %     colormap(brewer)
    %     figure(4),
    %     imagesc((abs((E_t)))),title(sprintf('abs(E_t) for p = %i',p)), colorbar, axis image
    %     colormap(brewer)
    %     figure(5),
    %     imagesc(abs(D_Ed)), title(sprintf('D_{Ed} at p = %i',p)),axis image, colormap(brewer), colorbar
    %     figure(6),
    %     imagesc(abs(D_Ef)), title(sprintf('D_{Ef} at p = %i',p)),axis image, colormap(brewer), colorbar
    figure(7),
    imagesc((A)), title(sprintf('A after p = %i',p)),axis image, colormap(brewer), colorbar
end
display('Done.')
