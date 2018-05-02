%%
for i=1:size(t,3)
    subplot(121), imagesc(t(:,:,i)),title(sprintf('%f',wav(i))),axis image,colorbar,colormap(brewer)
    subplot(122), imagesc(absImages(:,:,i)),title(sprintf('%f',wav(i))),axis image,colorbar,colormap(brewer)
    pause(1e-40)
end

%%
 int1 = awgn(int(:,:,end),45);
 int1(int1<0)=0;
 inc1 = awgn(inc(:,:,end),45);
 inc1(inc1<0)=0;
 t=-log10((int1)./(inc1));
 imagesc((t)), axis image, colormap(brewer),colorbar
 
%%
int1 = awgn(int(:,:,1),40);
int1(int1<0)=0;
inc1 = awgn(inc(:,:,1),40);
inc1(inc1<0)=0;
t=-log10((int1)./(inc1));
imagesc((t)), axis image, colormap(brewer),colorbar
ftirA = absImages(:,:,1);
caxis([min(ftirA(:)) max(ftirA(:))])


%%
int1 = awgn(int(:,:,end),40);
int1(int1<0)=0;
inc1 = awgn(inc(:,:,end),40);
inc1(inc1<0)=0;
t=-log10((int1)./(inc1));
imagesc((t)), axis image, colormap(brewer),colorbar
ftirA = absImages(:,:,end);
caxis([min(ftirA(:)) max(ftirA(:))])


%%
ftirA = absImages(:,:,1);
t = awgn(A_b(:,:,1),40);
imagesc((t)), axis image, colormap(brewer),colorbar
caxis([min(ftirA(:)) max(ftirA(:))])
%%
t = awgn(A_b(:,:,end),50);
imagesc((t)), axis image, colormap(brewer),colorbar
ftirA = absImages(:,:,end);
caxis([min(ftirA(:)) max(ftirA(:))])

%% add noise to intensity and incident field images
% Add white Gaussian noise E the the blurred image, scaled such that
% || e ||_2 / || A x ||_2 = 0.01.
nl = 0.9;
E = randn(128,128);
E = E / norm(E,'fro');
numWav = 661;

for i=1:1
    %I_e(:,:,i) = int(:,:,i) + nl*norm(int(:,:,i),'fro')*E;
    %I0_e(:,:,i) = inc(:,:,i) + nl*norm(inc(:,:,i),'fro')*E;
    A_e(:,:,i) = A_b(:,:,i) + nl*norm(A_b(:,:,i),'fro')*E;
end
%I_e(I_e<0)=0; 
%I0_e(I0_e<0)=0;
%ratio = (I_e)./(I0_e);
%ratio(ratio<0)=0;
%A_e = real(-log10((ratio)));

%A_e(A_e<0)=0;

imagesc((A_e(:,:,1))), axis image, colormap(brewer),colorbar



%%
nl = 0.01;
E = randn(size(int));


for i=1:size(E,3)
    E(:,:,i) = E(:,:,i) / norm(E(:,:,i),'fro');
    I_e(:,:,i) = int(:,:,i) + nl*norm(int(:,:,i),'fro')*E(:,:,i);
    I0_e(:,:,i) = inc(:,:,i) + nl*norm(inc(:,:,i),'fro')*E(:,:,i);
end
%I_e(I_e<0)=0; 
%I0_e(I0_e<0)=0;
ratio = (I_e)./(I0_e);
%ratio(ratio<0)=0;
A_e = real(-log10((ratio)));

%A_e(A_e<0)=0;

imagesc(real(A_e(:,:,1))), axis image, colormap(brewer),colorbar

%% test for fft of circle
clear all;

samples = 128;

%radius of the circle
r_fraction = 6;
rad = floor(samples/r_fraction);

%generate the meshgrid
[X, Y] = meshgrid(-(samples-1)/2 : (samples-1)/2);
RHO = sqrt(X.^2 + Y.^2);

%generate the circle
C = double(RHO < rad);

%filter
C_f = fft2(C);
disp(C_f( ceil(samples/2), ceil(samples/2)));
C_f2 = C_f .* ifftshift(C);
C2 = ifft2(C_f2);

%----1D case

%generate the meshgrid
RHO_1D = 0:(samples-1)/2;

%generate the circle
C_1D = double(RHO_1D < rad);
C_1D = C_1D';


% Figure 1(c)

% p=4 (4th order), R=3 (truncation radius), Nr=256 (# of sample points),
% eps_roots=1e-13 (not that it matters, but the bessel_zeros function is
% supposed to use this to produce more or less accurate results—doesn't seem to
% change anything).
pOrder = 4; % change to 1 to get Figure 1(a)
R = length(RHO_1D);
Nr = numel(C_1D);
h = hankel_matrix(pOrder, R, Nr, 1e-13);


% Samples of f1(r) are on the zeros of Bessel function
r = h.J_roots / (2*pi*h.vmax);
% Transformed vector's sample points
v = h.J_roots / (2*pi*h.rmax);

% The transforms
HT  = @(f1) (h.T  * (f1 ./ h.J * h.rmax)) .* h.J / h.vmax;
IHT = @(f2) (h.T' * (f2 ./ h.J * h.vmax)) .* h.J / h.rmax;

f2 = HT(C_1D);

plot(v, f2, '.-');
title(sprintf('HT p = %d', pOrder))
%xlim([0 20])
xlabel('spatial frequency \nu')
grid
disp(f2(1));

%%
%perform the Hankel transform
[H,k,r,I,K,R,h]=dht([],length(RHO_1D), length(RHO_1D));
H = I*(C_1D./R).*K;
%C_1D_h = dht(C_1D);

%%

C_1D_h = I*(C_1D./R).*K;
disp(C_1D_h(1));

%filter
C_1D_h2 = C_1D_h .* C_1D;



imagesc(C_1D_h);


%%
R = sqrt((x - 64).^2 + (y - 64).^2);

%create the circle
C = double(R < rad);

%take an FFT of the signal
Cf = fft2(C);

%multiply the signal by the circle
Cf = Cf .* ifftshift(C);

%invert the transform
Ci = ifft2(Cf);


%create the radial signal
rho = 0:63

imagesc(Ci);





%%


[x,y] = meshgrid(-64:63,-64:63);
idx = sqrt(x.^2 + y.^2);
C = zeros(128,128);
C(idx<32) = 1;
subplot(2,3,1), mesh(C),title('circle')
fftC = fft2(C);
subplot(2,3,2),mesh(abs(fftshift(fftC))),title('fft2 of circle')
cfilter = zeros(128,128);
cfilter(idx<16)=1;
subplot(2,3,3),mesh(cfilter), title('circle filter')
filterFC = fftshift(fftC).*filter;
subplot(2,3,4),mesh(abs(fftshift(fftC))), title('filtered fft2 of circle')
iFilterFC = ifft2((filterFC));
subplot(2,3,5),mesh(abs((iFilterFC))), title('ifft2 of filtered fft2 of circle')

%% test for hankel of box
figure
b = zeros(128,1);
idx = -64:63;
b(idx>-32 & idx<32) = 1;
subplot(2,3,1), plot(b), title('rect')
[H,k,rad,I,K,R,h]=dht([],64,128);
hankelB = I*(b./R).*K;
subplot(2,3,2), plot(hankelB), title('hankel of rect')
filterB = zeros(128,1);
filterB(idx>=0 & idx<16) = 1;
subplot(2,3,3), plot(filterB), title('rect filter')
filterHC = hankelB.*filterB;
subplot(2,3,4), plot(filterHC), title('filtered hankel of rect')
iFilterHC = idht(filterHC,I,K,R);
subplot(2,3,5), plot(iFilterHC), title('iHankel of filtered hankel of rect')

%%

fil = (gaussianbpf_1d(99,objMinParam,objMaxParam,(fmag)));
PSF = exp( -(fmag.^2) / (2*objMaxParam^2) );

%%
nx=99;
ny=99;
d1=120;
d0=30;
% Initialize filter.
filter1 = ones(2*nx-1,2*ny-1);
filter2 = ones(2*nx-1,2*ny-1);
filter3 = ones(2*nx-1,2*ny-1);

is = linspace(-nx,nx,2*nx-1);
[ix,iy]=meshgrid(is,is);

dist = ((ix).^2 + (iy).^2).^.5;
filter1 = exp(-dist.^2./(2*d1^2));
filter2 = exp(-dist.^2./(2*d0^2));
filter3 = 1.0 - filter2;
filter3 = filter1.*filter3;


%%
nx=99;
ny=99;
% d1=120;
% d0=30;
d1=0.0526;d0=0.017;
% Initialize filter.
filter1 = ones(2*nx-1,2*ny-1);
filter2 = ones(2*nx-1,2*ny-1);
filter3 = ones(2*nx-1,2*ny-1);

for i = 1:2*nx-1
    for j =1:2*ny-1
        dist = ((i-(nx+1))^2 + (j-(ny+1))^2)^.5;
        % Use Gaussian filter.
        filter1(i,j) = exp(-dist^2/(2*d1^2));
        filter2(i,j) = exp(-dist^2/(2*d0^2));
        filter3(i,j) = 1.0 - filter2(i,j);
        filter3(i,j) = filter1(i,j).*filter3(i,j);
    end
end

%%
% 
nx=256;ny=256;
%d1=0.0526;d0=0.017;
d0=30; d1 = 120;
for i = 1:2*nx-1
    for j =1:2*ny-1
        dist = ((i-(nx+1))^2 + (j-(ny+1))^2)^.5;
        % Use Gaussian filter.
        filter1(i,j) = exp(-dist^2/(2*d1^2));
        filter2(i,j) = exp(-dist^2/(2*d0^2));
        filter3(i,j) = 1.0 - filter2(i,j);
        filter3(i,j) = filter1(i,j).*filter3(i,j);
    end
end



%%

% Initialize filter.
fs = 1/(2*params.gridSize);
fc = linspace(-fs/2,fs/2,params.simRes);
absFC = abs(fc);

filter1 = ones(params.simRes,1);
filter2 = ones(params.simRes,1);
filter3 = ones(params.simRes,1);

dist = absFC;
filter1 = exp(-absFC./(2*objMaxParam^2));
filter2 = exp(-absFC./(2*objMinParam^2));
filter3 = 1 - filter2;
filter3 = filter1.*filter3;

% 
% for i = 1:2*nx-1
%     for j =1:2*ny-1
%         dist = ((i-(nx+1))^2 + (j-(ny+1))^2)^.5;
%         % Use Gaussian filter.
%         filter1(i,j) = exp(-dist^2/(2*d1^2));
%         filter2(i,j) = exp(-dist^2/(2*d0^2));
%         filter3(i,j) = 1.0 - filter2(i,j);
%         filter3(i,j) = filter1(i,j).*filter3(i,j);
%     end
% end


%%
fs = 1/(2*params.gridSize);
fc = linspace(-fs/2,fs/2,params.simRes);

fEt = fftshift(fft(E_t));
absFC = abs(fc);
kernelBP = absFC > objMaxParam & absFC < objMinParam;
bpfEt = kernelBP.*fEt;
figure, plot(fc, bpfEt)

ibpfEt = ifft(ifftshift(bpfEt));
figure, plot(gridPoints, real(ibpfEt))


%%
t = bsxfun(@times, hlkr_Plcostheta, reshape(params.B,[1 1 params.numOrd+1]));
t50 = squeeze(t(50,:,:));
hp50 = squeeze(hlkr_Plcostheta(50,:,:));
test_2d = zeros(size(t50));
for i=1:99
   test_2d(i,:) = hp50(i,:).*(squeeze(params.B)).'; 
end


%%
for i=1:99
   t(i,:) = hlkr_Plcostheta(i,:).*params.B'; 
end

%%
params.ps = [0 0 0 ];
params.a = 6.5;                     %radius of the sphere

params.fov = round(params.a)*3; %field of view in micrometers

padding = 1; 
params.gridSize = round(params.fov/2)*(2*padding + 1);
params.res = 33;
params.simRes = params.res*(2*padding + 1);


gridPoints = (2*params.gridSize)*(0:params.simRes-1)/params.simRes - params.gridSize;


[x_2d,z_2d] = meshgrid(gridPoints, gridPoints); % field slice in the x z plane
y_2d = ones(params.simRes,params.simRes)*(params.a);   %field plane y = 0

params.rVecs_2d = zeros(params.simRes*params.simRes, 3);
params.rVecs_2d(:,1) = x_2d(:); params.rVecs_2d(:,2) = y_2d(:); params.rVecs_2d(:,3) = z_2d(:); %r value at each pixel position
params.psVecs_2d = bsxfun(@minus, params.rVecs_2d, params.ps);
normPMinPs_2d = sqrt(sum(params.psVecs_2d.^2,2));
params.r_2d = reshape(normPMinPs_2d, params.simRes, params.simRes);

r50 = params.r_2d(50,:);
tvecs(:,1) = x_2d(50,:)'; tvecs(:,2) = y_2d(50,:)'; tvecs(:,3) = z_2d(50,:)';
ntvecs = sqrt(sum(tvecs.^2,2));

%%
E_t(1:49,:)=0;
E_t(51:end,:)=0;

%%
tt = zeros(size(E_t));

for i=0:0.5:180
    er = imrotate(E_t,i,'crop');
    subplot(1,2,1), imagesc(abs(er)),axis image,colormap(brewer),colorbar
    tt(er~=0)= er(er~=0);
    subplot(1,2,2),imagesc(abs(tt)),axis image,colormap(brewer),colorbar
    pause(1e-7)
end

