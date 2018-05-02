

%%

clear 

lp = 0.01; hp = 0.05;
samples = 128;

% p=4 (4th order), R=3 (truncation radius), Nr=256 (# of sample points),
% eps_roots=1e-13 (not that it matters, but the bessel_zeros function is
% supposed to use this to produce more or less accurate results—doesn't seem to
% change anything).
pOrder = 4; % change to 1 to get Figure 1(a)
R = samples/2;
Nr = samples/2;
h = hankel_matrix(pOrder, R, Nr, 1e-13);


% Samples of f1(r) are on the zeros of Bessel function
r = h.J_roots / (2*pi*h.vmax);
% Transformed vector's sample points
v = h.J_roots / (2*pi*h.rmax);

filter_1d = (v < lp | v > hp);

%radius of the circle
r_fraction = 6;
rad = floor(samples/r_fraction);

%generate the meshgrid
[X, Y] = meshgrid([-wrev(r).' r.']);
RHO = sqrt(X.^2 + Y.^2);

%generate the circle
C = double(RHO < rad);

df = 1/(samples);

[iv, iu] = meshgrid(0:samples-1, 0:samples-1);
u=zeros(size(iu)); v = u;
idx = find(iu <= samples/2);
u(idx) = iu(idx);
idx = find(iu > samples/2);
u(idx) = (iu(idx) - samples+1);
u=u.*df;
idx = find(iv <= samples/2);
v(idx) = iv(idx);
idx = find(iv > samples/2);
v(idx) = iv(idx) - samples + 1;
v=v.*df;

fmag = sqrt(u.*u + v.*v);    %compute the magnitude of the frequencies

filter_2d = (fmag < lp | fmag > hp);


%filter
C_f = fft2(C);
disp(C_f( ceil(samples/2), ceil(samples/2)));
C_f(filter_2d)=0;
C2 = ifft2(C_f);

%----1D case

%generate the meshgrid
RHO_1D = RHO(samples/2 + 1,samples/2 + 1:end);

%transform to polar coordinates
rho = fmag;
phi = atan(u./v);
u_p = rho.*cos(phi);
v_p = rho.*sin(phi);

fmag_p = sqrt(u_p.*u_p + v_p.*v_p);
shiftFmag_p = fftshift(fmag_p);

fmag_1d = shiftFmag_p(samples/2+1, samples/2+1:end);
fmag_1d(isnan(fmag_1d))=0;

%filter_1d = fmag_1d<lp | fmag_1d>hp;

%generate the circle
C_1D = double(RHO_1D < rad).';


% The transforms
HT  = @(f1) (h.T  * (f1 ./ h.J * h.rmax)) .* h.J / h.vmax;
IHT = @(f2) (h.T' * (f2 ./ h.J * h.vmax)) .* h.J / h.rmax;

C_1D_h = HT(C_1D);
C_1D_h(filter_1d)=0;
iC_1D_h_f = IHT(C_1D_h);

disp(C_1D_h(1));

Ci = interp1(RHO_1D, iC_1D_h_f, RHO,'linear','extrap');
C2(isnan(Ci))=0;
Ci(isnan(Ci))=0;
%plot 
figure
subplot(2,2,1), mesh(abs(C2)),title('ifft')
subplot(2,2,2), mesh(abs(Ci)),title('ihankel')
subplot(2,2,3), imagesc(abs(C2)), colorbar
subplot(2,2,4), imagesc(abs(Ci)), colorbar


%% evaluate hankel at r

clear 

lp = 0.01; hp = 0.05;
samples = 128;

% p=4 (4th order), R=3 (truncation radius), Nr=256 (# of sample points),
% eps_roots=1e-13 (not that it matters, but the bessel_zeros function is
% supposed to use this to produce more or less accurate results—doesn't seem to
% change anything).
pOrder = 4; % change to 1 to get Figure 1(a)
R = samples/2;
Nr = samples/2;
h = hankel_matrix(pOrder, R, Nr, 1e-13);


% Samples of f1(r) are on the zeros of Bessel function
r = h.J_roots / (2*pi*h.vmax);
% Transformed vector's sample points
v = h.J_roots / (2*pi*h.rmax);

filter_1d = (v < lp | v > hp);

%radius of the circle
r_fraction = 6;
rad = floor(samples/r_fraction);

%generate the meshgrid
[X, Y] = meshgrid([-wrev(r).' r.']);
RHO = sqrt(X.^2 + Y.^2);

%generate the circle
C = double(RHO < rad);

df = 1/(samples);

[iv, iu] = meshgrid(0:samples-1, 0:samples-1);
u=zeros(size(iu)); v = u;
idx = find(iu <= samples/2);
u(idx) = iu(idx);
idx = find(iu > samples/2);
u(idx) = (iu(idx) - samples+1);
u=u.*df;
idx = find(iv <= samples/2);
v(idx) = iv(idx);
idx = find(iv > samples/2);
v(idx) = iv(idx) - samples + 1;
v=v.*df;

fmag = sqrt(u.*u + v.*v);    %compute the magnitude of the frequencies

filter_2d = (fmag < lp | fmag > hp);


%filter
C_f = fft2(C);
disp(C_f( ceil(samples/2), ceil(samples/2)));
C_f(filter_2d)=0;
C2 = ifft2(C_f);

%----1D case

%generate the meshgrid
RHO_1D = r;

%generate the circle
C_1D = double(RHO_1D < rad);

% The transforms
HT  = @(f1) (h.T  * (f1 ./ h.J * h.rmax)) .* h.J / h.vmax;
IHT = @(f2) (h.T' * (f2 ./ h.J * h.vmax)) .* h.J / h.rmax;

C_1D_h = HT(C_1D);
C_1D_h(filter_1d)=0;
iC_1D_h_f = IHT(C_1D_h);

disp(C_1D_h(1));

Ci = interp1(RHO_1D, iC_1D_h_f, RHO);
C2(isnan(Ci))=0;
Ci(isnan(Ci))=0;
%plot 
figure
subplot(2,2,1), mesh(abs(C2)),title('ifft')
subplot(2,2,2), mesh(abs(Ci)),title('ihankel')
subplot(2,2,3), imagesc(abs(C2)), colorbar
subplot(2,2,4), imagesc(abs(Ci)), colorbar



%% hankel for a circle function - low pass filter
clear 

lp = 0.01; hp = 0.05;
samples = 128;

%radius of the circle
r_fraction = 6;
rad = floor(samples/r_fraction);

%generate the meshgrid
[X, Y] = meshgrid(-(samples-1)/2 : (samples-1)/2);
RHO = sqrt(X.^2 + Y.^2);

%generate the circle
C = double(RHO < rad);

df = 1/(samples);

[iv, iu] = meshgrid(0:samples-1, 0:samples-1);
u=zeros(size(iu)); v = u;
idx = find(iu <= samples/2);
u(idx) = iu(idx);
idx = find(iu > samples/2);
u(idx) = (iu(idx) - samples+1);
u=u.*df;
idx = find(iv <= samples/2);
v(idx) = iv(idx);
idx = find(iv > samples/2);
v(idx) = iv(idx) - samples + 1;
v=v.*df;

fmag = sqrt(u.*u + v.*v);    %compute the magnitude of the frequencies

filter_2d = (fmag < lp | fmag > hp);


%filter
C_f = fft2(C);
disp(C_f( ceil(samples/2), ceil(samples/2)));
C_f(filter_2d)=0;
C2 = ifft2(C_f);

%----1D case

%generate the meshgrid
RHO_1D = 0:(samples-1)/2;

%generate the circle
C_1D = double(RHO_1D < rad);
C_1D = C_1D';

fmag = fftshift(fmag);
fmag = fmag(samples/2+1,samples/2+1:end);

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
filter_1d = (v < lp | v > hp);

% The transforms
HT  = @(f1) (h.T  * (f1 ./ h.J * h.rmax)) .* h.J / h.vmax;
IHT = @(f2) (h.T' * (f2 ./ h.J * h.vmax)) .* h.J / h.rmax;

C_1D_h = HT(C_1D);
C_1D_h(filter_1d)=0;
iC_1D_h_f = IHT(C_1D_h);

disp(C_1D_h(1));

Ci = interp1(RHO_1D, iC_1D_h_f, RHO);
C2(isnan(Ci))=0;
Ci(isnan(Ci))=0;
%plot 
figure
subplot(2,2,1), mesh(abs(C2)),title('ifft')
subplot(2,2,2), mesh(abs(Ci)),title('ihankel')
subplot(2,2,3), imagesc(abs(C2)), colorbar
subplot(2,2,4), imagesc(abs(Ci)), colorbar

%% hankel for a circle function
clear 

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

C_1D_h = HT(C_1D);
C_1D_h_f = C_1D_h.*C_1D;
iC_1D_h_f = IHT(C_1D_h_f);

disp(C_1D_h(1));

Ci = interp1(RHO_1D, iC_1D_h_f, RHO);
C2(isnan(Ci))=0;
Ci(isnan(Ci))=0;
figure
subplot(2,2,1), mesh(abs(C2)),title('ifft')
subplot(2,2,2), mesh(abs(Ci)),title('ihankel')
subplot(2,2,3), imagesc(abs(C2)), colorbar
subplot(2,2,4), imagesc(abs(Ci)), colorbar
