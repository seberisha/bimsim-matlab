clear

% Figure 1(c)

% p=4 (4th order), R=3 (truncation radius), Nr=256 (# of sample points),
% eps_roots=1e-13 (not that it matters, but the bessel_zeros function is
% supposed to use this to produce more or less accurate results—doesn't seem to
% change anything).
pOrder = 4; % change to 1 to get Figure 1(a)
R = 3;
Nr = 256;
h = hankel_matrix(pOrder, R, Nr, 1e-13);

% parameter for the test function, a sinc
gamma = 5;
% Our sinc function includes pi: sinc(x) = x == 0 ? 1 : sin(pi * x) / (pi * x)
f1fun = @(r) sinc(2*gamma*r);

% Samples of f1(r) are on the zeros of Bessel function
r = h.J_roots / (2*pi*h.vmax);
% Transformed vector's sample points
v = h.J_roots / (2*pi*h.rmax);

% The transforms
HT  = @(f1) (h.T  * (f1 ./ h.J * h.rmax)) .* h.J / h.vmax;
IHT = @(f2) (h.T' * (f2 ./ h.J * h.vmax)) .* h.J / h.rmax;

f2 = HT(f1fun(r));

figure; plot(v, f2, '.-');
title(sprintf('HT p = %d', pOrder))
xlim([0 20])
xlabel('spatial frequency \nu')
grid