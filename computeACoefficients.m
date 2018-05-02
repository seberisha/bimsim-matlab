function A = computeACoefficients(k,a,lambda,n)
%Description:
%
%

%Assumptions:

%Input: 

%   k - magnitudes of vectors of the propogating plane wave
%   a - sphere radius
%   lambda - wavelength
%   n - refractive index

%Output:

%Author: S. Berisha
%Last modified: 09/15/15

%% the maximum order required for convergence
N_l = maxOrder(p, p_s, lambda);

%% precompute spherical Bessel functions, spherical Hankel functions,and their derivatives
nu = 0:N_l; 
ka = k*a;

j_ka = sphbesselj(nu,ones(size(nu))*ka);

%% first derivative of the spherical Bessel function of the first kind
 j = sym('sqrt(1/2*pi/x)*besselj(n+1/2,x)');
 dj = simplify(diff(j));
 dj = vectorize(inline(char(dj),'n','x'));
 kna=k*n*a;
 
 deriv_j_kna = dj(nu, kna);

 j_kna = sphbesselj(nu,ones(size(nu))*kna);
 deriv_j_ka = dj(nu, ka);
 
 
 h_ka = shank1(nu,ka);
 
 %% first derivative of the spherical Hankel function of the first kind
 h = sym('((pi/(2*x))^.5)*besselj(n+.5,x)+i*((pi/(2*x))^.5)*bessely(n+.5,x)');
 dh = simplify(diff(h));
 dh = vectorize(inline(char(dh),'n','x'));
 
 deriv_h_ka = dh(nu, ka);
 
 A = ((2*l + 1)*i.^nu ).*(j_ka.*deriv_j_ka - j_kna.*deriv_j_ka)/(j_kna.*deriv_h_ka - h_ka*deriv_j_kna*n);