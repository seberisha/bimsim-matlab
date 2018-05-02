function E_s = computeInternalField(r, k, p_s, p_f)
%Description:
%

%Assumptions:
%   -Incident light is coherent: |k_j|=k, for all j \in J.

%Input: 
%   r - 
%   p - position of the point in EM?
%   p_s - center of the sphere
%   p_f - center of the focal point
%   lambda - wavelength
%   NA_c - condenser aperture
%   k - magnitudes of vectors of the propogating plane wave
%   a - sphere radius
%   alpha_1 - the inner angle

%Output:
%   B - the coefficients that couple the internal and external scattered
%   fields

%Authors: S. Berisha
%Last modified: 09/15/15



%% phase shift
c = p_s - p_f;
e_ikc = exp(i*k.*c);
%% get the coefficients
A = computeACoefficients(k,a,lambda,n);

%% the maximum order required for convergence
N_l = maxOrder(p, p_s, lambda);

%% compute the hankle function
nu = 0:N_l;
h_ka = shank1(nu,ka);

%% precompute Legendre polynomials
cos_theta = k'*r/(norm(k)*norm(r));
P_ct = legendre(N_l,cos_theta);

