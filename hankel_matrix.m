function s_HT = hankel_matrix(order, rmax, Nr, eps_roots);
%HANKEL_MATRIX: Generates data to use for Hankel Transforms
%
%	s_HT = hankel_matrix(order, rmax, Nr, eps_roots)
%
%	s_HT	=	Structure containing data to use for the pQDHT
%	order	=	Transform order
%	rmax	=	Radial extent of transform
%	Nr		=	Number of sample points
%	eps_roots	=	Error in estimation of roots of Bessel function (optional)
%
%	s_HT:
%		order, rmax, Nr =	As above
%		J_roots		=	Roots of the pth order Bessel fn.
%		J_roots_N1	=	(N+1)th root
%		r			=	Radial co-ordinate vector
%		v			=	frequency co-ordinate vector
%		kr			=	Radial wave number co-ordinate vector
%		vmax		=	Limiting frequency
%					=	roots_N1 / (2*pi*rmax)
%		S			=	rmax * 2*pi*vmax (S product)
%		T			=	Transform matrix
%		J			=	Scaling vector
%					=	J_(order+1){roots}
%
%	The algorithm used is that from:
%		"Computation of quasi-discrete Hankel transforms of the integer
%		order for propagating optical wave fields"
%		Manuel Guizar-Sicairos and Julio C. Guitierrez-Vega
%		J. Opt. Soc. Am. A 21(1) 53-58 (2004)
%
%	The algorithm also calls the function:
%	zn = bessel_zeros(1, p, Nr+1, 1e-6),
%	where p and N are defined above, to calculate the roots of the bessel
%	function. This algorithm is taken from:
%  		"An Algorithm with ALGOL 60 Program for the Computation of the
%  		zeros of the Ordinary Bessel Functions and those of their
%  		Derivatives".
%  		N. M. Temme
%  		Journal of Computational Physics, 32, 270-279 (1979)
%
%	Example: Propagation of radial field
%
%		% Note the use of matrix and element products / divisions
%		H = hankel_matrix(0, 1e-3, 512);
%		DR0 = 50e-6;
%		Ur0 = exp(-(H.r/DR0).^2);
%		Ukr0 = H.T * (Ur0./H.J);
%		k0 = 2*pi/800e-9;
%		kz = realsqrt((k0^2 - H.kr.^2).*(k0>H.kr));
%		z = (-5e-3:1e-5:5e-3);
%		Ukrz = (Ukr0*ones(1,length(z))).*exp(i*kz*z);
%		Urz = (H.T * Ukrz) .* (H.J * ones(1,length(z)));
%
%	See also bessel_zeros, besselj

if (~exist('eps_roots', 'var')||isempty(eps_roots))
	s_HT.eps_roots = 1e-6;
else
	s_HT.eps_roots = eps_roots;
end

s_HT.order = order;
s_HT.rmax = rmax;
s_HT.Nr = Nr;

%	Calculate N+1 roots:
J_roots = bessel_zeros(1, s_HT.order, s_HT.Nr+1, s_HT.eps_roots);
s_HT.J_roots = J_roots(1:end-1);
s_HT.J_roots_N1 = J_roots(end);

%	Calculate co-ordinate vectors
s_HT.r = s_HT.J_roots * s_HT.rmax / s_HT.J_roots_N1;
s_HT.v = s_HT.J_roots / (2*pi * s_HT.rmax);
s_HT.kr = 2*pi * s_HT.v;
s_HT.vmax = s_HT.J_roots_N1 / (2*pi * s_HT.rmax);
s_HT.S = s_HT.J_roots_N1;

%	Calculate hankel matrix and vectors
%	I use (p=order) and (p1=order+1)
Jp = besselj(s_HT.order, (s_HT.J_roots) * (s_HT.J_roots.') / s_HT.S);
Jp1 = abs(besselj(s_HT.order+1, s_HT.J_roots));
s_HT.T = 2*Jp./(Jp1 * (Jp1.') * s_HT.S);
s_HT.J = Jp1;