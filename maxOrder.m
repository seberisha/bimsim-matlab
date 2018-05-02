function N_l = maxOrder(r, lambda)
%
%
%Input:
%   p - point in EM?
%   p_s - position of sphere 


N_l = ceil(2*pi.*r./lambda + 4*(2*pi.*r./lambda).^(1/3) + 2);

