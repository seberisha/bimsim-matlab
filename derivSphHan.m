function h_p = derivSphHan(order, x)
%   Description:    Evalaute the derivative of the hankel functions of the first
%                    kind.
%
%   Input:  -order - postive real number or vector of orders of the functions
%           -x - scalar, vector, or matrix at which the spherical bessel
%               functions will be evaluated
%
%   Output: -h_p - result of evaluating the derivative of the hankel functions 
%                   of the first kind at x and of orders 'order'. This can be a scalar, vector,
%                   matrix, or tensor.

orderVec=(0:order)';
sh_n = ((shank1(order,x,'multiple')));
sh_n_m_1 = shank1(orderVec-1,x,'one');
sh_n_p_1 = shank1(orderVec+1,x,'one');
h_p = (1/2).*(sh_n_m_1 - (sh_n + x.*sh_n_p_1)./x);