function j_p = derivSphBes(order, x)
%   Description:    Evalaute the derivative of the spherical bessel function
%                   of the first kind.
%
%   Input:  -order - postive real number or vector of orders of the functions
%           -x - scalar, vector, or matrix at which the spherical bessel
%               functions will be evaluated
%
%   Output: -j_p - result of evaluating the derivative of the spherical 3
%                  bessel functions at x and of orders 'order'. This can be a scalar, vector,
%                   matrix, or tensor.

js_n = (sphbesselj(order,x,'multiple'));
js_n_m_1 = (sphbesselj((0:order)'-1,x,'one'));
js_n_p_1 = (sphbesselj((0:order)'+1,x,'one'));
j_p = (1/2).*(js_n_m_1 - (js_n + x.*js_n_p_1)./x);   
   