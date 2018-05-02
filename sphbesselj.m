function js = sphbesselj(order,x,type)
%   Description:    Evalaute the spherical bessel functions of the first
%                    kind.
%
%   Input:  -order - postive real number or vector of orders of the functions
%           -x - scalar, vector, or matrix at which the spherical bessel
%               functions will be evaluated
%           -type - string: 
%                       -'multiple' - if order is a scalar and we want
%                           to evaluate spherical bessel functions of orders 0:order
%                       -'one' - if order is a scalar and we want to evaluate spherical 
%                               bessel functions only at that particular
%                               order
%
%   Output: -js - result of evaluating the spherical bessel functions at x
%                   and of orders 'order'. This can be a scalar, vector,
%                   matrix, or tensor.
switch type
    case 'multiple'
        if isscalar(x)
            js = sqrt(pi./(2.*x)).*besselj((0:order)'+0.5,x);
        elseif isvector(x)
            js = zeros([numel(x) order+1]);
            for i=0:order
                js(:,i+1) = sqrt(pi./(2*x)).*besselj(i+0.5,x);
            end
        else
            js = zeros([size(x) order+1]);
            for i=0:order
                js(:,:,i+1) = sqrt(pi./(2*x)).*besselj(i+0.5,x);
            end
        end
    case 'one'
        js = sqrt(pi./(2*x)).*besselj(order+0.5,x);
    otherwise
        display('no type specified - sphbesselj.m')
end
end