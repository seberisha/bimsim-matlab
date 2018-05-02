function sh1 = shank1(order,x,type)
%Description:   
%   Evaluate the hankel functions of the first kind.
%
%Input:         
%   -order = postive real number or vector of orders of the functions
%   -x = scalar, vector, or matrix at which the spherical bessel functions will be evaluated
%   -type = string:
%               -'multiple' = if order is a scalar and we want to evaluate spherical bessel functions of orders 0:order
%               -'one' = if order is a scalar and we want to evaluate 
%                        spherical  bessel functions only at that particular
%                        order
%
%Output: 
%   -sh1 = result of evaluating the hankel functions of the first kind at x
%          and of orders 'order'. This can be a scalar, vector, matrix, or tensor.

switch type
    case 'multiple'
        if(isscalar(x))
            orderVec=(0:order)';
            sh1=(sqrt(pi./(2.*x))).*(besselj(orderVec + 0.5, x) + 1i.*bessely(orderVec + 0.5, x));
        elseif (isvector(x))
            sh1 = zeros([numel(x) order + 1]);
            for i=0:order
                sh1(:,i+1)=(sqrt(pi./(2.*x))).*(besselj(i + 0.5, x) + 1i.*bessely(i + 0.5, x));
            end
        else
            sh1 = zeros([size(x) order + 1]);
            for i=0:order
                sh1(:,:,i+1)=(sqrt(pi./(2.*x))).*(besselj(i + 0.5, x) + 1i.*bessely(i + 0.5, x));
            end
        end
    case 'one'
        sh1=(sqrt(pi./(2.*x))).*(besselj(order + 0.5, x) + 1i.*bessely(order + 0.5, x));
    otherwise
        display('shank1: incorrect type');
end