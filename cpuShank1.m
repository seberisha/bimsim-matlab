function sh1 = cpuShank1(order,x,type)

switch type
    case 'multiple'
        if(isscalar(x))
            orderVec=(0:order)';
            sh1=(sqrt(pi./(2.*x))).*(besselj(orderVec + 0.5, x) + 1i.*bessely(orderVec + 0.5, x));   
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