function [x,y, theta] = pointSources(rad, numPoints)
%origPoints = numPoints;
if isscalar(rad)
    %rad = floor(samples/r_fraction);
    theta=linspace(0,2*pi,numPoints);
    rho=ones(1,numPoints)*rad;
    [x,y] = pol2cart(theta,rho);
else
    numRad = numel(rad);
    x = zeros(numRad,numPoints);
    y = x;
    theta=linspace(0,2*pi-0.1,numPoints);
    
    for i=1:numRad
        %numPoints = round(origPoints/halfFov);
        numPoints = numPoints*i;
        rho=ones(1,numPoints)*rad(i);
        %x(i,:) = sqrt(rho).*cos(theta);
        %y(i,:) = sqrt(rho).*sin(theta);
        %[x(i,:),y(i,:)] = pol2cart(theta,rho);
        [x,y] = pol2cart(theta,rho);
        %halfFov = halfFov-1;
    end
end


