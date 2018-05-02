function showPointSourceSamples(rho, numPoints, amplitude)

%show sampling of point sources
figure
axis off
scatter(0,0,[],[1 1 0],'filled')
hold on
numRad = numel(rho);
if isscalar(numPoints)
    numPoints = ones(numRad)*numPoints;
    numPoints(1) = 1;
end

for i=1:numRad
    theta = 2*pi/numPoints(i):2*pi/numPoints(i):2*pi;
    [x,y] = pol2cart(theta,rho(i));
    scatter(x,y,[],[1 1 0]* amplitude(i),'filled')
end
axis square
set(gca,'color','black')
