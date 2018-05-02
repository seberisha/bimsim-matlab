
clear
addpath(genpath('~/source/stim-matlab/'))
datDir = '~/data/bimsim/matlab/bulkPMMA_6.5';

%% setup parameters and call forward gpu BimSim
params.interpolate=1;
params.a = 2.05;                     %radius of the sphere
%params.fov = round(params.a)*4; %field of view in micrometers
highMag = 1;

if highMag
    pixRes = 1.1;
else
    pixRes = 5.5;
end

params.res = 128;
params.fov = round(params.res*pixRes);

%specify padding
params.padding = 1;

%simRes = params.res*(2*params.padding + 1);

params.rad = 1:4;
params.numPS=4;
params.s=rng;
%load pmma.mat
%load([datDir '/measPMMA_wavRi.mat'])
%params.material=material;
%params.numWav = size(material,1);
params.material = [1e4/3502.09657 1.47634 0.000258604];

%params.material = [2.8554    1.4735    0.0076];


params.numWav=1;
params.ps = [0 0 0];

params.samples=400;
params.orderEf=100;

params.E0 = 1;
params.NA_in = .2;
params.NA_out = 0.62;
%this gives the kVec direction from the y axis
params.theta=1.5708;
params.phi=0;
params.fraction=8;

%% point sources
psMat = zeros(params.res);
sz = size(psMat);

if params.numPS>=1
    numRad = numel(params.rad);
    %[pf_x, pf_z, pf_theta] = pointSources(params.rad, params.numPS);
    halfFov = params.fov/2;
    figure
    numPS = 1;
    for i=1:numRad
        numPoints = (params.numPS)*i;
        theta=linspace(0,2*pi-0.01,numPoints);
        rho=ones(1,numPoints)*params.rad(i);
        [x,y] = pol2cart(theta,rho);
     
%         mat_x = round( axes2pix(sz(2), [1 sz(2)], x) );
%         mat_y = round( axes2pix(sz(1), [1 sz(1)], y) );
%          for j=1:numel(x)
%              psMat(round(mat_x(j)), round(mat_y(j))) = 255;
%          end
        plot(x,y,'w*','linewidth',2);
        hold on
        halfFov = halfFov - 1;
        numPS = numPS + numPoints;
    end
    xlim([-params.res/2 params.res/2])
    ylim([-params.res/2 params.res/2])
    axis square
    set(gca,'color','black')
end


%% run bim sim
output = cpuComputeBimSim(params);

%% display results
brewer = brewermap(1000);
wavenumbers = params.material(:,1);
wavenumbers = 1e4./wavenumbers;%wavenumbers in microns

%%
% figure;ax=axes;
% plot(wavenumbers,output.absSpec)
% set(ax, 'Xdir','reverse')