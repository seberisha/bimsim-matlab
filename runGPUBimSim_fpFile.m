
clear
addpath(genpath('~/source/stim-matlab/'));
datDir = '~/data/bimsim/matlab/bulkPMMA_6.5'

%% setup parameters and call forward gpu BimSim
params.interpolate=0;
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
%specify padding
params.padding = 1;
params.E0 = 1;
params.NA_in = .2;
params.NA_out = 0.62;
%this gives the kVec direction from the y axis
params.theta=1.5708;
params.phi=0;
params.fraction=8;

%% point sources
load fps2FullCircle.mat
params.fp = fps;

%% run bim sim
output = gpuBimSim_fpFile(params);

%% display results
% brewer = brewermap(1000);
% wavenumbers = params.material(:,1);
% wavenumbers = 1e4./wavenumbers;%wavenumbers in microns

%%
% figure;ax=axes;
% plot(wavenumbers,output.absSpec)
% set(ax, 'Xdir','reverse')