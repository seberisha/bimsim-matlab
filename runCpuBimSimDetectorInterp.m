%% setup parameters and call forward gpu BimSim

clear
%gpuDevice(2)
addpath(genpath('~/source/stim-matlab/'));
datDir = '~/data/bimsim/matlab/bulkPMMA_6.5'
addpath(genpath('~/source/stim-matlab/'))

params.res = 32;
sigma = 11;
params.numPS = 94;

params.interpolate=1;
params.a = 2.25;                     %radius of the sphere
%params.fov = round(params.a)*4; %field of view in micrometers
highMag = 1;

if highMag
    pixRes = 1.1;
else
    pixRes = 5.5;
end

params.material = [1e4/3502.09657 1.47634 0.000258604];
params.spacing=1;
%params.spacing = floor(params.material(1));

params.fov = round(params.res*pixRes);

params.rad = 1:round(params.res/2)-1;
%params.rad = params.spacing:params.spacing:round(params.fov/2);
%params.numPS = round(params.fov/2);



params.s=rng;
%load pmma.mat
%load([datDir '/measPMMA_wavRi.mat'])
%params.material=material;
%params.numWav = size(material,1);

%params.material = [2.8554    1.4735    0.0076];


params.numWav=1;
params.ps = [0 0 0];

params.samples=400;
params.orderEf=50;
%specify padding
params.padding = 1;
params.E0 = 1;
params.NA_in = .34;
params.NA_out = 0.62;
%this gives the kVec direction from the y axis
params.theta=1.5708;
params.phi=0;
%params.fraction=8;

%% point sources
psMat = zeros(params.res);
sz = size(psMat);
numRad = numel(params.rad);

t = exp(-(([0:numRad].^2))/(2*(sigma^2)))/(2*pi*(sigma^2));
amplitude = exp(-(([1:numRad].^2))/(2*(sigma^2)))/(2*pi*(sigma^2));
amplitude=amplitude./max(amplitude(:));
%amplitude = changeRange(amplitude, 0, 1);
figure

rho = params.rad(end);
numPoints = round(((2*pi*rho-params.spacing)/params.spacing + 1));
maxRange = floor(2*pi*rho-params.spacing/2);
s = linspace(0,maxRange,numPoints);
theta = s./rho;

scatter(0,0,[],[1 1 0],'filled')
hold on
if params.numPS>1
    
    %[pf_x, pf_z, pf_theta] = pointSources(params.rad, params.numPS);
    %numPoints = round(2*pi*params.rad(end));
    radiusIdx=1;
    for i=params.spacing:params.spacing:numRad
        %theta=linspace(0,2*pi*i-0.1,numPoints);
        rho=params.rad(radiusIdx);
        %s = 0:numPoints;
        [x,y] = pol2cart(theta,rho);
        scatter(x,y,[],[1 1 0]* amplitude(i),'filled')
        radiusIdx = radiusIdx + 1;
        
    end
    axis square
    set(gca,'color','black')
end

params.amplitude = amplitude(1:params.spacing:end);

%% run bim sim 
% version -- interpolate the near fields

% version -- interpolate intensities at the detector
output = cpuBimSimInterpAtDetector(params);

return

%% add noise to intensity and incident field images
% Add white Gaussian noise E the the blurred image, scaled such that
% || e ||_2 / || A x ||_2 = 0.01.
% nl = 0.01;
E = randn(params.res,params.res);
E = E / norm(E,'fro');
numWav = 661;

for i=1:numWav
    I_e(:,:,i) = int(:,:,i) + nl*norm(int(:,:,i),'fro')*E;
    I0_e(:,:,i) = inc(:,:,i) + nl*norm(inc(:,:,i),'fro')*E;
end
%I_e(I_e<0)=0;
%I0_e(I0_e<0)=0;
ratio = (I_e)./(I0_e);
%ratio(ratio<0)=0;
A_e = real(-log10((ratio)));

%A_e(A_e<0)=0;

imagesc(real(A_e(:,:,1))), axis image, colormap(brewer),colorbar

%% display results
brewer = brewermap(1000);
wavenumbers = params.material(:,1);
wavenumbers = 1e4./wavenumbers;%wavenumbers in microns

%%
% figure;ax=axes;
% plot(wavenumbers,output.absSpec)
% set(ax, 'Xdir','reverse')