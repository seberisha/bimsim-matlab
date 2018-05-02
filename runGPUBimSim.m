%% setup parameters 

clear
gpuDevice(2)
addpath(genpath('~/source/stim-matlab/'));
datDir = '~/data/bimsim/matlab/bulkPMMA_6.5'
addpath(genpath('~/source/stim-matlab/'))

params.res = 32;   %resolution at the detector

% standard deviation for a Gaussian extended source --to be estimated
% from FTIR interferogram data - for Agilent FTIR select show raw data in
% lancer control before start of data acquisition (when save dialog box pops up)
sigma = 11;    

params.interpolate=1;   %option to interpolate for extended sources
params.a = 2.25;    %radius of the sphere
highMag = 1;    %for high mag mode simulation choose 1, 0 for low mag

if highMag
    pixRes = 1.1;   %pixel resolution in um for high mag
else
    pixRes = 5.5;  %pixel resolution in um for low mag
end

%material properties: [lambda(um) real(n) imag(n)]
params.material = [1e4/3502.09657 1.47634 0.000258604];

%params.material = [6 1.47634 0.000258604];
params.spacing=1;   %spacing between point sources in the sampling radius

params.fov = round(params.res*pixRes);  %FOV in um

%params.fov = 30;


% this is ok for low mage mode, e.g. for resolution of 128 (at the detector) the field of
% view will be 141 um => sampling radius for interpolation will go from 0
% to 63; actually to cover all of FOV it should go up to 69
% However, this gets worse in high mag mode. In this case the FOV would be
% 704 um while sampling radius range is [0 63].
% Possible solution: increase the size of the detector image at far field,
% or don't resample at all at far field. Do the resampling at the end after
% interpolation is completely done.
params.rad = 1:round(params.fov/2-1);

params.s=rng;   %control random number generation
%load pmma.mat
%load([datDir '/measPMMA_wavRi.mat'])
%params.material=material;
%params.numWav = size(material,1);

params.numWav=1;    %number of wavelengths to simulate
params.ps = [0 0 0];    %position of the sphere in x,y,z

params.samples=400; %number of Monte Carlo samples
params.orderEf=50;  %number of orders for the focused/incident field

params.padding = 1; %specify padding
params.E0 = 1;  %amplitude of the field

% condenser NA pair of values specifying an inner obscuration
params.condenserNA_in = .34; 
params.condenserNA_out = 0.62;

% objective NA pair of values specifying an inner obscuration
params.objectiveNA_in = .34; 
params.objectiveNA_out = 0.62;

%this gives the kVec direction from the y axis
params.theta=1.5708;    
params.phi=0;

%% point sources
numRad = numel(params.rad);

%t = exp(-(([0:numRad].^2))/(2*(sigma^2)))/(2*pi*(sigma^2));
%Gaussian amplitudes for the specified sigma and number of sampling points
%in the radius of interpolation
amplitude = exp(-(([1:numRad].^2))/(2*(sigma^2)))/(2*pi*(sigma^2));
amplitude=amplitude./max(amplitude(:)); %normalize the amplitudes
%amplitude = changeRange(amplitude, 0, 0.46);

%show sampling of point sources
figure

rho = params.rad(end);  %maximum radius
%number of points to sample in the outer most radius -- chose such that the
%length of the sectors is 1 -- this remains constant for all other radii
numPoints = ceil(((2*pi*rho-params.spacing)/params.spacing + 1));
maxRange = floor(2*pi*rho-params.spacing/2);
s = linspace(0,maxRange,numPoints); %sectors
theta = s./rho; %angles corresponding to sectors
params.numPS = numPoints;

theta = 2*pi/numPoints:2*pi/numPoints:2*pi;
scatter(0,0,[],[1 1 0],'filled')
hold on
totalPS = 1;
if params.numPS>1
    
    %[pf_x, pf_z, pf_theta] = pointSources(params.rad, params.numPS);
    numPoints = round(2*pi*params.rad(end));
    
    radiusIdx=1;
    for i=1:numRad
        %theta=linspace(0,2*pi*i-0.1,numPoints);
        rho=params.rad(radiusIdx);
        numPoints = ceil(((2*pi*rho-params.spacing)/params.spacing + 1));
        theta = 2*pi/numPoints:2*pi/numPoints:2*pi;
        
        %numPoints = round(((2*pi*rho-params.spacing)/params.spacing + 1));
        %numPoints = round(2*pi*rho);
        %maxRange = floor(2*pi*rho-params.spacing/2);
        %s = linspace(0,maxRange,numPoints);
        %s = 0:numPoints;
        %theta = s./rho;
        [x,y] = pol2cart(theta,rho);
        scatter(x,y,[],[1 1 0]* amplitude(i),'filled')
        radiusIdx = radiusIdx + 1;
        totalPS = totalPS + numPoints;
    end
    axis square
    set(gca,'color','black')
end

params.amplitude = amplitude(1:params.spacing:end);

%% run forward BimSIm
% version -- interpolate at the detector for extended sources

output = gpuBimSimInterpAtDetector(params);

return

%% add noise to intensity and incident field images
% Add white Gaussian noise E to Ax, scaled such that
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