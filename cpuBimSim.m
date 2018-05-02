%% set up and input parameters
%clear
addpath(genpath('~/source/stim-matlab/'));
datDir = '~/data/bimsim/matlab/bulkPMMA_6.5'
addpath(genpath('~/source/stim-matlab/'))

input.interpolate=2;   %options to interpolate for extended sources

input.constantNumSamples = 1;
input.fieldInsideSphere = 0;

input.resolution = 128;   %resolution at the detector

% standard deviation for a Gaussian extended source --to be estimated
% from FTIR interferogram data - for Agilent FTIR select show raw data in
% lancer control before start of data acquisition (when save dialog box pops up)
input.sigma = 11;

input.sphereRadius = 10;    %radius of the sphere
input.highMag = 0;    %for high mag mode simulation choose 1, 0 for low mag

%material properties: [lambda(um) real(n) imag(n)]
input.material = [1e4/3502.09657 1.47634 0.000258604];
input.fov = [];
input.spacing=1;   %spacing between point sources in the sampling radius

input.numWav=1;    %number of wavelengths to simulate
input.ps = [0 0 0];    %position of the sphere in x,y,z

input.samples=1000; %number of Monte Carlo samples
input.orderEf=100;  %number of orders for the focused/incident field

input.padding = 1; %specify padding
input.E0 = 1;  %amplitude of the field

% condenser NA pair of values specifying an inner obscuration
input.condenserNA_in = .34;
input.condenserNA_out = 0.62;

% objective NA pair of values specifying an inner obscuration
input.objectiveNA_in = .34;
input.objectiveNA_out = 0.62;

%this gives the kVec direction from the y axis
input.theta=1.5708;
input.phi=0;

%focal point
input.pf = [0 0 0];

% create a struct with input parameters for simulation
params = cpuBimSimParamsStruct(input);


%% run forward BimSim
%figure
if (input.interpolate == 1)
   tic 
   [D_Et, D_Ef] = cpuInterpolateDetector(params);
   dispvar('time in secs: ', toc)
elseif (input.interpolate == 2)
   tic
   [D_Et, D_Ef] = cpuInterpolateWholePointSourcesAtDetector(params);   
   dispvar('time in secs: ', toc)
else
    tic
    % run the near-field simulation
    [E_t, E_f] = cpuSimulateScattering(params);
    
    % run the far-field simulation
    [D_Et, D_Ef] = cpuSimulateImaging(E_t, E_f, params);
    dispvar('time in secs: ', toc)
end

%% save output
output = params;
output.A = -log10((D_Et)./(D_Ef));
output.absSpec = -log10(sum(D_Et(:))/sum(D_Ef(:)));
output.D_Et = D_Et;
output.D_Ef = D_Ef;
