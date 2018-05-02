function [D_Et, D_Ef] = cpuInterpolateDetector(params)
% This version computes the first point source at a certain distance from
% the location of the sphere. Then it approximates the rest of the point
% sources in the same ring by interpolating a 1D profile (radius) of the
% first point source at other focal points in all the rings (computing the
% contribution of a point source in the entire FOV).

R = numel(params.rad);
radius = params.rad(end);

D_x_ps = params.x;
D_z_ps = params.z;
[~, D_r_ps] = cart2pol(D_x_ps,D_z_ps);   %rho polar values at the above cartesian coordinates

totalEf_profile = 0;
totalEt_profile = totalEf_profile;

% center of fov in simulated resolution
uc = round(params.simRes/2) + 1;
vc = round(params.simRes/2) + 1;


figure
[X, Y] = meshgrid(1:params.simRes,1:params.simRes);


% compute and sum the 1d profiles of all point sources
for i=1:params.numRad
    
    dispvar('ring',i)
    %focal point for the first point source in the ring
    params.pf(1) = params.rad(i);
    params.pf(2) = 0;
    params.pf(3) = 0;
    
    params.E0 = 1;
    
    params.subA = 2 * pi * params.E0 * ( (1 - cos(params.alpha2)) - (1 - cos(params.alpha1)) );
    
    % run the near-field simulation
    [E_t, E_f] = cpuSimulateScattering(params);
    
    % backpropagation
    greenFun = exp(-1i*params.magKVector*params.sphereRadius);
    E_t = E_t.* greenFun;
    E_f = E_f .* greenFun;
    
    % compute far field later to avoid NaN values in interpolation since
    % point sources are sampled in the radius of FOV which is larger than
    % the radius of the resolution of the detector
    
    %angular spectrum and band pass filtering
    [Et_bpf, Ef_bpf] = applyBandPass(E_t, E_f, params.BPF);
    
    [D_Et_psi, D_Ef_psi] = cpuIntegrateDetector(params, Et_bpf, Ef_bpf);   
    
    
    temp_Ef = (simRing(i, uc, vc, params.numPoints, R, radius, D_Ef_psi, X, Y));

%     subplot(121)
%     plot(temp_Ef)
%     hold on
%     drawnow
    
    totalEf_profile = totalEf_profile + temp_Ef;
%     subplot(122)
%     plot(totalEf_profile)
%     hold on
%     drawnow
    
    temp_Et = (simRing(i, uc, vc, params.numPoints, R, radius, D_Et_psi, X, Y));
    
    totalEt_profile = totalEt_profile + temp_Et;
    
    %numPoints = ceil(((2*pi*rho(i+1)-params.spacing)/params.spacing + 1));
    dispvar()
    
end

D_Ef = interp1([params.rad], totalEf_profile,D_r_ps);
D_Ef(isnan(D_Ef))=0; D_Ef(isinf(D_Ef))=0;
D_Et = interp1([ params.rad], totalEt_profile, D_r_ps);
D_Et(isnan(D_Et))=0; D_Et(isinf(D_Et))=0;

%resample to get far field
[D_Et, D_Ef] = getFarField(D_Et, D_Ef,params.startIdx, params.endIdx);
