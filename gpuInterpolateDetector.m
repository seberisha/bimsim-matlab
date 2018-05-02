function [D_Et, D_Ef] = gpuInterpolateDetector(params)

R = numel(params.rad);
radius = params.rad(end);

D_x_ps = params.x(params.startIdx:params.endIdx,params.startIdx:params.endIdx) - params.ps(1);
D_z_ps = params.z(params.startIdx:params.endIdx,params.startIdx:params.endIdx) - params.ps(3);
[~, D_r_ps] = cart2pol(D_x_ps,D_z_ps);   %rho polar values at the above cartesian coordinates


totalEf_profile = zeros(1,R,'gpuArray');
totalEt_profile = totalEf_profile;

numPoints = 1;
rho = [params.rad 0];
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
    [E_t, E_f] = gpuSimulateScattering(params);
    
    % run the far-field simulation
    [D_Et_psi, D_Ef_psi] = gpuSimulateImaging(E_t, E_f, params);
    
    % detector x, z coordinates for Ef
    x_fp = params.x - params.pf(1);
    z_fp = params.z - params.pf(3);
    D_x_fp = x_fp(params.startIdx:params.endIdx,params.startIdx:params.endIdx);
    D_z_fp = z_fp(params.startIdx:params.endIdx,params.startIdx:params.endIdx);
    [~, D_r_fp] = cart2pol(D_x_fp,D_z_fp);   %rho polar values at the above cartesian coordinates
    
    
    
    params.pf_theta = 2*pi/numPoints:2*pi/numPoints:2*pi;
    
    temp_Ef = (simRing(numPoints, R, radius, D_Ef_psi, D_x_ps, D_z_ps));
    
    subplot(121)
    plot(temp_Ef)
    hold on
    drawnow
    
    totalEf_profile = totalEf_profile + temp_Ef;
    subplot(122)
    plot(totalEf_profile)
    hold on
    drawnow
    
    temp_Et = (simRing(numPoints, R, radius, D_Et_psi, D_x_ps, D_z_ps));
    
    totalEt_profile = totalEt_profile + temp_Et;
    
    numPoints = ceil(((2*pi*rho(i+1)-params.spacing)/params.spacing + 1));
    dispvar()
    
end

D_Ef = interp1([params.rad], totalEf_profile,D_r_ps);
D_Ef(isnan(D_Ef))=0; D_Ef(isinf(D_Ef))=0;
D_Et = interp1([ params.rad], totalEt_profile, D_r_ps);
D_Et(isnan(D_Et))=0; D_Et(isinf(D_Et))=0;