function [psProfileVec_Et, psProfileVec_Ef,output] = gpuComputePointSourcesInterpAtDetector(numPoints, D_x, D_y, materialIdx,D_Et, D_Ef, E_t, params, output)

if(params.interpolate)
    
    %focal point for the first point source in the ring
    params.pf(1) = params.pf_x(1);
    params.pf(2) = 0;
    params.pf(3) = params.pf_z(1);
    
    %compute E_f, E_t for the first point source
    E_f = gpuComputeEf(params, params.kVecs(materialIdx,:));
    E_f(isnan(E_f))=0; E_f(isinf(E_f))=0;
    
    % scattered and internal fields for the first point source
    [E_s, E_i] = gpuComputeScatteredFields(params, params.material(materialIdx,:));
    
    % total field for the first point source
    E_t(params.psMask<params.a) = E_i(params.psMask<params.a);
    E_t(params.psMask>=params.a) = E_f(params.psMask>=params.a) + E_s(params.psMask>=params.a);
    E_t(isnan(E_t))=0; E_t(isinf(E_t))=0;
    
    %apply bandpass filtering in the frequency domain
    [Et_params.BPF, Ef_params.BPF] = applyBandPass(E_t, E_f, params.BPF);
    
    %resample to get far field
    [Et_ff, Ef_ff] = getFarField(Et_params.BPF, Ef_params.BPF,params.startIdx, params.endIdx);
    
    %save intensities at the detector for the first point source
    D_Et_ps1 = (abs(Et_ff)).^2;
    D_Ef_ps1 = (abs(Ef_ff)).^2;
    
    psProfileVec_Ef = zeros(1,numel(params.rad)+1,'gpuArray');
    psProfileVec_Et = psProfileVec_Ef;
    psProfileVec_Ef(1) = interp2(D_x, D_y, D_Ef_ps1,0,0);
    psProfileVec_Et(1) = interp2(D_x, D_y, D_Et_ps1,0,0);
    numRad = numel(params.rad);
    % interpolate for all rings
    for ring=1:numRad
        
        %params.scale = 1/numPoints;
        params.scale=1;
        rho=params.rad(ring);  %radius
        %x and z coordinates for the point source ring
        [params.pf_x,params.pf_z] = pol2cart(params.pf_theta,rho);
        
        D_Ef_r = interp2(D_x, D_y, D_Ef_ps1,params.pf_x,params.pf_z);
        psProfileVec_Ef(ring+1) = params.scale*sum(D_Ef_r);
        D_Et_r = interp2(D_x, D_y, D_Et_ps1, params.pf_x, params.pf_z);
        psProfileVec_Et(ring+1) = params.scale*sum(D_Et_r);
    end
    
%     % averaged and scaled value of Ef and Et for each ring
%     output.R_Ef(rIdx+1) = sum(D_Ef_r)*params.scale;
%     output.R_Et(rIdx+1) = sum(D_Et_r)*params.scale;
else
    numPS = numel(params.pf_x);
    %all point sources
    for p=1:numPS
        
        pf(1) = params.pf_x(p);
        pf(2) = 0;
        pf(3) = params.pf_z(p);
        
        E_f = gpuComputeEf(params);
        E_f(isnan(E_f))=0; E_f(isinf(E_f))=0;
        
        [E_s, E_i] = gpuComputeScatteredFields(params, params.material(materialIdx,:));
        E_t(params.psMask<params.a) = E_i(params.psMask<params.a);
        E_t(params.psMask>=params.a) = E_f(params.psMask>=params.a) + E_s(params.psMask>=params.a);
        
        E_t(isnan(E_t))=0; E_t(isinf(E_t))=0;
        
        [Et_params.BPF, Ef_params.BPF] = applyBandPass(E_t, E_f, params.BPF);
        
        [Et_ff, Ef_ff] = getFarField(Et_params.BPF, Ef_params.BPF,params.startIdx, params.endIdx);
        
        %integrate and resample
        D_Et = D_Et + (abs(Et_ff)).^2;   %this is out_i in cuda bimsim -- the measured intesity
        D_Ef = D_Ef + (abs(Ef_ff)).^2;   %this out_inc in cuda bimsim -- the measured incident field
        
    end
    
end