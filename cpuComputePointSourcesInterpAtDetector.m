function [D_Et, D_Ef,output] = cpuComputePointSourcesInterpAtDetector(D_x, D_y, rIdx, materialIdx,D_Et, D_Ef, E_t, params, output)

if(params.interpolate)
    
    params.pf(2) = 0;
    params.pf(1) = params.pf_x(1);
    params.pf(3) = params.pf_z(1);
    
    %compute E_f, E_t for one point source in the ring
    
    E_f = cpuComputeEf(params, params.kVecs(materialIdx,:));
    E_f(isnan(E_f))=0; E_f(isinf(E_f))=0;
    
    [E_s, E_i] = cpuComputeScatteredFields(params, params.material(materialIdx,:));
    E_t(params.psMask<params.a) = E_i(params.psMask<params.a);
    E_t(params.psMask>=params.a) = E_f(params.psMask>=params.a) + E_s(params.psMask>=params.a);
    E_t(isnan(E_t))=0; E_t(isinf(E_t))=0;
    
    %apply bandpass fparams.iltering
    [Et_params.BPF, Ef_params.BPF] = applyBandPass(E_t, E_f, params.BPF);
    
    %resample to get far field
    [Et_ff, Ef_ff] = getFarField(Et_params.BPF, Ef_params.BPF,params.startIdx, params.endIdx);
    
    %save intensities for the first point source
    D_Et_ps1 = (abs(Et_ff)).^2;
    D_Ef_ps1 = (abs(Ef_ff)).^2;
    %[xq, yq] = pol2cart(params.pf_theta,ones(size(params.pf_theta)).*params.rad(rIdx));
    % interpolate for the rest of the point sources in the ring
    D_Ef_r = interp2(D_x, D_y, D_Ef_ps1,params.pf_x,params.pf_z);
    D_Et_r = interp2(D_x, D_y, D_Et_ps1,params.pf_x,params.pf_z);
    
    output.R_Ef(rIdx+1) = sum(D_Ef_r)*params.scale;
    output.R_Et(rIdx+1) = sum(D_Et_r)*params.scale;
    
else
    numPS = numel(params.pf_x);
    %all point sources
    for p=1:numPS
        
        pf(1) = params.pf_x(p);
        pf(2) = 0;
        pf(3) = params.pf_z(p);
        
        E_f = cpuComputeEf(params);
        E_f(isnan(E_f))=0; E_f(isinf(E_f))=0;
        
        [E_s, E_i] = cpuComputeScatteredFields(params, params.material(materialIdx,:));
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