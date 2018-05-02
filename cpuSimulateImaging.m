function [D_Et, D_Ef] = cpuSimulateImaging(E_t, E_f, params)

%angular spectrum and band pass filtering
[Et_bpf, Ef_bpf] = applyBandPass(E_t, E_f, params.BPF);

%resample to get far field
[Et_ff, Ef_ff] = getFarField(Et_bpf, Ef_bpf,params.startIdx, params.endIdx);

[D_Et, D_Ef] = cpuIntegrateDetector(params, Et_ff, Ef_ff);