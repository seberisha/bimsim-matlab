function [Et_bpf, Ef_bpf] = applyBandPass(E_t, E_f, BPF)


Et_d = (fft2(E_t));
fftEf = fft2(E_f);
%bandpass
Et_d(BPF) = 0;
fftEf(BPF) = 0;

Et_bpf = ifft2((Et_d));
Ef_bpf = ifft2(fftEf);