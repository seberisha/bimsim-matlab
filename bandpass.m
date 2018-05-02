function BPF = bandpass(fmag, NAin, NAout, lambda)

objMinParam = NAin/lambda; %min cut off of frequencies
objMaxParam = NAout/lambda; %max cut off of frequencies
BPF = (fmag < objMinParam | fmag > objMaxParam);   %compute the band pass filter