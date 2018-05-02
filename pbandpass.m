function BPF = pbandpass(df, simRes, NAin, NAout, lambda)

[iv, iu] = meshgrid(0:simRes-1, 0:simRes-1);
u=zeros(size(iu)); v = u;
idx = find(iu <= simRes/2);
u(idx) = iu(idx);
idx = find(iu > simRes/2);
u(idx) = (iu(idx) - simRes+1);
u=u.*df;
idx = find(iv <= simRes/2);
v(idx) = iv(idx);
idx = find(iv > simRes/2);
v(idx) = iv(idx) - simRes + 1;
v=v.*df;
fmag = sqrt(u.*u + v.*v);    %compute the magnitude of the frequencies
objMinParam = NAin/lambda; %min cut off of frequencies
objMaxParam = NAout/lambda; %max cut off of frequencies
BPF = (fmag < objMinParam | fmag > objMaxParam);   %compute the band pass filter