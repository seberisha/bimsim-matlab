%% setup
M=256;N=256;
lambda=3;
k = 2*pi/lambda;
n_0 = rand; NA_o=1;
f_c = (2*pi/lambda)*NA_o; S = 100; a=3; deltaU = 1/S;
res=256; lambda=3; N_l=100; alpha_1=0; NA=NA_o; E_0=1; a=3; n=3; c=0; gridIdx=40; 

%measured data
A_m = rand(M,N);
%imaginary part of the complex refractive index
k_j_v=rand;
%% algorithm
j=1;
scale = log(10)/(4*pi*v(2*R));

while(1)
    %calculate the predicted data A_j(r,v,n(v),k(v), v)
    simulateFullField
    %evaluate the difference E_j = A - A_j
    E_j = A_m - A;
    %update the imaginary part of the complex refractive index
    k_j_v = k_j_v + scale.*E;
    k_j(kv<0)=0;
    n = n_0 + KK(kv);
end