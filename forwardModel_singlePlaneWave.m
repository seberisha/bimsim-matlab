function [E_t, E_f, E_s, E_i] = forwardModel_singlePlaneWave(N_l, E_0, a, r, j_kr, P_ct, B, h_kr, A, j_knr, e_ikc)
%% setup

%params
%res = 256; lambda=3; N_l = 100; alpha_1=0; NA=.9; E_0=1; a=3; n=3;
%c=0; gridIdx=20;
% [x,y, z] = meshgrid(linspace(-gridIdx,gridIdx,res), linspace(-gridIdx,gridIdx,res), linspace(-gridIdx,gridIdx,res));
% [theta, r] = cart2pol(x,y,z);


%[theta, r] = cart2pol(x,y);

%% Compute E_f
E_f = Ef_singlePlaneWave(j_kr, P_ct, E_0, N_l);

%% 
%Compute E_s
E_s = computeEs(B, h_kr, P_ct, N_l, e_ikc);

%Compute E_i
E_i = computeEi(A, j_knr, P_ct, e_ikc, N_l);

E_i(r>a) = 0;
E_s(r<a) = E_i(r<a);
E_i(abs(r-a)<=1e-8) = E_f(r==a) + E_s(r==a);



E_t = E_s + E_f;





