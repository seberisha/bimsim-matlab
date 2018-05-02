%%
clear

%% setup
NA_in = 0;
NA_out= 1;
E_0=1;
res=256;
a=1;lambda=1;
n=1.4;
k_vec = [1 0 0];

gridIdx=15;

alpha_1=asin(NA_in);


alpha_2=asin(NA_out);

order= computeN_l(a, lambda);

k=2*pi/lambda; ka= k*a;  kna=k*n*a;
[x,y] = meshgrid(linspace(-gridIdx,gridIdx,res), linspace(-gridIdx,gridIdx,res));
z=zeros(size(x));
[theta, phi, r] = cart2sph(x, y, z);

% r=ones(res,res);
% theta = meshgrid(linspace(0,2*pi,res),linspace(0,2*pi,res));
% [x,y,z] = sph2cart(theta,0,r);

idx=1;
rVecs=zeros(3,res*res);
for j=1:res
    for i=1:res
        temp = [x(i,j); y(i,j); z(i,j)];
        rVecs(:,idx)=temp./norm(temp);
        idx=idx+1;
    end
end

cos_theta=k_vec*rVecs;
cos_theta = reshape(cos_theta, res, res);

% first derivative of the spherical Bessel function of the first kind
j = sym('sqrt(1/2*pi/x)*besselj(n+1/2,x)');
dj = simplify(diff(j));
dj = vectorize(inline(char(dj),'n','x'));

% first derivative of the spherical Hankel function of the first kind
h = sym('((pi/(2*x))^.5)*besselj(n+.5,x)+i*((pi/(2*x))^.5)*bessely(n+.5,x)');
dh = simplify(diff(h));
dh = vectorize(inline(char(dh),'n','x'));


j_ka = squeeze(sphbesselj(order,ka,'multiple'));
orderVec = (0:order)';


%h_ka_p = dh((0:order)',ka);
h_ka_p = derivSphHan(order, ka);

%j_ka_p = dj((0:order)',ka);

j_ka_p = derivSphBes(order, ka);
h_ka = squeeze(shank1(order,ka,'multiple'));

j_kna = squeeze(sphbesselj(order,kna,'multiple'));
%j_kna_p = dj((0:order)',kna);

j_kna_p = derivSphBes(order, kna);


num_a = j_ka.*h_ka_p - j_ka_p.*h_ka;
den = j_kna.*h_ka_p - h_ka.*j_kna_p*n;

A = (2*(0:order)' + 1).*(1i.^((0:order)')).*(num_a./den);


knr=k*n*r;
j_knr = squeeze(sphbesselj(order,knr,'multiple'));



kr=k*r;


num_b = j_ka.*j_kna_p.*n - j_kna.*j_ka_p;
B = (2*(0:order)'+1).*(1i.^((0:order)')).*(num_b./den);
h_kr = squeeze(shank1(order,kr,'multiple'));


%% compute E_f
order=100;
j_kr = squeeze(sphbesselj(order,kr,'multiple'));

P_costh = myLegendre(order,cos_theta);

jkr_Pcosth = j_kr.*P_costh;
    
for j=0:order
    jkr_Pcosth(:,:,j+1) = jkr_Pcosth(:,:,j+1)*(2*j+1)*1i^j;
end

E_f = sum(jkr_Pcosth,3);

brewer = brewermap(1000);
figure, imagesc((abs((E_f)))),title('E_f'), colorbar, axis image
colormap(brewer)

%% compute E_s and E_i via Monte Carlo integration
samples=1;
NA_in=0; NA_out=1;
k_j = monteCarlo(samples,k_vec, NA_in, NA_out);
E_s = zeros(res,res);
E_i = E_s;
h = waitbar(0, 'Monte Carlo integration...');
%compute the amplitude that makes it through the condenser
amplitude=1;
subA = 2 * pi * amplitude * ( (1 - cos(asin(alpha_2))) - (1 - cos(asin(alpha_1))) );
subA=1;
for i=1:samples
    kRvecs=(k_j(:,i)'*rVecs);
    %kRvecs=(k_vec*rVecs);
    cos_theta=reshape(kRvecs, res,res);
    % phase shift
    % c = p_s - p_f;
    % e_ikc = exp(1i*k_j(:,i)'*c);
    
    % Legendre polynomials
    P_ct = myLegendre(order,cos_theta);
    
    E_s = E_s + subA.*computeEs(B, h_kr, P_ct, order, 1);
    
    %Compute E_i
    E_i = E_i + subA.*computeEi(A, j_knr, P_ct, order, 1);
    
    E_i(r>a) = 0;
    E_s(r<a) = 0;%E_i(r<a);
    waitbar(i/(samples+1), h)
end
E_s = (1/samples).*E_s;
E_i = (1/samples).*E_i;
% E_s = 2*pi*(1-cos(alpha_2))*E_s;
% %
% E_i = 2*pi*(1-cos(alpha_2))*E_i;
close(h)

figure
subplot(1,2,1), imagesc((abs(E_s))),title('E_s'), colorbar, axis image
subplot(1,2,2), imagesc((abs(E_i))),title('E_i'), colorbar, axis image
colormap(brewer)
E_t = zeros(res,res);
E_t(r>=a) = E_s(r>=a) + E_f(r>=a);
E_t(r<a)=E_i(r<a);
figure, imagesc((abs(E_t))),title('E_t'), colorbar, axis image
colormap(brewer)


