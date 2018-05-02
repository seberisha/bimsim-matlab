%%
clear
% first derivative of the spherical Bessel function of the first kind
j = sym('sqrt(1/2*pi/x)*besselj(n+1/2,x)');
dj = simplify(diff(j));
dj = vectorize(inline(char(dj),'n','x'));

% first derivative of the spherical Hankel function of the first kind
h = sym('((pi/(2*x))^.5)*besselj(n+.5,x)+i*((pi/(2*x))^.5)*bessely(n+.5,x)');
dh = simplify(diff(h));
dh = vectorize(inline(char(dh),'n','x'));

%%
type = 1;
res=1;
order=1;
%r=linspace(0,1,res);
a=1;lambda=1;
k=2*pi/lambda; ka= k*a; n=1.4; kna=k*n*a;
k_vec = [1 0 0];
% [x,y,z] = sph2cart(0,0,r);

r=ones(res,1);
theta = linspace(0,2*pi,res)';
[x,y,z] = sph2cart(theta,0,r);
% x = r .* cos(zeros(size(r))) .* cos(theta);
% y = r .* cos(zeros(size(r))) .* sin(theta);
% z = r .* sin(zeros(size(r)));

for j=1:res
    rVecs(:,j)=[x(j);y(j);z(j)];
    rVecs(:,j)=rVecs(:,j)/norm(rVecs(:,j));
end

%r_vec=[x;y;z];

cos_theta=k_vec*rVecs;

j_ka = squeeze(sphbesselj(order,ka,'one'));
orderVec = (0:order)';

if type==1
    h_ka_p = dh((order),ka);
    j_ka_p = dj((order)',ka);
    h_ka = shank1((order)',ka,'one');
else
    h_ka_p = dh((0:order)',ka);
    j_ka_p = dj((0:order)',ka);
    h_ka = shank1((0:order)',ka,'one');
end

j_kna = squeeze(sphbesselj(order,kna,'one'));
if type==1
    j_kna_p = dh((order),kna);
else
    j_kna_p = dh((0:order)',kna);
end
num_a = j_ka.*h_ka_p - j_ka_p.*h_ka;
den = j_kna.*h_ka_p - h_ka.*j_kna_p*n;

if type==1
    A = (2*(order) + 1).*1i.^((order)).*(num_a./den);
else
    A = (2*(0:order)' + 1).*1i.^((0:order)').*(num_a./den);
end

knr=k*n*r;
j_knr = squeeze(sphbesselj(order,knr,'one'));

P_costh = myLegendre(order,cos_theta);
if type==1
    if isscalar(r)
        jknr_Pcosth = j_knr.*P_costh(end);
    else
        
        jknr_Pcosth = j_knr.*P_costh(:,end);
    end
else
    jknr_Pcosth = j_knr.*P_costh;
end
if isscalar(r)
    jknr_Pcosth = jknr_Pcosth.';
end
if type==1
    jknr_Pcosth =jknr_Pcosth*A;
else
    for j=0:order
        jknr_Pcosth(:,j+1)=jknr_Pcosth(:,j+1)*A(j+1);
    end
end
E_i = sum(jknr_Pcosth,2);

kr=k*r;
j_kr = squeeze(sphbesselj(order,kr,'one'));
if (type==1)
    if isscalar(r)
        jkr_Pcosth = j_kr.*P_costh(end);
    else
        jkr_Pcosth = j_kr.*P_costh(:,end);
    end
else
    jkr_Pcosth = j_kr.*P_costh;
end
if isscalar(r)
    jkr_Pcosth=jkr_Pcosth.';
end
if type==1
    jkr_Pcosth = jkr_Pcosth*(2*order+1)*1i^order;
else
    
    for j=0:order
        jkr_Pcosth(:,j+1) = jkr_Pcosth(:,j+1)*(2*j+1)*1i^j;
    end
end
E_f = sum(jkr_Pcosth,2);

num_b = j_ka.*j_kna_p.*n - j_kna.*j_ka_p;
if type==1
    B = (2*(order)+1).*1i.^(order).*(num_b./den);
else
    B = (2*(0:order)'+1).*1i.^(0:order)'.*(num_b./den);
end
h_kr = squeeze(shank1(order,kr,'one'));

if type==1
    if isscalar(r)
        hkr_Pcosth = h_kr.*P_costh(end);
    else
        
        hkr_Pcosth = h_kr.*P_costh(:,end);
    end
else
    hkr_Pcosth = h_kr.*P_costh;
end
if isscalar(r)
    hkr_Pcosth=hkr_Pcosth.';
end
if type==1
    hkr_Pcosth = hkr_Pcosth*B;
else
    
    for j=0:order
        hkr_Pcosth(:,j+1) = hkr_Pcosth(:,j+1)*B(j+1);
    end
end

E_s = sum(hkr_Pcosth,2);
E_t = E_f+E_s;
norm(E_i - (E_t))
figure, plot((abs((E_f)))),title('E_f')
figure, plot((abs((E_s)))),title('E_s')
figure, plot((abs((E_t)))),title('E_t E_i')
hold on
plot((abs((E_i))),'r')


