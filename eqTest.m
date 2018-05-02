%first derivative of the spherical Bessel function of the first kind
j = sym('sqrt(1/2*pi/x)*besselj(n+1/2,x)');
dj = simplify(diff(j));
dj = vectorize(inline(char(dj),'n','x'));

first derivative of the spherical Hankel function of the first kind
h = sym('((pi/(2*x))^.5)*besselj(n+.5,x)+i*((pi/(2*x))^.5)*bessely(n+.5,x)');
dh = simplify(diff(h));
dh = vectorize(inline(char(dh),'n','x'));

%
order=1; r=1; a=1;lambda=1;
k=2*pi/lambda; ka= k*a; n=1.4; kna=k*n*a;
k_vec = [1 0 0]; 
[x,y,z] = sph2cart(0,0,r);
r_vec=[x;y;z];

cos_theta=k_vec*r_vec;

j_ka = sphbesselj(order,ka,'one');
h_ka_p = dh(order,ka);
j_ka_p = dj(order,ka);
h_ka = shank1(order,ka,'one');
j_kna = sphbesselj(order,kna,'one');
j_kna_p = dh(order,kna);
num_a = j_ka*h_ka_p - j_ka_p*h_ka;
den = j_kna*h_ka_p - h_ka*j_kna_p*n;

A = (2*order + 1)*1i^order*(num_a/den);

knr=k*n*r;
j_knr = sphbesselj(order,knr,'one');

P_costh = myLegendre(order,cos_theta);
E_i = A*j_knr*P_costh;

kr=k*r;
j_kr = sphbesselj(order,kr,'one');
E_f = (2*order+1)*1i^order*j_kr*P_costh;

num_b = j_ka*j_kna_p*n - j_kna*j_ka_p;
B = (2*order+1)*1i^order*(num_b/den);
h_kr = shank1(order,kr,'one');
E_s = B*h_kr*P_costh;
E_t = E_f+E_s;
norm(E_i - (E_t))
figure, plot((abs((E_f)))),title('E_f')
figure, plot((abs((E_s)))),title('E_s')
figure, plot((abs((E_t)))),title('E_t E_i')
hold on
plot((abs((E_i))),'r')


