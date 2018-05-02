EsEf(r>=a)=E_s(r>=a) + E_f(r>=a);
EsEf=zeros(res,res);
EsEf(r>=a)=E_s(r>=a) + E_f(r>=a);
figure,imagesc(abs(EsEf))
colormap(brewer), axis image, colorbar
mEsEf=zeros(res,res);
mEsEf(r>=a)=outN(r>=a);
figure,imagesc(abs(mEsEf))
colormap(brewer), axis image, colorbar
norm(mEsEf(:)-EsEf(:))

%%



%%
y = sym('((pi/(2*x))^.5)*bessely(n+.5,x)');
dy = simplify(diff(y));
dy = vectorize(inline(char(dy),'n','x'));

deriv_y_ka = dy((0:N_l)', ka);


%%
hp = deriv_j_ka + 1i*deriv_y_ka;

%%
%deriv of spherical bessel j
N_l=order
js_nu_m_1 = squeeze(sphbesselj((0:N_l)'-1,ka,'one'));

js_nu_p_1 = squeeze(sphbesselj((0:N_l)'+1,ka,'one'));

d_j_s = 1/2*(js_nu_m_1 - (j_ka + ka*js_nu_p_1)/ka);