function [E_s , E_i] = parComputeScatteredFields(a,psMask, n,rVecs,wavNum,r,numOrd,simRes,ps,pf,samples,k_j,subA,B,A)

%normalize the position vectors
normRvecs = bsxfun(@rdivide, rVecs,  sqrt(sum(rVecs.^2,2)));

kr = wavNum.*r;
hl_kr = shank1(numOrd, kr, 'multiple');

%The internal field specifying the field inside of a sphere
%for an incident plane wave

%compute the argument for the spherical bessel function k*n*r
knr = wavNum * n .* r;

%compute the spherical bessel function
jl_knr = sphbesselj(numOrd,knr,'multiple');

E_s = zeros(simRes,simRes);
E_i = E_s;
%h = waitbar(0, 'Monte Carlo integration...');


% r = sqrt(normKvec(1)^2 + normKvec(2)^2 + normKvec(3)^2);
% theta =  atan2(normKvec(2), normKvec(1));
% phi = acos(normKvec(2)/r);
% kSpherical = [r theta phi]

%the vector from the focal point to the center of the sphere
c = ps - pf;
display('Starting Monte Carlo simulation...')


parfor i=1:samples
    cos_theta = (k_j(:,i)'*normRvecs')';
    cos_theta = reshape(cos_theta, simRes, simRes);
    %compute the legendre polynomials needed for computing E_s and E_i
    Pl_costheta = myLegendre(numOrd,cos_theta);
    hlkr_Plcostheta = hl_kr.*Pl_costheta;
    
    %multiply by the legendre polynomial
    jlknr_Plcostheta = jl_knr.*Pl_costheta;
    
    %multiply by the scattering coefficients B
    %sum all orders
    phase = exp(1i.*wavNum.*k_j(:,i)'*c');
    E_s = E_s + (1/samples).*subA.*phase.*sum(bsxfun(@times, hlkr_Plcostheta, reshape(B,[1 1 numOrd+1])),3);
    %E_s(r<A) = 0;%E_i(r<a);
    %subplot(1,2,1), imagesc((abs((E_s)))),title('E_s'), colorbar, axis image, colormap(brewer)
    E_i = E_i + (1/samples).*subA.*phase.*sum(bsxfun(@times, jlknr_Plcostheta, reshape(A,[1 1 numOrd+1])),3);
    %E_i(r>A) = 0;
    %subplot(1,2,2), imagesc((abs((E_i)))),title('E_i'), colorbar, axis image, colormap(brewer)
    %pause(0.1)
    %waitbar(i/samples);
end
display('Finished Monte Carlo simulation...')

%close(h)

E_s(psMask<a) = 0;%E_i(r<a);

E_i(psMask>a) = 0;