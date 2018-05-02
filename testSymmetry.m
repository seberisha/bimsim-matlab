%%
r = params.r;
t=r(193,193:end);
jt = sphbesselj(0,params.wavNum.*params.r,'multiple');
tt = sphbesselj(0,params.wavNum.*t,'multiple');

jm = zeros(384,384);
for i=1:192
    temp = find(abs(r-t(i))<eps);
    r(temp)
    jm(temp) = tt(i);
    imagesc(jm),axis image
    pause(1e-7)

end

%%
r = params.r;
t=r(193,193:end);

rm = zeros(384,384);
for i=1:192
    temp = find(abs(r-t(i))< .1);
    %rm(temp)
    rm(temp) = t(i);
    imagesc(rm),axis image
    pause(1e-7)

end

%%
r = params.r;
rm = zeros(384,384);
C = unique(r);
%[tf, loc] = ismember(C, r);
loc = find(ismember(r,C));
rm(loc) = r(loc);
imagesc(rm),axis image

%%
r = params.r;
rm = zeros(384,384);

C = unique(r);
tic
jt = sphbesselj(0,params.wavNum.*C,'multiple');

for i=1:numel(C)
    loc = find(ismember(r,C(i)));
    jm(loc) = jt(i);
    %imagesc(jm),axis image
   % pause(1e-7)
end
toc

%[tf, loc] = ismember(C, r);
% loc = find(ismember(r,C));
% jm(loc) = jt(loc);
 imagesc(jm),axis image

%%
tic
jt = sphbesselj(0,params.wavNum.*r,'multiple');
toc

%%

r = params.r;
rm = zeros(384,384);

C = unique(r);

for i=1:numel(C)
    loc = find(ismember(r,C(i)))
    rm(loc) = r(loc);
    imagesc(rm),axis image
    pause(1e-7)
end

