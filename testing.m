% //default plane size in microns
% s = 40;
% pos = 0;
% 
% 
% //calculate the plane corners and normal based on the size and position
% pMin = bsPoint(-s/2, -s/2, pos);
% pMax = bsPoint(s/2, s/2, pos);
% normal = bsVector(0, 0, 1);

pMin=[-5; 0; -5];
pMax=[5; 0; 5];
normal=[0; 1; 0];
a=pMin;b=pMax;c=normal;

A = a;		
Y = b - a;
X = c - a - Y;

%%
res=256;
rVecs=zeros(3,res*res);
idx=1;
[u,v]=meshgrid(linspace(-20,20,res),linspace(-20,20,res));
for j=1:res
    for i=1:res
        su=u(i,j)/res;
        sv=v(i,j)/res;
        test(:,idx)=[su;sv;0];
        rVecs(:,idx) = A + X * su + Y * sv;
        idx=idx+1;
    end
end

%%
d = sqrt(sum(rVecs.^2,1));
d = reshape(d,256,256);
td = sqrt(sum(test.^2,1));
td = reshape(td,256,256);
