%%
numPs = 200;
myFPs = zeros(numPs,3);
r = 6.5;
for i=1:numPs
    myFPs(i,:) = -r + rand(1,3)*2*r;
    myFPs(i,2)= params.a/2;
end
