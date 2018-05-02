%%
clear
test_1 = rand(768,768,50,'gpuArray') + rand(768,768,50,'gpuArray')*1i;
test_2 = rand(768,768  ,50, 'gpuArray') + rand(768,768,50,'gpuArray')*1i;

for i=1:400
    testA  = test_1.*test_2;
end

test_3 = rand(768,768,50) + rand(768,768,50)*1i;
test_4 = rand(768,768  ,50) + rand(768,768,50)*1i;

for i=1:400
    testB  = test_3.*test_4;
end