function [y] = sinc(x)
% SINC evaluates t == 0 ? 1 : sin(pi*t)/(pi*t)
y = sin(pi * x) ./ (pi * x);
y(x==0) = 1;

end