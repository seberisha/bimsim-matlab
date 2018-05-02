function gbpf = gaussianBPF(data,lowerCutOff,upperCutOff)
% Construct a gaussian bandpass filter given a vector of frequency magnitudes
% and the lower and upper cut off frequencies.

l = size(data,2);

switch l
case 1
  x = data(:,1);
  y = 0;
case 2
  x = data(:,1);
  y = data(:,2);
otherwise
  error('illegal data dimension')
end

numerator = (x.^2 + y.^2);
gauss1 = exp( -numerator / (2*upperCutOff^2) ) ;%/ ( sqrt((2*pi)^l)*upperCutOff^l );
%fg1 = fft((gauss1));
gauss2 = exp( -numerator / (2*lowerCutOff^2) );% / ( sqrt((2*pi)^l)*lowerCutOff^l );
%fg2 = fft((gauss2));
gbpf = gauss1.* (1 - gauss2);
%fgbpf = fg1.*(1-fg2);