function filter3 = gaussianbpf(nx,d0,d1,dist)
% Butterworth Bandpass Filter
% This simple  function was written for my Digital Image Processing course
% at Eastern Mediterranean University taught by
% Assoc. Prof. Dr. Hasan Demirel
% for the 2010-2011 Spring Semester
% for the complete report:
% http://www.scribd.com/doc/51981950/HW4-Frequency-Domain-Bandpass-Filtering
%
% Written By:
% Leonardo O. Iheme (leonardo.iheme@cc.emu.edu.tr)
% 24th of March 2011
%
%   I = The input grey scale image
%   d0 = Lower cut off frequency
%   d1 = Higher cut off frequency
%
% The function makes use of the simple principle that a bandpass filter
% can be obtained by multiplying a lowpass filter with a highpass filter
% where the lowpass filter has a higher cut off frquency than the high pass filter.
%
% Usage GAUSSIANBPF(I,DO,D1)
% Example
% ima = imread('grass.jpg');
% ima = rgb2gray(ima);
% filtered_image = gaussianbpf(ima,30,120);
% Gaussian Bandpass Filter



% Initialize filter.
filter1 = ones(nx,1);
filter2 = ones(nx,1);
filter3 = ones(nx,1);

for i = 1:nx
         %dist = ((i-round(nx/2))^2)^.5;
         %dist=dist.*fs;
        % Use Gaussian filter.
        filter1(i) = exp(-dist(i)^2/(2*d1^2));%./( sqrt((2*pi))*d1 );
        filter2(i) = exp(-dist(i)^2/(2*d0^2));%./( sqrt((2*pi))*d0 );
        filter3(i) = 1.0 - filter2(i);
        filter3(i) = filter1(i).*filter3(i);
end
%plot(filter3)
%title('Frequency Domain Filter Function Image')