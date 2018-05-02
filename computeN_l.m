function N_l = computeN_l(a, lambda)
%   Description:    Compute the maximum order of the field.
%
%   Input:          -a - radius of the sphere
%                   -lambda - wavelength in micrometers
%
%   Output:         N_l - maximum order of the field
%
%   References:     C. F. Bohren and D. R. Huffman. Absorption and Scattering of Light
%                   by Small Particles. John Wiley & Sons, Sept. 2008.
%
%   Authors:        S. Berisha
%
%   Last modified:  11/24/15
N_l = ceil(2*pi.*a./lambda + 4*(2*pi.*a./lambda).^(1/3) + 2);

