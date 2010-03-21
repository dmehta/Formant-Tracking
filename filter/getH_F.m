function [H] = getH_F(frmVals, bwVals, numFormants, cepOrder, fs)

% Return linear approximation of the nonlinear observation function h()
% around the conditional mean. This function assumes that only formants
% are random and bandwidths are being provided.

% Author: Daniel Rudoy
% Created: 02/12/10
% Modified: 02/15/10

H = zeros(cepOrder, numFormants);
for ii = 1:cepOrder
    for jj = 1:numFormants
        bw = exp(-pi * ii * bwVals(jj)/fs);
        H(ii,jj) =  -4*pi/fs*bw*sin(2*pi*ii*frmVals(jj)/fs);
    end
end