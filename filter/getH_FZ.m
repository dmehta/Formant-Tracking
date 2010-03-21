function [H] = getH_FZ(frmVals, bwVals, numFormants, cepOrder, fs)

% INPUT
%  frmVals = [formantFrq zeroFreq]
%  bwVals  = [formantBW zeroBW]
%
% Return linear approximation of the nonlinear observation function h()
% around the conditional mean. This function assumes that only formants
% are random and bandwidths are being provided.

% Author: Daniel Rudoy
% Created: 02/12/10
% Modified: 02/15/10

H = zeros(cepOrder, length(frmVals));
for i = 1:cepOrder
    for j = 1:numFormants
        bw = exp(-pi * i * bwVals(j)/fs);
        H(i,j) =  -4*pi/fs*bw*sin(2*pi*i*frmVals(j)/fs);
    end
    for j = numFormants + 1:length(frmVals)
        bw = exp(-pi * i * bwVals(j)/fs);
        H(i,j) =  4*pi/fs*bw*sin(2*pi*i*frmVals(j)/fs);
    end
end