function [H] = getH_FZ(freqVals, bwVals, numFormants, cepOrder, fs)

% INPUT
%  freqVals = [formantFrq zeroFreq]
%  bwVals  = [formantBW zeroBW]
%
% Return linear approximation of the nonlinear observation function h()
% around the conditional mean. This function assumes that only formants
% are random and bandwidths are being provided.

% Author: Daniel Rudoy
% Created: 02/12/10
% Modified: 02/15/10

H = zeros(cepOrder, length(freqVals));
for i = 1:cepOrder
    % This loop updates the matrix based on pole values
    for j = 1:numFormants
        bw = exp(-pi * i * bwVals(j)/fs);
        H(i,j) =  -4*pi/fs*bw*sin(2*pi*i*freqVals(j)/fs);
    end
    % This loop updates the value above based on zero values
    for j = numFormants + 1:length(freqVals)
        bw = exp(-pi * i * bwVals(j)/fs);
        H(i,j) =  4*pi/fs*bw*sin(2*pi*i*freqVals(j)/fs);
    end
end