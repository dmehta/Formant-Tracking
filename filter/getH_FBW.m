function [H] = getH_FBW(frmVals, bwVals, numFormants, cepOrder, fs)

% Return linear approximation of the nonlinear observation function h()
% around the conditional mean. This function assumes that both formants
% and bandwidths are being tracked.

% Author: Daniel Rudoy
% Created: 02/12/10
% Modified: 02/15/10

H = zeros(cepOrder, 2*numFormants);

for ii = 1:cepOrder
    for jj = 1:numFormants
        bw = exp(-pi * ii * bwVals(jj)/fs);
        H(ii,jj) =  -4*pi/fs*bw*sin(2*pi*ii*frmVals(jj)/fs);
    end

    for jj = numFormants + 1 : 2*numFormants
        bw = exp(-pi * ii * bwVals(jj-numFormants)/fs);
        H(ii,jj) =  -2*pi/fs*bw*cos(2*pi*ii*frmVals(jj-numFormants)/fs);
    end
end