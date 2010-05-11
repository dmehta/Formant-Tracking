function H = getH_FZBW(freqVals, bwVals, numFormants, cepOrder, fs)

% Return linear approximation of the nonlinear observation function h()
% around the conditional mean. This function assumes that both formants
% and bandwidths are random and thus tracked.
% 
% INPUT
%  freqVals = [formantFrq zeroFreq]
%  bwVals  = [formantBW zeroBW]

% Author: Daniel Rudoy, Daryush Mehta
% Created: 02/12/10
% Modified: 05/08/10 (adding random bandwidths)

H = zeros(cepOrder, length(freqVals)+length(bwVals));
numAntiformants = length(freqVals)-numFormants;

for ii = 1:cepOrder
    % Thiis loop updates the matrix based on resonance center frequencies
    for jj = 1:numFormants
        bw = exp(-pi * ii * bwVals(jj)/fs);
        H(ii,jj) =  -4*pi/fs*bw*sin(2*pi*ii*freqVals(jj)/fs);
    end
    % This loop updates the matrix based on resonance bandwidths
    for jj = numFormants + 1:2*numFormants
        bw = exp(-pi * ii * bwVals(jj-numFormants)/fs);
        H(ii,jj) =  -2*pi/fs*bw*cos(2*pi*ii*freqVals(jj-numFormants)/fs);
    end
    % This loop updates the value above based on anti-resonance center
    % frequencies
    for jj = 2*numFormants + 1:2*numFormants + numAntiformants
        bw = exp(-pi * ii * bwVals(jj-numFormants)/fs);
        H(ii,jj) =  4*pi/fs*bw*sin(2*pi*ii*freqVals(jj-numFormants)/fs);
    end
    % This loop updates the value above based on anti-resonance bandwidths
    for jj = 2*numFormants + numAntiformants + 1:length(freqVals)
        bw = exp(-pi * ii * bwVals(jj-2*numFormants)/fs);
        H(ii,jj) =  2*pi/fs*bw*cos(2*pi*ii*freqVals(jj-2*numFormants)/fs);
    end    
end