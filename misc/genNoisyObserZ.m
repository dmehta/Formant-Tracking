function [C, oNoiseVar] = genNoisyObserZ(SNR, f_frqs, f_bws, z_frqs, z_bws, cepOrder,fs)

% Generate observation sequence (LPCC) from a set of formants and zeros
% along with their respective bandwidths
%
% Author: Daniel Rudoy
% Modified: 03/21/2010

% Calculate Observations for state progression
for kk = 1:length(f_frqs)
    C(:,kk) = fb2cp(f_frqs(:,kk),f_bws(:,kk),cepOrder,fs);             % LPCC coeffs due to poles
    C(:,kk) = C(:,kk) - fb2cp(z_frqs(:,kk),z_bws(:,kk),cepOrder,fs)';  % Adjust due to effect of zeros
end

% Add observation noise
[C, oNoiseVar] = addONoise(C, SNR);
