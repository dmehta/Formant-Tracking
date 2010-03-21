%
%
% Generate observation sequence (LPCC) from a set of formants and zeros
% along with their respective bandwidths
%


function [C, oNoiseVar] = genNoisyObserZ(SNR, f_frqs, f_bws, z_frqs, z_bws, cepOrder,fs)
% Generate noisy observations given true formant state data



% Calculate Observations for state progression
for kk = 1:length(f_frqs)
    C(:,kk) = fb2cp(f_frqs(:,kk),f_bws(:,kk),cepOrder,fs);
    C(:,kk) = C(:,kk) - fb2cp(z_frqs(:,kk),z_bws(:,kk),cepOrder,fs)';
end

% Add observation noise
[C, oNoiseVar] = addONoise(C, SNR);
