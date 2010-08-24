function [C, oNoiseVar, C_clean] = genNoisyObser(SNR, frm, bw, cepOrder,fs)
% Generate noisy observations given true formant state data

% Get Desired lengths
[rows columns] = size(frm);

% Calculate Observations for state progression
for kk = 1:columns
    C(:,kk) = fb2cp(frm(:,kk),bw(:,kk),cepOrder,fs);
end

C_clean = C;

% Add observation noise
[C, oNoiseVar] = addONoise(C_clean, SNR);
