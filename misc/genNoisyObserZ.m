function [C, oNoiseVar, C_clean] = genNoisyObserZ(SNR, f_frqs, f_bws, z_frqs, z_bws, cepOrder,fs)

% Generate observation sequence (LPCC) from a set of formants and zeros
% along with their respective bandwidths
%
% Author: Daniel Rudoy, Daryush Mehta
% Modified: 03/21/2010, 08/23/2010 (added C_clean), 08/31/10 (handled
% no-zero or no-pole cases)

% Calculate Observations for state progression

numObs = max(size(f_frqs, 2), size(z_frqs, 2));
C = zeros(cepOrder, numObs);
C1 = zeros(cepOrder, numObs);
C2 = zeros(cepOrder, numObs);

for kk = 1:numObs
    if ~isempty(f_frqs)
        C1(:,kk) = fb2cp(f_frqs(:,kk),f_bws(:,kk),cepOrder,fs);   % LPCC coeffs due to poles
    else
        C1(:,kk) = zeros(cepOrder,1);
    end
    
    if ~isempty(z_frqs)
        C2(:,kk) = fb2cp(z_frqs(:,kk),z_bws(:,kk),cepOrder,fs)';  % coeffs due to zeros
    else
        C2(:,kk) = zeros(cepOrder,1);
    end
    
    C(:,kk) = C1(:,kk) - C2(:,kk);  % Adjust due to effect of zeros
end

C_clean = C;

% Add observation noise
[C, oNoiseVar] = addONoise(C_clean, SNR);