% function [C_final, oNoiseVar] = addONoise(C, SNR)
%
% Computes necessary noise variance for a given SNR, and adds it
%
% Author: Daniel Spendley, Daniel Rudoy
% Created:  03/06/2007
% Modified: 12/13/2007, 03/06/2007

function [C_final, oNoiseVar] = addONoise(C, SNR)

% Number of Observations
num_obser = length(C(1,:));

% Cepstral Order
cep_order = length(C(:,1));

% Determine Noise Power to Add for given SNR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct the Volume of the C's and average them over all observations
for i = 1:num_obser
    volume_C(i) = (sum(C(:,i).^2));
end

avg_vol_C = sum(volume_C)/num_obser;

% For the given SNR ratio, compute the neccessary average noise:
% 10*log10(avg_signal_power/noise_power) = SNR_dB
avg_noise_volume  = avg_vol_C/(10^(SNR/10));
oNoiseVar = avg_noise_volume/cep_order;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add noise back into the cepstrum
for i = 1:num_obser
    C_final(:,i) = C(:,i) + sqrt(oNoiseVar)*randn(cep_order,1);
end