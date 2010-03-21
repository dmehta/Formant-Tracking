% CepTest.m - Cepstral order test
clear all

% Set fixed parameters
snr_dB = 15; % input SNR (dB)
fs = 16000; % sample rate
num_formants = 6; % Number of formants to generate
num_MC_runs = 25; % Monte Carlo runs
cepstral_index = 15:5:40; % Cepstral order; must be at least twice num_formants...

for i = 1:length(cepstral_index),
    for j = 1:num_MC_runs;
        if j==1, % New random seed
            label{i} = num2str(cepstral_index(i)); % for plot legend
        end
        E = TrackExp('Synth',snr_dB,cepstral_index(i),fs,num_formants);
        avg_rmse(i,j,:) = mean(E.rmse,2);
    end
end

mc_avg_avg_rmse = squeeze(mean(avg_rmse,2));
%mc_var_avg_rmse = squeeze(var(avg_rmse,2));
figure(3), clf 
plot(mc_avg_avg_rmse','o-')
xlim([0.5 num_formants+0.5])
xlabel('Resonance Number')
ylabel('Root MSE (Hz)')
title(['Average Root MSE for ' int2str(num_MC_runs) ' Trials of various cepstral orders'])
grid
legend(label)
