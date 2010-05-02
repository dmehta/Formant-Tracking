%% RMSE vs number of cepstral coefficients
clear
close all

%% parameters
F = [500 1500 2500 3500]';
Fbw = [100 100 100 100]';
Z = [700 1300]'; Zbw = [50 50]';
dur = .5; % in s
pNoiseVar = 10;
snr_dB = 25;
cepOrder_vect = 50:-1:8;
fs = 16e3;
plot_flag = 0;
algFlag = [1 0]; % Select 1 to run, 0 not to; [EKF EKS]
x0 = [F; Z]+100;

rmse = zeros(1, length(cepOrder_vect));

% loop through
for ii = 1:length(cepOrder_vect)
    cepOrder = cepOrder_vect(ii);
    tic
    rmse(ii) = runSynth_ARMApq(F, Fbw, Z, Zbw, dur, pNoiseVar, snr_dB, ...
        cepOrder, fs, plot_flag, algFlag, x0);
    toc
end

% plot RMSE vs frequency spacing
figure
plot(cepOrder_vect, rmse, 'b')
xlabel('# cepstral coefficents')
ylabel('Average RMSE')