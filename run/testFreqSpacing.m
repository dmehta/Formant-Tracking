%% RMSE vs spacing between center frequencies
clear
close all

%% parameters for pole/pole
F1 = 1000; F2vect = 2000:-100:1000;
Fbw = [100 100]';
Z = []; Zbw = [];
dur = .5; % in s
pNoiseVar = 10;
snr_dB = 25;
cepOrder = 15;
fs = 16e3;
plot_flag = 0;
algFlag = [0 1]; % Select 1 to run, 0 not to; [EKF EKS]

rmse = zeros(1, length(F2vect));

% loop through F2
for ii = 1:length(F2vect)
    F = [F1 F2vect(ii)]';
    x0 = [F; Z]+0;
    rmse(ii) = runSynth_ARMApq(F, Fbw, Z, Zbw, dur, pNoiseVar, snr_dB, ...
        cepOrder, fs, plot_flag, algFlag, x0);
end

% plot RMSE vs frequency spacing
figure
plot(F2vect-F1, rmse)
xlabel('Spacing (Hz)')
ylabel('Average RMSE')
title('Two resonances')

%% parameters for zero/zero
F = 500; Fbw = 100;
Z1 = 1000; Z2vect = 2000:-100:1000;
Zbw = [100 100]';
dur = .5; % in s
pNoiseVar = 10;
snr_dB = 25;
cepOrder = 15;
fs = 16e3;
plot_flag = 0;
algFlag = [0 1]; % Select 1 to run, 0 not to; [EKF EKS]

rmse = zeros(1, length(Z2vect));

% loop through Z2
for ii = 1:length(Z2vect)
    Z = [Z1 Z2vect(ii)]';
    x0 = [F; Z]+0;
    rmse(ii) = runSynth_ARMApq(F, Fbw, Z, Zbw, dur, pNoiseVar, snr_dB, ...
        cepOrder, fs, plot_flag, algFlag, x0);
end

% plot RMSE vs frequency spacing
figure
plot(Z2vect-Z1, rmse)
xlabel('Spacing (Hz)')
ylabel('Average RMSE')
title('Two anti-resonances')

%% parameters for pole/zero
F = 1000;
Fbw = [100]';
Zvect = 2000:-100:1000; Zbw = [50];
dur = .5; % in s
pNoiseVar = 10;
snr_dB = 25;
cepOrder = 15;
fs = 16e3;
plot_flag = 0;
algFlag = [0 1]; % Select 1 to run, 0 not to; [EKF EKS]

rmse = zeros(1, length(Zvect));

% loop through
for ii = 1:length(Zvect)
    Z = Zvect(ii);
    x0 = [F; Z]+0;
    rmse(ii) = runSynth_ARMApq(F, Fbw, Z, Zbw, dur, pNoiseVar, snr_dB, ...
        cepOrder, fs, plot_flag, algFlag, x0);
end

% plot RMSE vs frequency spacing
figure
plot(Zvect-F, rmse)
xlabel('Spacing (Hz)')
ylabel('Average RMSE')
title('One resonance, one anti-resonance')
