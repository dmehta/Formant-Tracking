%% RMSE vs spacing between center frequencies
clear

%% parameters for pole/pole
clear

F1 = 1000; F2vect = 2000:-10:1000;
Fbw = [100 100]';
Z = []; Zbw = [];
dur = .5; % in s
pNoiseVar = 10;
snr_dB = 25;
cepOrder = 15;
fs = 16e3;
trackBW = 0;
plot_flag = 0;
algFlag = [0 1]; % Select 1 to run, 0 not to; [EKF EKS]

numTrials = 10;
rmse = zeros(numTrials, length(F2vect));

% loop through F2
for ii = 1:length(F2vect)
    for jj = 1:numTrials
        F = [F1 F2vect(ii)]';
        x0 = [F; Z]+0;
        rmse(jj, ii) = runSynth_ARMApq(F, Fbw, Z, Zbw, dur, pNoiseVar, snr_dB, ...
            cepOrder, fs, trackBW, plot_flag, algFlag, x0);
    end
end

%% plot RMSE vs frequency spacing
xdata = F2vect-F1;
figure, box off, hold on
[ydata_lower ydata_upper ydata] = findCI(rmse, 95);
fill([xdata xdata(end:-1:1)], [ydata_lower ydata_upper(end:-1:1)], [0.9 0.9 0.9], 'EdgeColor', 'none')
plot(xdata, ydata, 'b-', 'MarkerFace', 'b', 'MarkerSize', 1, 'LineWidth', 1)
xlabel('Spacing (Hz)')
ylabel('Average RMSE')
title('Two resonances')

% save('../results/testFreqSpacing_pole-pole.mat')

%% parameters for zero/zero
clear

F = 500; Fbw = 100;
Z1 = 1000; Z2vect = 2000:-10:1000;
Zbw = [100 100]';
dur = .5; % in s
pNoiseVar = 10;
snr_dB = 25;
cepOrder = 15;
fs = 16e3;
trackBW = 0;
plot_flag = 0;
algFlag = [0 1]; % Select 1 to run, 0 not to; [EKF EKS]

numTrials = 10;
rmse = zeros(numTrials, length(Z2vect));

% loop through Z2
for ii = 1:length(Z2vect)
    for jj = 1:numTrials    
        Z = [Z1 Z2vect(ii)]';
        x0 = [F; Z]+0;
        rmse(jj, ii) = runSynth_ARMApq(F, Fbw, Z, Zbw, dur, pNoiseVar, snr_dB, ...
            cepOrder, fs, trackBW, plot_flag, algFlag, x0);
    end
end

%% plot RMSE vs frequency spacing
xdata = Z2vect-Z1;
figure, box off, hold on
[ydata_lower ydata_upper ydata] = findCI(rmse, 95);
fill([xdata xdata(end:-1:1)], [ydata_lower ydata_upper(end:-1:1)], [0.9 0.9 0.9], 'EdgeColor', 'none')
plot(xdata, ydata, 'b-', 'MarkerFace', 'b', 'MarkerSize', 1, 'LineWidth', 1)
xlabel('Spacing (Hz)')
ylabel('Average RMSE')
title('Two anti-resonances')

% save('../results/testFreqSpacing_zero-zero.mat')

%% parameters for pole/zero
clear

F = 1000;
Fbw = [100]';
Zvect = 2000:-10:1000; Zbw = [50];
dur = .5; % in s
pNoiseVar = 10;
snr_dB = 25;
cepOrder = 15;
fs = 16e3;
trackBW = 0;
plot_flag = 0;
algFlag = [0 1]; % Select 1 to run, 0 not to; [EKF EKS]

numTrials = 10;
rmse = zeros(1, length(Zvect));

% loop through
for ii = 1:length(Zvect)
    for jj = 1:numTrials
        Z = Zvect(ii);
        x0 = [F; Z]+0;
        rmse(jj, ii) = runSynth_ARMApq(F, Fbw, Z, Zbw, dur, pNoiseVar, snr_dB, ...
            cepOrder, fs, trackBW, plot_flag, algFlag, x0);
    end
end

%% plot RMSE vs frequency spacing
xdata = Zvect-F;

figure, box off, hold on
[ydata_lower ydata_upper ydata] = findCI(rmse, 95);
fill([xdata xdata(end:-1:1)], [ydata_lower ydata_upper(end:-1:1)], [0.9 0.9 0.9], 'EdgeColor', 'none')
plot(xdata, ydata, 'b-', 'MarkerFace', 'b', 'MarkerSize', 1, 'LineWidth', 1)
xlabel('Spacing (Hz)')
ylabel('Average RMSE')
title('One resonance, one anti-resonance')

% save('../results/testFreqSpacing_pole-zero.mat')