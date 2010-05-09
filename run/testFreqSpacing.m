%% RMSE vs spacing between center frequencies
clear
close all

%% parameters for pole/pole
clear
close all

F1 = 1000; F2vect = 2000:-1:1000;
Fbw = [100 100]';
Z = []; Zbw = [];
dur = .5; % in s
pNoiseVar = 10;
snr_dB = 25;
cepOrder = 15;
fs = 16e3;
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
            cepOrder, fs, plot_flag, algFlag, x0);
    end
end

% plot RMSE vs frequency spacing
xdata = F2vect-F1;
ydata = mean(rmse, 1);
yerror = std(rmse, 0, 1);

ydata_upper = ydata + yerror;
ydata_lower = ydata - yerror;

figure, fill([xdata xdata(end:-1:1)], [ydata_lower ydata_upper(end:-1:1)], [0.9 0.9 0.9], 'EdgeColor', 'none')
hold on, plot(xdata, ydata, 'b-', 'MarkerFace', 'b', 'MarkerSize', 1, 'LineWidth', 1)
box off

xlabel('Spacing (Hz)')
ylabel('Average RMSE')
title('Two resonances')

save('../testFreqSpacing_results/pole-pole.mat')

%% parameters for zero/zero
clear
close all

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

numTrials = 2;
rmse = zeros(numTrials, length(Z2vect));

% loop through Z2
for ii = 1:length(Z2vect)
    for jj = 1:numTrials    
        Z = [Z1 Z2vect(ii)]';
        x0 = [F; Z]+0;
        rmse(jj, ii) = runSynth_ARMApq(F, Fbw, Z, Zbw, dur, pNoiseVar, snr_dB, ...
            cepOrder, fs, plot_flag, algFlag, x0);
    end
end

% plot RMSE vs frequency spacing
xdata = Z2vect-Z1;
ydata = mean(rmse, 1);
yerror = std(rmse, 0, 1);

ydata_upper = ydata + yerror;
ydata_lower = ydata - yerror;

figure, fill([xdata xdata(end:-1:1)], [ydata_lower ydata_upper(end:-1:1)], [0.9 0.9 0.9], 'EdgeColor', 'none')
hold on, plot(xdata, ydata, 'b-', 'MarkerFace', 'b', 'MarkerSize', 1, 'LineWidth', 1)
box off

xlabel('Spacing (Hz)')
ylabel('Average RMSE')
title('Two anti-resonances')

save('../testFreqSpacing_results/zero-zero.mat')

%% parameters for pole/zero
clear
close all

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

numTrials = 10;
rmse = zeros(1, length(Zvect));

% loop through
for ii = 1:length(Zvect)
    for jj = 1:numTrials
        Z = Zvect(ii);
        x0 = [F; Z]+0;
        rmse(jj, ii) = runSynth_ARMApq(F, Fbw, Z, Zbw, dur, pNoiseVar, snr_dB, ...
            cepOrder, fs, plot_flag, algFlag, x0);
    end
end

% plot RMSE vs frequency spacing
figure
xdata = Zvect-F;
ydata = mean(rmse, 1);
yerror = std(rmse, 0, 1);

ydata_upper = ydata + yerror;
ydata_lower = ydata - yerror;

figure, fill([xdata xdata(end:-1:1)], [ydata_lower ydata_upper(end:-1:1)], [0.9 0.9 0.9], 'EdgeColor', 'none')
hold on, plot(xdata, ydata, 'b-', 'MarkerFace', 'b', 'MarkerSize', 1, 'LineWidth', 1)
box off

xlabel('Spacing (Hz)')
ylabel('Average RMSE')
title('One resonance, one anti-resonance')

save('../testFreqSpacing_results/pole-zero.mat')