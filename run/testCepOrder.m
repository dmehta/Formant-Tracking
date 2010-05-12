%% RMSE vs number of cepstral coefficients
clear

%% parameters
F = [500 1500 2500 3500]';
Fbw = [100 100 100 100]';
Z = [700 1300 2700]'; Zbw = [50 50 50]';
dur = .5; % in s
pNoiseVar = 10;
snr_dB = 25;
cepOrder_vect = 25:-1:5;
fs = 16e3;
trackBW = 0;
plot_flag = 0;
algFlag = [0 1]; % Select 1 to run, 0 not to; [EKF EKS]
x0 = [F; Z]+0;

numTrials = 10;
rmse = zeros(numTrials, length(cepOrder_vect));

% loop through
for ii = 1:length(cepOrder_vect)
    for jj = 1:numTrials
        cepOrder = cepOrder_vect(ii);
        rmse(jj, ii) = runSynth_ARMApq(F, Fbw, Z, Zbw, dur, pNoiseVar, snr_dB, ...
            cepOrder, fs, trackBW, plot_flag, algFlag, x0);
    end
end

% plot RMSE vs cepstral order
xdata = cepOrder_vect;
ydata = mean(rmse, 1);
yerror = std(rmse, 0, 1);

ydata_upper = ydata + yerror;
ydata_lower = ydata - yerror;

figure, fill([xdata xdata(end:-1:1)], [ydata_lower ydata_upper(end:-1:1)], [0.9 0.9 0.9], 'EdgeColor', 'none')
hold on, plot(xdata, ydata, 'b-', 'MarkerFace', 'b', 'MarkerSize', 1, 'LineWidth', 1)
box off

xlabel('# cepstral coefficents')
ylabel('Average RMSE (Hz)')

save('../results/testCepOrder_1-40.mat')