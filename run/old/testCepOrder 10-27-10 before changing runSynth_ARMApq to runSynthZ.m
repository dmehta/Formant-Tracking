%% RMSE vs number of cepstral coefficients
clear

%% parameters
F = [500 1500 2500 3500]';
Fbw = [100 100 100 100]';
Z = [700 1300 2700]'; Zbw = [50 50 50]';
dur = .5; % in s
pNoiseVar = 10;
snr_dB = 25;
cepOrder_vect = 25:-1:23;
fs = 16e3;
trackBW = 0;
plot_flag = 0;
algFlag = [0 1]; % Select 1 to run, 0 not to; [EKF EKS]
x0 = [F; Z]+0;

% analysis parameters
aParams.wType = 'hamming';      % window type
aParams.wLengthMS = 20;         % Length of window (in milliseconds)
aParams.wOverlap = 0.5;         % Factor of overlap of window
aParams.lpcOrder = length(F)*2; % Number of AR coefficients
aParams.zOrder = length(Z)*2;   % Number of MA coefficients
aParams.peCoeff = 0;            % Pre-emphasis factor (0.9, 0.7, etc.)
aParams.fs = fs;                % sampling rate (in Hz)

numTrials = 10;
rmse = zeros(numTrials, length(cepOrder_vect));

% loop through
for ii = 1:length(cepOrder_vect)
    for jj = 1:numTrials
        fprintf('.')
        cepOrder = cepOrder_vect(ii);
        rmse(jj, ii) = runSynth_ARMApq(F, Fbw, Z, Zbw, dur, pNoiseVar, snr_dB, cepOrder, fs, trackBW, plot_flag, algFlag, x0, aParams);
    end
end

%% plot RMSE vs cepstral order
xdata = cepOrder_vect;
figure, box off, hold on
[ydata_lower ydata_upper ydata] = findCI(rmse, 68);
fill([xdata xdata(end:-1:1)], [ydata_lower ydata_upper(end:-1:1)], [0.9 0.9 0.9], 'EdgeColor', 'none')
plot(xdata, ydata, 'b-', 'MarkerFace', 'b', 'MarkerSize', 1, 'LineWidth', 1)
xlabel('# cepstral coefficents')
ylabel('Average RMSE (Hz)')

% save('../results/testCepOrder_1-40.mat')