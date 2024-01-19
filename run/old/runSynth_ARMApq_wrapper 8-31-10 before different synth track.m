%% parameters

% synthesis parameters
dur = .05; % in s
fs = 8e3; % in Hz

F = [500 1000]'; % in Hz, time-invariant formant center frequency (F)
Fbw = [50 50]'; % in Hz, time-invariant formant bandwidth (Fbw)
Z = []'; % in Hz, time-invariant formant center frequency (Z)
Zbw = []'; % in Hz, time-invariant formant bandwidth (Zbw)

% analysis parameters
aParams.wType = 'hamming';      % window type
aParams.wLengthMS = 20;         % Length of window (in milliseconds)
aParams.wOverlap = 0.5;         % Factor of overlap of window
aParams.lpcOrder = length(F)*2; % Number of AR coefficients
aParams.zOrder = length(Z)*2;   % Number of MA coefficients
aParams.peCoeff = 0;            % Pre-emphasis factor (0.9, 0.7, etc.)
aParams.fs = fs;                % sampling rate (in Hz)

% tracker parameters
pNoiseVar = 4;
snr_dB = 30;
cepOrder = 2; %max(aParams.lpcOrder, aParams.zOrder);
trackBW = 1;
plot_flag = 0; % plot transfer functions
algFlag = [0 1]; % Select 1 to run, 0 not to; [EKF EKS]
offset = 0; % set initial state offset, in Hz

% Monte Carlo analysis parameters
numTrials = 1;
trial = 1; % pick a trial for which plot spectrogram

%% misc calculations
if trackBW
    x0 = [F; Fbw; Z; Zbw]+offset;
else
    x0 = [F; Z]+offset;
end

%% synthesize and track
rmse = zeros(numTrials, 1);
x_est = cell(numTrials, 1);
x_errVar = cell(numTrials, 1);
x = cell(numTrials, 1);

for jj = 1:numTrials
    disp(['Processing Trial #', num2str(jj), '...'])
    [rmse{jj}, x_est{jj}, x_errVar{jj}, x{jj}] = runSynth_ARMApq(F, Fbw, Z, Zbw, dur, pNoiseVar, snr_dB, ...
            cepOrder, fs, trackBW, plot_flag, algFlag, x0, aParams);
end

numObs = size(x_est{1}, 2);
numStates = size(x_est{1}, 1);

if trackBW
    trueState = [repmat(F, 1, numObs); repmat(Fbw, 1, numObs); ...
        repmat(Z, 1, numObs); repmat(Zbw, 1, numObs)];
else
    trueState = repmat([F; Z], 1, numObs);
end

%% plot tracks individually in grid style with confidence intervals
titleCell(1,2) = {'EKS'}; % hard coded for now
titleCell(2,2) = {'b:'};
titleCell(1,1)  = {'True'};   % Keeps track of trackers used for plotter
titleCell(2,1)  = {'r'};      % Color for true state plot
nP = size(F,1);

plotStateTracksFZ_CI(trueState,x_est,titleCell(:,[1 2]), nP, trackBW);

disp(['Mean RMSE: ', num2str(mean(rmse))])

%% Super-impose over a spectrogram
nZ = size(Z,1);
aParams.wLength = floor(aParams.wLengthMS/1000*aParams.fs);
aParams.wLength = aParams.wLength + (mod(aParams.wLength,2)~=0); % Force even

plotSpecTracks2BW(x{trial}, x_est{trial}, aParams, nZ, trackBW);
firstline = get(get(gca,'Title'),'String');
title({['Trial #', num2str(trial)]; firstline})
axis tight

%% plot frequency vs frame with confidence intervals
% repackage
x_estPerFreq = zeros(numTrials, numObs, numStates);
x_errVarPerFreq = zeros(numTrials, numObs, numStates);

for kk = 1:numStates
    for jj = 1:numTrials
        x_estPerFreq(jj, :, kk) = x_est{jj}(kk, :);
        x_errVarPerFreq(jj, :, kk) = x_errVar{jj}(kk, kk, :);
    end
end

obs = 1:numObs;
means = zeros(numStates, numObs);
vOfm = zeros(numStates, numObs);
mOfv = zeros(numStates, numObs);

figure(42), box off, hold on
figure(43), box off, hold on

for kk = 1:numStates
    figure(42)
    [L U ave v] = findCI(x_estPerFreq(:,:,kk), 95);
    fill([obs obs(end:-1:1)], [L U(end:-1:1)], [0.9 0.9 0.9], 'EdgeColor', 'none')
    plot(obs, ave, 'b-', 'MarkerFace', 'b', 'MarkerSize', 1, 'LineWidth', 1)

    means(kk,:) = ave;
    vOfm(kk,:) = v;
    mOfv(kk,:) = mean(x_errVarPerFreq(:,:,kk), 1);
    
    figure(43)
    plot(vOfm(kk,:))
    plot(mOfv(kk,:), 'r')
end

figure(42)
title('Tracks with 95 % CI')
xlabel('Frame')
ylabel('Frequency (Hz)')

figure(43)
title('Comparison of confidence interval and KF covariances')
xlabel('Frame')
ylabel('Variance (Hz^2)')
legend('Variance of means', 'Mean of variances')