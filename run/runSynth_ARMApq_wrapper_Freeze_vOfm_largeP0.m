%% runSynth_ARMApq_wrapper.m
% set parameters for runSynth_ARMApq.m
% 
% INPUT:
%    F:         center frequencies of the resonances (col vector), in Hz
%    Fbw:       corresponding bandwidths of the resonances (col vector), in Hz
%    Z:         center frequencies of the anti-resonances (col vector), in Hz
%    Zbw:       corresponding bandwidths of the anti-resonances (col vector), in Hz
%    dur:       duration of signal, in s
%    pNoiseVar: process noise variance
%    snr_dB:    observation noise, in dB
%    cepOrder:  Number of cepstal coefficients to compute
%    fs:        sampling rate of waveform, in Hz
%    trackBW:   track bandwidths if 1
%    plot_flag: plot figures if 1
%    algFlag:   select 1 to run, 0 not to for [EKF EKS]
%    x0:        initial state of formant trackers [F;Z], in Hz
% 
% OUTPUT:
%    Depends on algFlag. For each algorithm, three outputs generated--
%       rmse_mean:  average RMSE across all tracks
%       x_est:      estimated tracks
%       x_errVar:   covariance matrix of estimated tracks
%       So that if two algorithms run, the following are output:
%       [rmse_meanEKF, x_estEKF, x_errVarEKF, rmse_meanEKS, x_estEKS, x_errVarEKS]

clear 

%% parameters
F = [500 1000]';
Fbw = [50 50]';
Z = []'; Zbw = []';

dur = .25; % in s
pNoiseVar = 4;
snr_dB = 25;
cepOrder = 20;
fs = 10e3;
trackBW = 1;
plot_flag = 0;
algFlag = [0 1]; % Select 1 to run, 0 not to; [EKF EKS]

if trackBW
    x0 = [F; Fbw; Z; Zbw];
else
    x0 = [F; Z];
end

%%    
numTrials = 100;
rmse = zeros(numTrials, 1);
x_est = cell(numTrials, 1);
x_errVar = cell(numTrials, 1);

for jj = 1:numTrials
    [rmse(jj), x_est{jj}, x_errVar{jj}] = runSynth_ARMApq(F, Fbw, Z, Zbw, dur, pNoiseVar, snr_dB, ...
            cepOrder, fs, trackBW, plot_flag, algFlag, x0);
end

numObs = size(x_est{1}, 2);
numStates = size(x_est{1}, 1);

if trackBW
    trueState = [repmat(F, 1, numObs); repmat(Fbw, 1, numObs); ...
        repmat(Z, 1, numObs); repmat(Zbw, 1, numObs)];
else
    trueState = repmat([F; Z], 1, numObs);
end

%% plot all simulations against truth
figure, hold on
for jj = 1:numTrials
    plot(x_est{jj}(1:size(F,1), :)', 'b')
    plot(x_est{jj}(size(F,1)+1:end, :)', 'r')
end
plot(trueState', 'LineWidth', 3, 'Color', 'k')
xlabel('Frame')
ylabel('Frequency (Hz)')
title('Estimated EKS trajectories (Formant-blue, Antiformant-red, True-black)')

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
xlabel('Frame')
ylabel('Frequency (Hz)')

figure(43)
xlabel('Frame')
ylabel('Variance (Hz^2)')
legend('Variance of means', 'Mean of variances')

%% plot tracks individually in grid style with confidence intervals
titleCell(1,2) = {'EKS'}; % hard coded for now
titleCell(2,2) = {'b:'};
titleCell(1,1)  = {'True'};   % Keeps track of trackers used for plotter
titleCell(2,1)  = {'r'};      % Color for true state plot
nP = size(F,1);

plotStateTracksFZ_CI(trueState,x_est,titleCell(:,[1 2]), nP, trackBW);

disp(['Mean RMSE: ', num2str(mean(rmse))])