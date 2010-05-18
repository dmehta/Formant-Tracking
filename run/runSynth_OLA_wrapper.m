% set parameters for runSynth_OLA.m
% 
% INPUT:
%    Fcontour:  center frequencies of the resonances (numFrames x numFormants), in Hz
%    Fbw:       corresponding bandwidths of the resonances (numFormants x 1), in Hz
%    Z:         center frequencies of the anti-resonances (numFrames x numAntiformants), in Hz
%    Zbw:       corresponding bandwidths of the anti-resonances (numAntiformants x 1), in Hz
%    N:         length of signal, in samples
%    snr_dB:    observation noise, in dB
%    cepOrder:  Number of cepstal coefficients to compute
%    fs:        sampling rate of waveform, in Hz
%    plot_flag: plot figures if 1
%    algFlag:   select 1 to run, 0 not to for [EKF EKS]
%    x0:        initial states (center frequencies) of formant trackers [(numFormants + numAntiformants) x 1], in Hz
% 
% OUTPUT:
%    Depends on algFlag. For each algorithm, two outputs generated--
%       rmse_mean1: average RMSE across all tracks
%       x_est1:  estimated tracks
%       So that if two algorithms run, the following are output:
%       [rmse_mean1, x_est1, rmse_mean2, x_est2]

%% parameters
dur = .25; % in s
snr_dB = 25;
cepOrder = 15;
fs = 10e3;
trackBW = 0;
plot_flag = 0;
algFlag = [0 1]; % Select 1 to run, 0 not to; [EKF EKS]
wType = 'hanning';  % window type
wLengthMS  = 20;    % Length of window (in milliseconds)
wOverlap = 0;     % Factor of overlap of window

wLength = floor(wLengthMS/1000*fs);
wLength = wLength + (mod(wLength,2)~=0); % Force even
win = feval(wType,wLength);
N = floor(dur*fs);
N = N - mod(N,wLength);
numFrames = floor(N/((1-wOverlap)*wLength))-1;

%% Piecewise linear trajectory
Fbeg = [500 1500 2500]';
Fend = [1000 1100 3000]';
Fcontour = interp1( [1 numFrames]', [Fbeg, Fend]', 1:numFrames );
Fcontour = Fcontour + repmat(10*randn(size(Fcontour,1),1), 1, size(Fcontour,2));
Fbw = [100 100 100]';

Zbeg = [4000]';
Zend = [4000]';
Zcontour = interp1( [1 numFrames]', [Zbeg, Zend]', 1:numFrames )'; % if one anti-formant, need to take transpose here
Zcontour = Zcontour + repmat(10*randn(size(Zcontour,1),1), 1, size(Zcontour,2));
Zbw = [100]';

% Zcontour = []'; Zbw = []';

%% sinewave modulated trajectory
Fbeg = [500 1500 2500]';
Fbw = [100 100 100]';
Fpert = 75*repmat(cos(2*pi*4*[0:numFrames-1]/numFrames*dur)', 1, size(Fbeg, 1));
Fcontour = interp1( [1 numFrames]', [Fbeg Fbeg]', 1:numFrames );
Fcontour = Fcontour + Fpert + repmat(10*randn(size(Fcontour,1),1), 1, size(Fcontour,2));

Zbeg = [800]';
Zbw = [100]';
Zpert = 100*repmat(cos(2*pi*4*[0:numFrames-1]/numFrames*dur)', 1, size(Zbeg, 1));
Zcontour = interp1( [1 numFrames]', [Zbeg Zbeg]', 1:numFrames )';
Zcontour = Zcontour + Zpert + repmat(10*randn(size(Zcontour,1),1), 1, size(Zcontour,2));

% Fcontour = []'; Fbw = []';
% Zcontour = []'; Zbw = []';

%% initial state
if ~isempty(Fcontour) && ~isempty(Zcontour)
    x0 = [Fcontour(1, :)'; Zcontour(1, :)']+100;
end

if ~isempty(Fcontour) && isempty(Zcontour)
    x0 = Fcontour(1, :)'+100;
end

if isempty(Fcontour) && ~isempty(Zcontour)
    x0 = Zcontour(1, :)'+100;
end

%%
numTrials = 10;
rmse = zeros(numTrials, 1);
x_est = cell(numTrials,1);

for jj = 1:numTrials
    [rmse(jj), x_est{jj}] = runSynth_OLA(Fcontour, Fbw, Zcontour, Zbw, N, snr_dB, ...
            cepOrder, fs, trackBW, plot_flag, algFlag, x0);
end

% save('../results/lines_20sims.mat')
% save('../results/sines_20sims.mat')

%% plot all simulations against truth
trueState = [Fcontour'; Zcontour'];

figure, hold on
for jj = 1:numTrials
    plot(x_est{jj}(1:size(Fcontour,2), :)', 'b')
    plot(x_est{jj}(size(Fcontour,2)+1:end, :)', 'r')
end
plot(trueState', 'LineWidth', 3, 'Color', 'k')
xlabel('Frame')
ylabel('Frequency (Hz)')
title('Estimated EKS trajectories (Formant-blue, Antiformant-red, True-black)')
mean(rmse)

%% plot frequency vs frame
% repackage
x_estPerFreq = zeros(numTrials, numFrames, size(trueState,1));
for kk = 1:size(trueState,1)
    for jj = 1:numTrials
        x_estPerFreq(jj, :, kk) = x_est{jj}(kk, :);
    end
end

xdata = 1:numFrames;
figure, box off, hold on
for kk = 1:size(trueState,1)
    [ydata_lower ydata_upper ydata] = findCI(x_estPerFreq(:,:,kk), 95);
    fill([xdata xdata(end:-1:1)], [ydata_lower ydata_upper(end:-1:1)], [0.9 0.9 0.9], 'EdgeColor', 'none')
    plot(xdata, ydata, 'b-', 'MarkerFace', 'b', 'MarkerSize', 1, 'LineWidth', 1)
    xlabel('Frame')
    ylabel('Frequency (Hz)')
end

%% plot tracks individually in grid style
titleCell(1,2) = {'EKS'}; % hard coded for now
titleCell(2,2) = {'b:'};
titleCell(1,1)  = {'True'};   % Keeps track of trackers used for plotter
titleCell(2,1)  = {'r'};      % Color for true state plot
nP = size(Fcontour,2);

plotStateTracksFZ_CI(trueState,x_est,titleCell(:,[1 2]), nP);