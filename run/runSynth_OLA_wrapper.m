% set parameters for runSynth_OLA.m
% 
% INPUT:
%    F:         center frequencies of the resonances (numFrames x numFormants), in Hz
%    Fbw:       corresponding bandwidths of the resonances (numFrames x numFormants), in Hz
%    Z:         center frequencies of the anti-resonances (numFrames x numAntiformants), in Hz
%    Zbw:       corresponding bandwidths of the anti-resonances (numFrames x numAntiformants), in Hz
%    N:         length of signal, in samples
%    snr_dB:    observation noise, in dB
%    cepOrder:  Number of cepstal coefficients to compute
%    fs:        sampling rate of waveform, in Hz
%    trackBW:   track bandwidths if 1
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

clear

%% parameters
dur = .25; % in s
snr_dB = 25;
cepOrder = 15;
fs = 10e3;
trackBW = 1;
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
F = interp1( [1 numFrames]', [Fbeg, Fend]', 1:numFrames );
F = F + repmat(10*randn(size(F,1),1), 1, size(F,2));
Fbw = [100 100 100]';

Zbeg = [4000]';
Zend = [4000]';
Z = interp1( [1 numFrames]', [Zbeg, Zend]', 1:numFrames )'; % if one anti-formant, need to take transpose here
Z = Z + repmat(10*randn(size(Z,1),1), 1, size(Z,2));
Zbw = [100]';

% Z = []'; Zbw = []';

%% sinewave modulated trajectory
Fbeg = [500 1500 2500]';
Fbw = [100 100 100]';
Fpert = 75*repmat(cos(2*pi*4*[0:numFrames-1]/numFrames*dur)', 1, size(Fbeg, 1));
F = interp1( [1 numFrames]', [Fbeg Fbeg]', 1:numFrames );
F = F + Fpert + repmat(10*randn(size(F,1),1), 1, size(F,2));

Zbeg = [800]';
Zbw = [100]';
Zpert = 100*repmat(cos(2*pi*4*[0:numFrames-1]/numFrames*dur)', 1, size(Zbeg, 1));
Z = interp1( [1 numFrames]', [Zbeg Zbeg]', 1:numFrames )';
Z = Z + Zpert + repmat(10*randn(size(Z,1),1), 1, size(Z,2));

% F = []'; Fbw = []';
% Z = []'; Zbw = []';

%% sinewave modulated trajectory with BANDWIDTHS
Fbeg = [500 1000]';
Fpert = 75*repmat(cos(2*pi*8*[0:numFrames-1]/numFrames*dur)', 1, size(Fbeg, 1));
F = interp1( [1 numFrames]', [Fbeg Fbeg]', 1:numFrames ); % if one, need to take transpose here
F = F + Fpert + repmat(10*randn(size(F,1),1), 1, size(F,2));

Fbwbeg = [80 80]';
Fbwpert = 80*repmat(cos(2*pi*10*[0:numFrames-1]/numFrames*dur)', 1, size(Fbwbeg, 1));
Fbw = interp1( [1 numFrames]', [Fbwbeg Fbwbeg]', 1:numFrames ); % if one, need to take transpose here
Fbw = Fbw + Fbwpert + repmat(10*randn(size(Fbw,1),1), 1, size(Fbw,2));

Zbeg = [800]';
Zpert = 100*repmat(cos(2*pi*7*[0:numFrames-1]/numFrames*dur)', 1, size(Zbeg, 1));
Z = interp1( [1 numFrames]', [Zbeg Zbeg]', 1:numFrames )';
Z = Z + Zpert + repmat(10*randn(size(Z,1),1), 1, size(Z,2));

Zbwbeg = [100]';
Zbwpert = 20*repmat(cos(2*pi*6*[0:numFrames-1]/numFrames*dur)', 1, size(Zbwbeg, 1));
Zbw = interp1( [1 numFrames]', [Zbwbeg Zbwbeg]', 1:numFrames )';
Zbw = Zbw + Zbwpert + repmat(10*randn(size(Zbw,1),1), 1, size(Zbw,2));

% F = []'; Fbw = []';
% Z = []'; Zbw = []';

%% initial state
if ~isempty(F) && ~isempty(Z)
    if trackBW
        x0 = [F(1, :)'; Fbw(1, :)'; Z(1, :)'; Zbw(1, :)']+100;
    else
        x0 = [F(1, :)'; Z(1, :)']+100;
    end
end

if ~isempty(F) && isempty(Z)
    if trackBW
        x0 = [F(1, :)'; Fbw(1, :)']+100;
    else
        x0 = F(1, :)'+100;
    end
end

if isempty(F) && ~isempty(Z)
    if trackBW
        x0 = [Z(1, :)'; Zbw(1, :)']+100;
    else
        x0 = Z(1, :)'+100;
    end
end

%%
numTrials = 10;
rmse = zeros(numTrials, 1);
x_est = cell(numTrials,1);

for jj = 1:numTrials
    [rmse(jj), x_est{jj}] = runSynth_OLA(F, Fbw, Z, Zbw, N, snr_dB, ...
            cepOrder, fs, trackBW, plot_flag, algFlag, x0);        
end

if trackBW
    trueState = [F'; Fbw'; Z'; Zbw'];
else
    trueState = [F'; Z']; % 
end

% save('../results/lines_20sims.mat')
% save('../results/sines_20sims.mat')

%% plot all simulations against truth
figure, hold on
for jj = 1:numTrials
    plot(x_est{jj}(1:size(F,2), :)', 'b')
    plot(x_est{jj}(size(F,2)+1:end, :)', 'r')
end
plot(trueState', 'LineWidth', 3, 'Color', 'k')
xlabel('Frame')
ylabel('Frequency (Hz)')
title('Estimated EKS trajectories (Formant-blue, Antiformant-red, True-black)')
mean(rmse)

%% plot frequency vs frame with confidence intervals
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

%% plot tracks individually in grid style with confidence intervals
titleCell(1,2) = {'EKS'}; % hard coded for now
titleCell(2,2) = {'b:'};
titleCell(1,1)  = {'True'};   % Keeps track of trackers used for plotter
titleCell(2,1)  = {'r'};      % Color for true state plot
nP = size(F,2);

plotStateTracksFZ_CI(trueState,x_est,titleCell(:,[1 2]), nP, trackBW);