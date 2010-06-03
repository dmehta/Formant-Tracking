% set parameters for runSynth_OLA.m
% 
% INPUT:
%    F:         center frequencies of the resonances (numFormants x numFrames), in Hz
%    Fbw:       corresponding bandwidths of the resonances (numFormants x numFrames), in Hz
%    Z:         center frequencies of the anti-resonances (numAntiformants x numFrames), in Hz
%    Zbw:       corresponding bandwidths of the anti-resonances (numAntiformants x numFrames), in Hz
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
if trackBW
    Fbeg = [500 1500 2500]';
    Fend = [1000 1100 3000]';
    F = interp1([1 numFrames]', [Fbeg, Fend]', 1:numFrames)'; if length(Fbeg) == 1, F=F'; end;
    F = F + 10*randn(size(F));

    Fbwbeg = [100 100 100]';
    Fbwend = [200 200 200]';
    Fbw = interp1([1 numFrames]', [Fbwbeg, Fbwend]', 1:numFrames)'; if length(Fbwbeg) == 1, Fbw=Fbw'; end;
    Fbw = Fbw + 10*randn(size(Fbw));

    Zbeg = [4000]';
    Zend = [4500]';
    Z = interp1([1 numFrames]', [Zbeg, Zend]', 1:numFrames)'; if length(Zbeg) == 1, Z=Z'; end;
    Z = Z + 10*randn(size(Z));

    Zbwbeg = [100]';
    Zbwend = [200]';
    Zbw = interp1([1 numFrames]', [Zbwbeg, Zbwend]', 1:numFrames)'; if length(Zbwbeg) == 1, Zbw=Zbw'; end;
    Zbw = Zbw + 10*randn(size(Zbw));
    
else
    Fbeg = [500 1500 2500]';
    Fend = [1000 1100 3000]';
    F = interp1([1 numFrames]', [Fbeg, Fend]', 1:numFrames)'; if length(Fbeg) == 1, F=F'; end;
    F = F + 10*randn(size(F));
    Fbw = [100 100 100]';

    Zbeg = [4000]';
    Zend = [4500]';
    Z = interp1([1 numFrames]', [Zbeg, Zend]', 1:numFrames)'; if length(Zbeg) == 1, Z=Z'; end;
    Z = Z + 10*randn(size(Z));
    Zbw = [100]';
end

%% sinewave modulated trajectory
if trackBW
    Fbeg = [500 1000]';
    Fpert = 75*repmat(cos(2*pi*8*[0:numFrames-1]/numFrames*dur), size(Fbeg, 1), 1);
    F = interp1([1 numFrames]', [Fbeg Fbeg]', 1:numFrames)'; if length(Fbeg) == 1, F=F'; end;
    F = F + Fpert + 10*randn(size(F));

    Fbwbeg = [80 80]';
    Fbwpert = 10*repmat(cos(2*pi*10*[0:numFrames-1]/numFrames*dur), size(Fbwbeg, 1), 1);
    Fbw = interp1([1 numFrames]', [Fbwbeg Fbwbeg]', 1:numFrames)'; if length(Fbwbeg) == 1, Fbw=Fbw'; end;
    Fbw = Fbw + Fbwpert + 10*randn(size(Fbw));

    Zbeg = [2000 4000]';
    Zpert = 100*repmat(cos(2*pi*7*[0:numFrames-1]/numFrames*dur), size(Zbeg, 1), 1);
    Z = interp1([1 numFrames]', [Zbeg Zbeg]', 1:numFrames)'; if length(Zbeg) == 1, Z=Z'; end;
    Z = Z + Zpert + 10*randn(size(Z));

    Zbwbeg = [100 100]';
    Zbwpert = 10*repmat(cos(2*pi*6*[0:numFrames-1]/numFrames*dur), size(Zbwbeg, 1), 1);
    Zbw = interp1([1 numFrames]', [Zbwbeg Zbwbeg]', 1:numFrames)'; if length(Zbwbeg) == 1, Zbw=Zbw'; end;
    Zbw = Zbw + Zbwpert + 10*randn(size(Zbw));
else
    Fbeg = [500 1500 2500]';
    Fbw = [100 100 100]';
    Fpert = 75*repmat(cos(2*pi*4*[0:numFrames-1]/numFrames*dur), size(Fbeg, 1), 1);
    F = interp1([1 numFrames]', [Fbeg Fbeg]', 1:numFrames)'; if length(Fbeg) == 1, F=F'; end;
    F = F + Fpert + 10*randn(size(F));

    Zbeg = [800]';
    Zbw = [100]';
    Zpert = 10*repmat(cos(2*pi*4*[0:numFrames-1]/numFrames*dur), size(Zbeg, 1), 1);
    Z = interp1([1 numFrames]', [Zbeg Zbeg]', 1:numFrames)'; if length(Zbeg) == 1, Z=Z'; end;
    Z = Z + Zpert + 10*randn(size(Z));
end

% F = []'; Fbw = []';
% Z = []'; Zbw = []';

%% initial state
if ~isempty(F) && ~isempty(Z)
    if trackBW
        x0 = [F(:, 1); Fbw(:, 1); Z(:, 1); Zbw(:, 1)]+0;
    else
        x0 = [F(:, 1); Z(:, 1)]+0;
    end
end

if ~isempty(F) && isempty(Z)
    if trackBW
        x0 = [F(:, 1); Fbw(:, 1)]+0;
    else
        x0 = F(:, 1)+0;
    end
end

if isempty(F) && ~isempty(Z)
    if trackBW
        x0 = [Z(:, 1); Zbw(:, 1)]+0;
    else
        x0 = Z(:, 1)+0;
    end
end

%%
numTrials = 10;
rmse = zeros(numTrials, 1);
x_est = cell(numTrials, 1);

for jj = 1:numTrials
    tic
    [rmse(jj), x_est{jj}] = runSynth_OLA(F, Fbw, Z, Zbw, N, snr_dB, ...
            cepOrder, fs, trackBW, plot_flag, algFlag, x0);        
    toc
end

if trackBW
    trueState = [F; Fbw; Z; Zbw];
else
    trueState = [F; Z]; % 
end

% save('../results/lines_20sims.mat')
% save('../results/sines_20sims.mat')

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
nP = size(F,1);

plotStateTracksFZ_CI(trueState,x_est,titleCell(:,[1 2]), nP, trackBW);

disp(['Mean RMSE: ', num2str(mean(rmse))])