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
dur = .5; % in s
snr_dB = 25;
cepOrder = 15;
fs = 10e3;
trackBW = 0;
plot_flag = 1;
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
% Fbeg = [500 1500 2500]';
% Fend = [1000 1100 3000]';
% Fcontour = interp1( [1 numFrames]', [Fbeg, Fend]', 1:numFrames );
% Fbw = [100 100 100]';
% 
% Zbeg = [4000]';
% Zend = [4200]';
% Zcontour = interp1( [1 numFrames]', [Zbeg, Zend]', 1:numFrames )'; % if one anti-formant, need to take transpose here
% Zbw = [100]';

% Zcontour = []'; Zbw = []';

%% sinewave modulated trajectory
Fbeg = [500 1500 2500]';
Fbw = [100 100 100]';
Fpert = 75*repmat(cos(2*pi*4*[1:numFrames]/numFrames*dur)', 1, size(Fbeg, 1));
Fcontour = interp1( [1 numFrames]', [Fbeg Fbeg]', 1:numFrames );
Fcontour = Fcontour + Fpert;

Zbeg = [800]';
Zbw = [100]';
Zpert = 100*repmat(cos(2*pi*4*[1:numFrames]/numFrames*dur)', 1, size(Zbeg, 1));
Zcontour = interp1( [1 numFrames]', [Zbeg Zbeg]', 1:numFrames );
Zcontour = Zcontour' + Zpert;

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
[rmse_EKS, x_estEKS] = runSynth_OLA(Fcontour, Fbw, Zcontour, Zbw, N, snr_dB, ...
        cepOrder, fs, trackBW, plot_flag, algFlag, x0);

%%
figure, hold on
plot(x_estEKS(1:size(Fcontour,2), :)', 'b')
plot(x_estEKS(size(Fcontour,2)+1:end, :)', 'r')
xlabel('Frame')
ylabel('Frequency (Hz)')
title('Estimated EKS trajectories')
rmse_EKS