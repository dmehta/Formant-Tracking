%% parameters

% synthesis parameters
dur = .25; % in s
fs = 8000; % in Hz
sParams.wType = 'hanning';      % window type
sParams.wLengthMS  = 20;       % Length of window (in milliseconds)
sParams.wOverlap = 0.5;         % Factor of overlap of window

% linear trajectory
Fbeg = [500]'; Fend = [2000]';
Fbwbeg = [100]'; Fbwend = [100]';
Zbeg = [1000]'; Zend = [1000]';
Zbwbeg = [100]'; Zbwend = [100]';

% analysis parameters
aParams.wType = sParams.wType;          % window type
aParams.wLengthMS = sParams.wLengthMS;  % Length of window (in milliseconds)
aParams.wOverlap = sParams.wOverlap;    % Factor of overlap of window
aParams.lpcOrder = length(Fbeg)*2;      % Number of AR coefficients
aParams.zOrder = length(Zbeg)*2;        % Number of MA coefficients
aParams.peCoeff = 0;                    % Pre-emphasis factor (0.9, 0.7, etc.)
aParams.fs = fs;                        % sampling rate (in Hz) to resample to

% tracker parameters
snr_dB = 20;
cepOrder = 20;
trackBW = 0;
plot_flag = 0; % plot transfer functions
algFlag = [0 1]; % Select 1 to run, 0 not to; [EKF EKS]
offset = 100; % set initial state offset, in Hz)

% Monte Carlo analysis parameters
numTrials = 20;
trial = 1; % pick a trial for which plot spectrogram

%% misc calculations
wLength = floor(aParams.wLengthMS/1000*fs);
wLength = wLength + (mod(wLength,2)~=0); % Force even
win = feval(aParams.wType,wLength);
N = floor(dur*aParams.fs);
N = N - mod(N,wLength);
numFrames = floor(N/((1-aParams.wOverlap)*wLength))-1;

Fsinusoid = 0*repmat(cos(2*pi*8*[0:numFrames-1]/numFrames*dur), size(Fbeg, 1), 1);
F = interp1([1 numFrames]', [Fbeg Fend]', 1:numFrames)'; if length(Fbeg) == 1, F=F'; end;
F = F + Fsinusoid + 10*randn(size(F));

Fbwsinusoid = 0*repmat(cos(2*pi*10*[0:numFrames-1]/numFrames*dur), size(Fbwbeg, 1), 1);
Fbw = interp1([1 numFrames]', [Fbwbeg Fbwend]', 1:numFrames)'; if length(Fbwbeg) == 1, Fbw=Fbw'; end;
Fbw = Fbw + Fbwsinusoid + 10*randn(size(Fbw));

Zsinusoid = 0*repmat(cos(2*pi*7*[0:numFrames-1]/numFrames*dur), size(Zbeg, 1), 1);
Z = interp1([1 numFrames]', [Zbeg Zend]', 1:numFrames)'; if length(Zbeg) == 1, Z=Z'; end;
Z = Z + Zsinusoid + 10*randn(size(Z));

Zbwsinusoid = 0*repmat(cos(2*pi*6*[0:numFrames-1]/numFrames*dur), size(Zbwbeg, 1), 1);
Zbw = interp1([1 numFrames]', [Zbwbeg Zbwend]', 1:numFrames)'; if length(Zbwbeg) == 1, Zbw=Zbw'; end;
Zbw = Zbw + Zbwsinusoid + 10*randn(size(Zbw));

if ~isempty(F) && ~isempty(Z)
    if trackBW
        x0 = [F(:, 1); Fbw(:, 1); Z(:, 1); Zbw(:, 1)]+offset;
    else
        x0 = [F(:, 1); Z(:, 1)]+offset;
    end
end

if ~isempty(F) && isempty(Z)
    if trackBW
        x0 = [F(:, 1); Fbw(:, 1)]+offset;
    else
        x0 = F(:, 1)+offset;
    end
end

if isempty(F) && ~isempty(Z)
    if trackBW
        x0 = [Z(:, 1); Zbw(:, 1)]+offset;
    else
        x0 = Z(:, 1)+offset;
    end
end

%%
rmse = zeros(numTrials, 1);
x_est = cell(numTrials, 1);
x_errVar = cell(numTrials, 1);
x = cell(numTrials, 1);

for jj = 1:numTrials
    disp(['Processing Trial #', num2str(jj), '...'])
    [rmse(jj), x_est{jj}, x_errVar{jj}, x{jj}] = runSynth_OLA(F, Fbw, Z, Zbw, N, snr_dB, ...
            cepOrder, fs, trackBW, plot_flag, algFlag, x0, aParams, sParams);
end

if trackBW
    trueState = [F; Fbw; Z; Zbw];
else
    trueState = [F; Z]; % 
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
aParams.wLength = floor(aParams.wLengthMS/1000*fs);

plotSpecTracks2(x{trial}, x_est{trial}, aParams, nZ, trackBW);
firstline = get(get(gca,'Title'),'String');
title({['Trial #', num2str(trial)]; firstline})
axis tight