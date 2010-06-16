clear
randn('state', 5)

%% parameters
dur = 2; % in s
snr_dB = 20;
cepOrder = 20;
fs = 10e3;
trackBW = 1;
plot_flag = 0;
algFlag = [0 1]; % Select 1 to run, 0 not to; [EKF EKS]
wType = 'hanning';  % window type
wLengthMS  = 100;    % Length of window (in milliseconds)
wOverlap = 0.5;     % Factor of overlap of window

%%
wLength = floor(wLengthMS/1000*fs);
wLength = wLength + (mod(wLength,2)~=0); % Force even
win = feval(wType,wLength);
N = floor(dur*fs);
N = N - mod(N,wLength);
numFrames = floor(N/((1-wOverlap)*wLength))-1;

%% trajectory
Fbeg = [500]';
Fend = [900]';
Fpert = 0*repmat(cos(2*pi*8*[0:numFrames-1]/numFrames*dur), size(Fbeg, 1), 1);
F = interp1([1 numFrames]', [Fbeg Fend]', 1:numFrames)'; if length(Fbeg) == 1, F=F'; end;
F = F + Fpert + 10*randn(size(F));

Fbwbeg = [100]';
Fbwend = [100]';
Fbwpert = 0*repmat(cos(2*pi*10*[0:numFrames-1]/numFrames*dur), size(Fbwbeg, 1), 1);
Fbw = interp1([1 numFrames]', [Fbwbeg Fbwend]', 1:numFrames)'; if length(Fbwbeg) == 1, Fbw=Fbw'; end;
Fbw = Fbw + Fbwpert + 10*randn(size(Fbw));

Zbeg = [900]';
Zend = [500]';
Zpert = 0*repmat(cos(2*pi*7*[0:numFrames-1]/numFrames*dur), size(Zbeg, 1), 1);
Z = interp1([1 numFrames]', [Zbeg Zend]', 1:numFrames)'; if length(Zbeg) == 1, Z=Z'; end;
Z = Z + Zpert + 10*randn(size(Z));

Zbwbeg = [100]';
Zbwend = [100]';
Zbwpert = 0*repmat(cos(2*pi*6*[0:numFrames-1]/numFrames*dur), size(Zbwbeg, 1), 1);
Zbw = interp1([1 numFrames]', [Zbwbeg Zbwend]', 1:numFrames)'; if length(Zbwbeg) == 1, Zbw=Zbw'; end;
Zbw = Zbw + Zbwpert + 10*randn(size(Zbw));

% F = []'; Fbw = []';
% Z = []'; Zbw = []';

%% initial state
if ~isempty(F) && ~isempty(Z)
    if trackBW
        x0 = [F(:, 1); Fbw(:, 1); Z(:, 1); Zbw(:, 1)]+100;
    else
        x0 = [F(:, 1); Z(:, 1)]+100;
    end
end

if ~isempty(F) && isempty(Z)
    if trackBW
        x0 = [F(:, 1); Fbw(:, 1)]+100;
    else
        x0 = F(:, 1)+100;
    end
end

if isempty(F) && ~isempty(Z)
    if trackBW
        x0 = [Z(:, 1); Zbw(:, 1)]+100;
    else
        x0 = Z(:, 1)+100;
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