%% plot center frequency and bandwidth tracks on spectrogram assuming
%% Wavesurfer-type input (only formants and bandwidths, no antiformants)

function [] = plotSpecTracksWS(audio, f, bw, aParams)

fs = aParams.fs;
wLength = aParams.wLength;
wOverlap = aParams.wOverlap;

% Spectrogram Code
winlen = floor(fs/1000*8);  % 8 ms to make plots clear
winlen = winlen + (mod(winlen,2)~=0); % force even
winoverlap = winlen/2; % 50pct overlap

figure;
specgram(audio,winlen,fs,hamming(winlen),winoverlap);
axis xy
colormap(flipud(pink(256)))
axis tight, box off
xlabel('Time (s)'), ylabel('Frequency (Hz)'), title('Tracks of formants (yellow) and antiformants (cyan)')
hold on;

xScale = 1/fs*wLength*wOverlap;
xStart = xScale/2;
xInd = xStart:xScale:xStart+xScale*(size(f, 2)-1);
len = size(f, 2);

xdata = xInd(1:len);
numF = size(f,1);
for ii = 1:numF
    trackF = f(ii,:);
    trackBW = bw(ii,:);
    L = trackF-trackBW/2;
    U = trackF+trackBW/2;
    fill([xdata xdata(end:-1:1)], [L U(end:-1:1)], 'y', 'EdgeColor', 'none', 'FaceAlpha', 0.3)
    plot(xdata,trackF,'y');
end 
set(findobj('Type','line'),'LineWidth',2)