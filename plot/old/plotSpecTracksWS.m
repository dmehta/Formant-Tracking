%% plot center frequency and bandwidth tracks on spectrogram assuming
%% Wavesurfer-type input (only formants and bandwidths, no antiformants)

function [] = plotSpecTracksWS(audio, f, bw, aParams, trackBW)

if ~exist('trackBW', 'var')
    trackBW = 1;
end

fs = aParams.fs;
wLength = aParams.wLength;
wOverlap = aParams.wOverlap;
peCoeff = 0.97;

audio2 = normalize(filter([1 -peCoeff], 1, audio));

% Spectrogram Code
if wOverlap == 0
    winlen = wLength;
    winlen = winlen + (mod(winlen,2)~=0); % force even
    winoverlap = wOverlap;
else
    winlen = floor(fs/1000*4);  % 8 ms to make plots clear
    winlen = winlen + (mod(winlen,2)~=0); % force even
    winoverlap = 0.5*floor(winlen);
end
lineW = 1.5;
Fcolor = 'b';

%%
[B, freq, t] = specgram(audio2,512,fs,hamming(winlen),winoverlap);
imagesc(t,freq,20*log10(abs(B))), axis xy
% [S,F,T,P] = spectrogram(audio2,hamming(winlen),winoverlap,256,fs);
% surf(T,F,10*log10(abs(P)), 'EdgeColor', 'None');
% view(0, 90)
% grid off
axis tight, box off
colormap(flipud(pink(256)))
% colorbar EastOutside
rangemax = max(max(20*log10(abs(B))));
set(gca, 'CLim', [rangemax-60 rangemax])
xlabel('Time (s)'), ylabel('Frequency (Hz)'), title('Tracks of formants (blue) and anti-formants (red)')
hold on;

%%
if wOverlap == 0, xScale = wLength/fs; else xScale = 1/fs*wLength*wOverlap; end
xStart = xScale/2;
xInd = xStart:xScale:xStart+xScale*(size(f,2)-1);
len = size(f,2);

xdata = xInd(1:len);
numF = size(f,1);
for ii = 1:numF
    trackF = f(ii,:);
    if trackBW
        BW = bw(ii,:);
        L = trackF-BW/2;
        U = trackF+BW/2;
        fill([xdata xdata(end:-1:1)], [L U(end:-1:1)], Fcolor, 'EdgeColor', 'none')
    end
    plot(xdata,trackF, Fcolor, 'LineWidth', lineW);
end

format_plot