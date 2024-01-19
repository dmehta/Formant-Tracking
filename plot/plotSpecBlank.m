% blank spectrogram

function plotSpecBlank(audio, aParams)

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

%%
% figure;
[B, freq, t] = specgram(audio2,512,fs,hamming(winlen),winoverlap);
imagesc(t,freq,20*log10(abs(B))), axis xy
% [S,F,T,P] = spectrogram(audio2,hamming(winlen),winoverlap,256,fs);
% surf(T,F,10*log10(abs(P)), 'EdgeColor', 'None');
% view(0, 90)
% grid off
colormap(flipud(pink(256)))
% colorbar EastOutside
rangemax = max(max(20*log10(abs(B))));
set(gca, 'CLim', [rangemax-60 rangemax])
xlabel('Time (s)'), ylabel('Frequency (Hz)'), title('Tracks (blue) and truth (red)')
format_plot
hold on;