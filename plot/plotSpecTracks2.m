function [] = plotSpecTracks2(audio, tracks, aParams, varargin)

fs = aParams.fs;
wLength = aParams.wLength;
wOverlap = aParams.wOverlap;

% Spectrogram Code
winlen = floor(fs/1000*25);  % 8 ms to make plots clear...changed to 25 ms, 6/7/10
winlen = winlen + (mod(winlen,2)~=0); % force even
winoverlap = winlen/2; % 50pct overlap

figure;
specgram(audio,winlen,fs,hamming(winlen),winoverlap);
% spectrogram(audio, hamming(winlen), winoverlap, winlen, fs)
% B = specgram(audio,winlen,fs,hamming(winlen),winoverlap);
% [B,f,t] = specgram(audio,winlen,fs);
% imagesc(t, f, 20*log10(abs(B))), 
axis xy, colormap(flipud(pink(256)))
xlabel('Time (s)'), ylabel('Frequency (Hz)'), title('Tracks of formants (yellow) and antiformants (magenta)')
hold on;

xScale = 1/fs*wLength*wOverlap;
xStart = xScale/2;
xInd = xStart:xScale:xStart+xScale*(size(tracks(1,:,1),2)-1);
len = size(tracks,2);

if ~isempty(varargin)
    numAntiF = varargin{1};
    plot(xInd(1:len),tracks(1:end-numAntiF,1:len,1)','y*');
    
    if numAntiF
        plot(xInd(1:len),tracks(end-numAntiF+1:end,1:len,1)','m*');
    end
else
    plot(xInd(1:len),tracks(:,1:len,1)',['b*']);
end