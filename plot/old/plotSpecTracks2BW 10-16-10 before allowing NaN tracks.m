% same as plotSpecTracks2, except bandwidths now plotted as shaded regions around center
% frequency tracks
% varargin -- numAntiF, trackBW

function [] = plotSpecTracks2BW(audio, tracks, aParams, varargin)

fs = aParams.fs;
wLength = aParams.wLength;
wOverlap = aParams.wOverlap;

% Spectrogram Code
if wOverlap == 0
    winlen = wLength;
    winlen = winlen + (mod(winlen,2)~=0); % force even
    winoverlap = wOverlap;
else
    winlen = floor(fs/1000*8);  % 8 ms to make plots clear
    winlen = winlen + (mod(winlen,2)~=0); % force even
    winoverlap = round(0.95*winlen);
end
lineW = 1.5;
Fcolor = 'b';
Zcolor = 'r--';

%%
figure;
[B, freq, t] = specgram(audio,512,fs,hamming(winlen),winoverlap);
imagesc(t,freq,20*log10(abs(B))), axis xy
% [S,F,T,P] = spectrogram(audio,hamming(winlen),winoverlap,256,fs);
% surf(T,F,10*log10(abs(P)), 'EdgeColor', 'None');
% view(0, 90)
% grid off
colormap(flipud(gray(256)))
colorbar EastOutside
rangemax = max(max(20*log10(abs(B))));
set(gca, 'CLim', [rangemax-60 rangemax])
xlabel('Time (s)'), ylabel('Frequency (Hz)'), title('Tracks of formants (blue) and anti-formants (red)')
format_plot
hold on;

%%
if wOverlap == 0, xScale = wLength/fs; else xScale = 1/fs*wLength*wOverlap; end
xStart = xScale/2;
xInd = xStart:xScale:xStart+xScale*(size(tracks(1,:,1),2)-1);
len = size(tracks,2);

xdata = xInd(1:len);

switch length(varargin)
    case 0
        plot(xdata,tracks(:,1:len,1)'',Fcolor);
    case 1
        numAntiF = varargin{1};
        plot(xdata,tracks(1:end-numAntiF,1:len,1)'',Fcolor, 'LineWidth', lineW);
    
        if numAntiF
            plot(xdata,tracks(end-numAntiF+1:end,1:len,1)',Zcolor, 'LineWidth', lineW);
        end
    case 2
        numAntiF = varargin{1};
        trackBW = varargin{2};

        if trackBW
            numF = (size(tracks,1)-2*numAntiF)/2;
            
            if numAntiF
                for ii = 1:numF
                    trackF = tracks(ii,1:len,1);
                    trackBW = tracks(numF+ii,1:len,1);
                    L = trackF-trackBW/2;
                    U = trackF+trackBW/2;
                    fill([xdata xdata(end:-1:1)], [L U(end:-1:1)], Fcolor, 'EdgeColor', 'none')%, 'FaceAlpha', 0.3)
                    plot(xdata,trackF','k', 'LineWidth', lineW);
                end
                
                for ii = 1:numAntiF
                    trackF = tracks(2*numF+ii,1:len,1);
                    trackBW = tracks(2*numF+numAntiF+ii,1:len,1);
                    L = trackF-trackBW/2;
                    U = trackF+trackBW/2;
                    fill([xdata xdata(end:-1:1)], [L U(end:-1:1)], Zcolor, 'EdgeColor', 'none')%, 'FaceAlpha', 0.3)
                    plot(xdata,trackF,'k', 'LineWidth', lineW);
                end
                
            else
                for ii = 1:numF
                    trackF = tracks(ii,1:len,1);
                    trackBW = tracks(numF+ii,1:len,1);
                    L = trackF-trackBW/2;
                    U = trackF+trackBW/2;
                    fill([xdata xdata(end:-1:1)], [L U(end:-1:1)], Fcolor, 'EdgeColor', 'none')%, 'FaceAlpha', 0.3)
                    plot(xdata,trackF,'k', 'LineWidth', lineW);
                end
            end            
        else
            numF = size(tracks, 1)-numAntiF;

            if numAntiF
                for ii = 1:numF
                    trackF = tracks(ii,1:len,1);
                    plot(xdata,trackF','k', 'LineWidth', lineW);
                end
                
                for ii = 1:numAntiF
                    trackF = tracks(numF+ii,1:len,1);
                    plot(xdata,trackF,'k', 'LineWidth', lineW);
                end
                
            else
                for ii = 1:numF
                    trackF = tracks(ii,1:len,1);
                    plot(xdata,trackF','k', 'LineWidth', lineW);
                end
            end            
        end
    otherwise
        disp('Incorrect inputs into plotSpecTracks2().')
end

format_plot