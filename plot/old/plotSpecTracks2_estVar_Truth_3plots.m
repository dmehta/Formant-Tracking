% same as plotSpecTracks2, except formant center frequency and UNCERTAINTY now plotted as shaded regions around center
% frequency tracks; no bandwidth information plotted
% varargin -- numAntiF, trackBW

function [] = plotSpecTracks2_estVar_Truth(audio, x_est, x_errVar, trueState, aParams, varargin)

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
    winlen = floor(fs/1000*4);  % 4 ms to make plots clear
    winlen = winlen + (mod(winlen,2)~=0); % force even
    winoverlap = winlen - 1;
end
lineW = 1.5;
Fcolor = 'b';
Fcolor_Truth = 'r';
Zcolor = 'g';
Zcolor_Truth = 'g--';

%% clean spectrogram
figure; subplot(311)
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
xlabel('Time (s)'), ylabel('Frequency (Hz)'), title('Tracks of formants (blue) and anti-formants (green), with truth formants (red)')
format_plot
hold on;

%%
subplot(312)
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
xlabel('Time (s)'), ylabel('Frequency (Hz)'), title('Tracks of formants (blue) and anti-formants (green), with truth formants (red)')
format_plot
hold on;

%%
if wOverlap == 0, xScale = wLength/fs; else xScale = 1/fs*wLength*wOverlap; end
xStart = xScale;
xInd = xStart:xScale:xStart+xScale*(size(x_est(1,:,1),2)-1);
len = size(x_est,2);

if length(varargin) > 1
    trackBW = varargin{2};
    if trackBW
        numStates = size(x_est,1)/2; % don't plot BW data
    else
        numStates = size(x_est,1);
    end
else
    numStates = size(x_est,1);
end

xdata = xInd(1:len);

switch length(varargin)
    case 0
        plot(xdata,x_est(:,1:len,1)'',Fcolor);
    case {1, 2}
        numAntiF = varargin{1};
        
        % plot uncertainties
        for ff = 1:numStates-numAntiF
            means = x_est(ff,1:len);
            variances = squeeze(x_errVar(ff,ff,:))';
            low = means - sqrt(variances);
            high = means + sqrt(variances);
            fill([xdata xdata(end:-1:1)], [low high(end:-1:1)], [0.9 0.9 0.9], 'EdgeColor', 'none')
        end
        
        if numAntiF
            for ff = numStates-numAntiF+1:numStates
                means = x_est(ff,1:len);
                variances = squeeze(x_errVar(ff,ff,:))';
                low = means - sqrt(variances);
                high = means + sqrt(variances);
                fill([xdata xdata(end:-1:1)], [low high(end:-1:1)], [0.9 0.9 0.9], 'EdgeColor', 'none')
            end
            
            for ff = numStates-numAntiF+1:numStates
                means = x_est(ff,1:len);
                means_Truth = trueState(ff,1:len);
                plot(xdata,means,Zcolor, 'LineWidth', lineW);
                plot(xdata,means_Truth,Zcolor_Truth, 'LineWidth', lineW);
            end
            
        end
        
        % plot means last so they're always on top
        for ff = 1:numStates-numAntiF
            means = x_est(ff,1:len);
            means_Truth = trueState(ff,1:len);
            plot(xdata,means,Fcolor, 'LineWidth', lineW);
            plot(xdata,means_Truth,Fcolor_Truth, 'LineWidth', lineW);
        end
        
    otherwise
        disp('Incorrect inputs into plotSpecTracks2_estVar_Truth().')
end

format_plot