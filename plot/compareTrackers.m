function compareTrackers(E, method)
% function compareTrackers(E)

if(strcmp(method,'EKF'))
    methodIndex = 1;
else if (strcmp(method,'EKS'))
        methodIndex = 2;
    else
        methodIndex = 3;
    end
end

E.estTracks = E.estTracks(:,:,methodIndex);

% Errorbar Plot Code
figure(3),clf; plot(E.trueState','-.')
figure(3),hold on;
plot(E.wsState','-')
plot(E.estTracks',':')

% plotEBars(E)
% plotCoast(E)
grid on;
title('Wavesurfer (solid) vs EKF (dotted)')

% Spectrogram Code
winlen = floor(E.fs/1000*8);  % 8 ms to make plots clear
winlen = winlen + (mod(winlen,2)~=0); % force even
winoverlap = winlen/2; % 50pct overlap

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4); clf, subplot(4,1,1)
specgram(E.wav,winlen,E.fs,hamming(winlen),winoverlap);
hold on;

% Overlay Kalman Resonant Tracks onto Spectrogram
plotSpecTracks(E,E.estTracks,'k')
plotSpecCoast(E)
title('Tracker comparison: EKF (black) vs. VTR (white)')

% Overlay Wavesurfer Resonant Tracks onto Spectrogram
figure(4); subplot(4,1,2)
specgram(E.wav,winlen,E.fs,hamming(winlen),winoverlap);
hold on;
plotSpecTracks(E,E.wsState,'b')
title('Tracker comparison: Wavesurfer (blue) vs. VTR (white)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Repeat, swap order
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5); clf, subplot(4,1,2)
specgram(E.wav,winlen,E.fs,hamming(winlen),winoverlap);
hold on;

% Overlay Kalman Resonant Tracks onto Spectrogram
plotSpecTracks(E,E.estTracks,'k')
plotSpecCoast(E)
title('Tracker comparison: EKF (black) vs. VTR (white)')

% Overlay Wavesurfer Resonant Tracks onto Spectrogram
figure(5); subplot(4,1,1),
specgram(E.wav,winlen,E.fs,hamming(winlen),winoverlap);
hold on;
plotSpecTracks(E,E.wsState,'b')
title('Tracker comparison: Wavesurfer (blue) vs. VTR (white)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

debug = 0;
