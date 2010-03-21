function [] = plotSpecTracks(E,tracks,color)
%function [] = plotSpecTracks(E,tracks,color)
%
% Function overlays estimated formant tracks onto spectrogram 
% fs, wLength, and wOverlap are neccesary to convert frame index to time
% index as plotted by spectrogram
% Created : 3/23/2007
% Last Updated : 3/23/2007

% Overlay Resonant Tracks onto Spectrogram
xScale = 1/E.fs*E.wLength*E.wOverlap;
xStart = xScale/2;
xInd = xStart:xScale:xStart+xScale*(size(tracks(1,:,1),2)-1);
len = min([size(tracks,2) ...
          size(E.trueState,2) ...
          size(E.estTracks,2) ...
          size(E.wsState,2)]);
plot(xInd(1:len),tracks(:,1:len,1)',[color '*']);
plot(xInd(1:len),E.trueState(:,1:len,1)','w*');