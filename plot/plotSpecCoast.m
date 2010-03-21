function [] = plotSpecCoast(E)
% If coasting was used, function indicates at what samples coasting was
% implemented and plots them onto spectrogram
% Created : 3/23/2007, Dan S.
% Last Updated : 3/23/2007

% For converting indices to time series
xScale = 1/E.fs*E.wLength*E.wOverlap;
xStart = xScale/2;
xInd = xStart:xScale:xStart+xScale*(size(E.estTracks(1,:,1),2)-1);

for i = 1:size(E.formantInds,2)
    
    % Find non-zero indices of mask
    inds = find(E.formantInds(:,i) == 0);
    
    % Convert indices to time scale
    timeInds = (inds-1).*xScale+xStart;

    % Overlay data onto plots
    plot(timeInds,E.estTracks(i,(inds),1)','m*');
end