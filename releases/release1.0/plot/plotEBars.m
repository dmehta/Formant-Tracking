function [] = plotEBars(E)
% Plots errorbars with std_dev at every point
% Created : 3/22/2007, Dan S.
% Last Updated : 3/22/2007

% Plot Errorbars
for i = 1:size(E.estTracks,1)
    errorbar(E.estTracks(i,:,1)',sqrt(E.estVar(1,1,:)),'.')
end