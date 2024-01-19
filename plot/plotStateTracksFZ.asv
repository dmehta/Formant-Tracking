%
%
% Author: Patrick J. Wolfe, Daniel Rudoy
%
% Created: 03/13/2007 
% Modified: 12/13/2007

function plotStateTracksFZ(trueState,estTracks,titleCell,nP)
% Plot estimated formant tracks for poles/zeros vs. ground truth

% TODO: Make the loop dependent on state size

% Number of track estimates
numEst = size(estTracks,3);
numForms = size(trueState,1);

%Titles for plot legends
ii = 1:(numEst+1);
S = titleCell(1,ii);

% Assumes tracking all four formants
m = ceil(sqrt(numForms));
n = ceil(numForms/m);
figure; clf,
for ff = 1:numForms
    subplot(m,n,ff)
    plot(trueState(ff,:), char(titleCell(2,1)));
    grid on;
    hold on;
    for tt = 1:numEst
        plot(estTracks(ff,:,tt), char(titleCell(2,tt+1)));
    end
    title(['Resonance ' int2str(ff)]);
    legend(S{1},S{2});
    if ff==(n*(m-1)+1)
        xlabel('Time Block');
        ylabel('Frequency (Hz)');
    end
end
