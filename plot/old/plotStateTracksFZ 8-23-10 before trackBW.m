%
%
% Author: Patrick J. Wolfe, Daniel Rudoy, Daryush Mehta
%
% Created: 03/13/2007 
% Modified: 05/09/2010

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
yrangemax = 0;
figure; clf,
for ff = 1:numForms
    subplot(m,n,ff)
    plot(trueState(ff,:), char(titleCell(2,1)));
    %grid on;
    box off
    hold on;
    for tt = 1:numEst
        plot(estTracks(ff,:,tt), char(titleCell(2,tt+1)));
    end
    yrange = get(gca, 'YLim');
    yrangemax = max(yrangemax, yrange(2)-yrange(1));
    
    if ff > nP
        title(['Anti-resonance ' int2str(ff-nP)]);
    else
        title(['Resonance ' int2str(ff)]);
    end
    
    if ff==1
        legend(S{1},S{2});
        xlabel('Time Block');
        ylabel('Frequency (Hz)');
    end
end

for ff = 1:numForms
    subplot(m,n,ff)
    yrange = get(gca, 'YLim');
    set(gca, 'YLim', [yrange(1) yrange(1)+yrangemax])
end