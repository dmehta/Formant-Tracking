%
%
% Author: Patrick J. Wolfe, Daniel Rudoy, Daryush Mehta
%
% Created: 03/13/2007 
% Modified: 05/09/2010, 08/23/2010 (track bandwidths and zeros, but no CIs)

function plotStateTracksFZ(trueState,estTracks,titleCell,nP,trackBW)
% Plot estimated formant tracks for poles/zeros vs. ground truth

% Number of track estimates
numEst = size(estTracks,3);
numStates = size(trueState,1);

if trackBW
    nZ = numStates/2-nP;
else
    nZ = numStates-nP;
end


%Titles for plot legends
ii = 1:(numEst+1);
S = titleCell(1,ii);

m = ceil(sqrt(numStates));
n = ceil(numStates/m);
yrangemax = 0;
figure; clf,
for ff = 1:numStates
    subplot(m,n,ff)
    %grid on;
    box off
    hold on;
    plot(trueState(ff,:), char(titleCell(2,1)));
    for tt = 1:numEst
        plot(estTracks(ff,:,tt), char(titleCell(2,tt+1)));
    end
    yrange = get(gca, 'YLim');
    yrangemax = max(yrangemax, yrange(2)-yrange(1));
    
    if ~nZ
        if ~trackBW
            title(['Formant ' int2str(ff)]);
        else
            if ff > nP
                title(['Formant BW ' int2str(ff-nP)]);
            else
                title(['Formant ' int2str(ff)]);
            end            
        end
    else% zeros also
        if ~trackBW
            if ff > nP
                title(['Anti-formant ' int2str(ff-nP)]);
            else
                title(['Formant ' int2str(ff)]);
            end
        else
            if ff <= nP
                title(['Formant ' int2str(ff)]);
            end

            if ff > nP && ff <= 2*nP
                title(['Formant BW ' int2str(ff-nP)]);
            end

            if ff > 2*nP && ff <= (2*nP + nZ)
                title(['Anti-formant ' int2str(ff-2*nP)]);
            end
            
            if ff > (2*nP + nZ) && ff <= (2*nP + 2*nZ)
                title(['Anti-formant BW ' int2str(ff-2*nP-nZ)]);
            end
        end
    end
    
    if ff==1
        legend(S{1},S{2});
        xlabel('Frame number');
        ylabel('Frequency (Hz)');
    end
    
    format_plot
end

for ff = 1:numStates
    subplot(m,n,ff)
    yrange = get(gca, 'YLim');
    set(gca, 'YLim', [yrange(1) yrange(1)+yrangemax])
end