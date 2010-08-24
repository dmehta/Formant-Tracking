% x_est is a cell array, where x_est{ii} is one trials of a 2D matrix
% of numStates x numFrames
%
% Plot estimated formant tracks for poles/zeros and bandwidths with
% estimator and 95 % confidence interval
% 
% Author: Patrick J. Wolfe, Daniel Rudoy, Daryush Mehta
%
% Created: 05/12/2010
% Modified: 05/12/2010, 06/01/2010 (trackBW and zeroes), 07/14/2010
% (estimator variance)

function plotStateTracksFZ_EstVar(x_est,x_errVar,nP,trackBW)

% Number of track estimates
numStates = size(x_est,1);
numObs = size(x_est,2);
xdata = 1:numObs;

if trackBW
    nZ = numStates/2-nP;
else
    nZ = numStates-nP;
end

m = ceil(sqrt(numStates));
n = ceil(numStates/m);
yrangemax = 0;
figure; clf,
for ff = 1:numStates
    subplot(m,n,ff)
    %grid on;
    box off
    hold on;
    means = x_est(ff,:);
    plot(means)
    variances = squeeze(x_errVar(ff,ff,:))';
    low = means + tinv(0.5-95/100/2, length(means)-1)*sqrt(variances)/sqrt(length(means));
    high = means + tinv(0.5+95/100/2, length(means)-1)*sqrt(variances)/sqrt(length(means));
    fill([xdata xdata(end:-1:1)], [low high(end:-1:1)], [0.9 0.9 0.9], 'EdgeColor', 'none')
    plot(xdata, means, 'LineWidth', 1)
    yrange = get(gca, 'YLim');
    yrangemax = max(yrangemax, yrange(2)-yrange(1));
    
    if ~nZ
        if ~trackBW
            title(['Resonance ' int2str(ff)]);
        else
            if ff > nP
                title(['Bandwidth ' int2str(ff-nP)]);
            else
                title(['Resonance ' int2str(ff)]);
            end            
        end
    else% zeros also
        if ~trackBW
            if ff > nP
                title(['Anti-resonance ' int2str(ff-nP)]);
            else
                title(['Resonance ' int2str(ff)]);
            end
        else
            if ff <= nP
                title(['Resonance ' int2str(ff)]);
            end

            if ff > nP && ff <= 2*nP
                title(['Resonance BW ' int2str(ff-nP)]);
            end

            if ff > 2*nP && ff <= (2*nP + nZ)
                title(['Anti-resonance ' int2str(ff-2*nP)]);
            end
            
            if ff > (2*nP + nZ) && ff <= (2*nP + 2*nZ)
                title(['Anti-resonance BW ' int2str(ff-2*nP-nZ)]);
            end
        end
    end
    
    if ff==1
        xlabel('Time Block');
        ylabel('Frequency (Hz)');
    end
end

for ff = 1:numStates
    subplot(m,n,ff)
    yrange = get(gca, 'YLim');
    set(gca, 'YLim', [yrange(1) yrange(1)+yrangemax])
end