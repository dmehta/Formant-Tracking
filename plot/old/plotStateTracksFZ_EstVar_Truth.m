% x_est is a cell array, where x_est{ii} is one trials of a 2D matrix
% of numStates x numFrames
%
% Plot estimated formant tracks for poles/zeros and bandwidths with
% estimator and posterior interval as plus/minus one standard deviation
% around point estimates. With truth tracks as well (this function is good
% when using generative models)
% 
% Author: Patrick J. Wolfe, Daniel Rudoy, Daryush Mehta
%
% Created: 05/12/2010
% Modified: 05/12/2010, 06/01/2010 (trackBW and zeroes), 07/14/2010
% (estimator variance as 95 % interval), 08/26/2010 (estimator variance as stdev, truth)

function plotStateTracksFZ_EstVar_Truth(trueState,x_est,x_errVar,titleCell,nP,trackBW)

% Number of track estimates
numEst = size(x_est,3);
numStates = size(trueState,1);
numObs = size(x_est,2);
xdata = 1:numObs;

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
    means = x_est(ff,:);
    variances = squeeze(x_errVar(ff,ff,:))';
    low = means - sqrt(variances);
    high = means + sqrt(variances);

    % find indices for NaN boundaries
    notnan = diff(find(~isnan(low)));
    [r, c] = find(notnan > 1);
    for ii = 1:length(c)
        itp = c(ii);
        fill([xdata(1:itp) xdata(itp:-1:1)], [low(1:itp) high(itp:-1:1)], [0.9 0.9 0.9], 'EdgeColor', 'none')        
        if ii == length(c)
            itp2 = c(ii)+notnan(c(ii));
            fill([xdata(itp2:end) xdata(end:-1:itp2)], [low(itp2:end) high(end:-1:itp2)], [0.9 0.9 0.9], 'EdgeColor', 'none')
        end
    end
    
    fill([xdata xdata(end:-1:1)], [low high(end:-1:1)], [0.9 0.9 0.9], 'EdgeColor', 'none')
    plot(xdata, means, char(titleCell(2,2)), 'LineWidth', 1)
    plot(trueState(ff,:), char(titleCell(2,1)));
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
        legend('+/- 1 SD', S{2},S{1})
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