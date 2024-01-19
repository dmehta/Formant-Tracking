function [] = EKFvsPFSynth()

doRuns  = 0;
doPlots = 1;

if(doRuns)

    numParticles = [10 50 100 200 500 1000];
    numIterations = 25;

    for n = 1:numIterations
        for numP = 1:length(numParticles)
            rmse = runSynth('Synth',15, 15, 10000, 4, 0, numParticles(numP), 100, 30^2);
            meanError(n,numP,:) = mean(rmse);
        end
    end

    save '../results/EKFvsPF.mat';
end

if(doPlots)
    load '../results/EKFvsPF.mat';
    close all;
    N = numIterations;
    errs   = squeeze(mean(meanError));  % Compute mean error
    varN   = squeeze(var(meanError))/N; % Compute variance
    stdErr = squeeze(std(meanError));   % Compute standard deviation
    
    pct = prctile(meanError, [2.5 50 97.5]);
    
    figure;
    plot(numParticles, errs(:,1));
    hold on;
    
    plot(numParticles, errs(:,2), 'r');
    plot(numParticles, pct(1,:,2),'--k');
    plot(numParticles, pct(3,:,2),'--k');
    
    axis tight;
    xlim([-25 1100]);
    ylim([23 50]);
    legend('EKF', 'PF');
    xlabel('Number of particles');
    ylabel('RMSE (Hz)')
    fmakep5(2,4);
end