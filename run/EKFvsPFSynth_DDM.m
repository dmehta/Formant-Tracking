function [] = EKFvsPFSynth_DDM()

doRuns  = 0;
doPlots = 1;

% don't forget to set algorithm flags
% Decision Flags for Tracking Processes, and parameters
% algFlag = [1 0 0 1 0]; % Select 1 to run, 0 not to
% EKF = 1; EKS = 2; EKS_EM = 3; PF = 4; RBPF = 5;

if(doRuns)

    numParticles = [10 50 100 200 500 1000];
    numIterations = 25;

    for numP = 1:length(numParticles)
        for n = 1:numIterations
            rand('state', n+100); randn('state', n+100); % Set seeds
            
            % runSynth(testMethod, snr_dB, cepOrder, fs, numFormants, trackBW, numParticles, doPlots, varargin)
            rmse = runSynth('Synth',15, 15, 10000, 4, 0, numParticles(numP), 0, 100, 30^2);
            meanError(n,numP,:) = mean(rmse);
        end
    end
        
    %save '../results/EKFvsPF.mat'; % Dan's
    %save '../results/EKFvsPF2.mat'; % DDM trying to replicate
    %save '../results/EKFvsPF3.mat'; % DDM trying to replicate, success! fix seeds on generation
    save '../results/EKFvsPF4.mat'; % set different seed for each iteration but fixed within change in numParticles: JASA resubmission
end

if(doPlots)
    %load '../results/EKFvsPF.mat'; % Dan's
    %load '../results/EKFvsPF2.mat'; % DDM trying to replicate
    %load '../results/EKFvsPF3.mat'; % DDM trying to replicate, success! fix seeds on generation
    load '../results/EKFvsPF4.mat'; % set different seed for each iteration but fixed within change in numParticles: JASA resubmission
    
    per = 95; % 95% confidence interval
    error1 = meanError(:,:,1);
    error2 = meanError(:,:,2);
    
    [low1 high1 ave1] = findCI(error1, per);
    [low2 high2 ave2] = findCI(error2, per);

    figure; hold on
    plot(numParticles, ave1, 'b')
    plot(numParticles, ave2, 'r', numParticles, low2, '--k', numParticles, high2, '--k')
    axis tight;
    xlim([0 1000]);
    ylim([25 50]);
    legend('EKF', 'PF');
    xlabel('Number of particles');
    ylabel('RMSE (Hz)')
end