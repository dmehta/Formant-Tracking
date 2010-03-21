function plotRMSE(trueState,estTracks,titleCell)
%function plotRMSE(trueState,estTracks,titleCell)
% Plots RMSE for formant tracks vs. ground truth

% TODO: Make the loop dependent on state size

% Number of track estimates
numEst = size(estTracks,3);
numObser = size(trueState,2);
numForms = size(trueState,1);

%Titles for plot legends
ii = 2:(numEst+1);
S = titleCell(1,ii);

% Perfrom MSE Calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2); clf,
for ff = 1:numForms
    for ii = 1:numEst
        rmse(ii,ff) = norm((estTracks(ff,:,ii)-trueState(ff,:)))/sqrt(numObser);
        relRmse(ii,ff) = (rmse(ii,ff)/norm(trueState(ff,:)))*sqrt(numObser);
    end
end

subplot(1,2,1)
for ii = 1:numEst
    plot(rmse(ii,:),[char(titleCell(2,ii+1)) 'o']); 
    hold on; grid on; xlim([0.5 numForms+0.5])
end
title('Root MSE')
ylabel('RMSE (Hz)')
xlabel('Formant Number')

subplot(1,2,2)
for ii = 1:numEst
    plot(relRmse(ii,:),[char(titleCell(2,ii+1)) 'o']); 
    hold on; grid on; xlim([0.5 numForms+0.5])
end
title('Relative Root MSE')
ylabel('RMSE Ratio')
xlabel('Formant Number')
legend(S)
