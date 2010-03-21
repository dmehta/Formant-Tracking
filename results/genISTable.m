function genISTable(matFName)
%function genISTable(matFName) 
% genISTable('rmse_3500_cj1_qt0_uc1_se100');

%mse = mseAvg('rmse_3500_cj1_qt15_uc1_se0');
%mse = mseAvg('rmse_3500_cj1_qt15_uc1_se50');
%mse = mseAvg('rmse_3500_cj1_qt15_uc1_se100');

% Load best-fit dataset
load(matFName)

% Average over utterance (all of VTR database)
mean_rmse = mean(rmse,3); 
mean_rmseAll = mean(rmseAll,3);

ws = mean_rmse(:,3);
wsAll = mean_rmseAll(:,3);

% Compute reduction in MSE, ignoring smoother
All = diff(mean_rmseAll(:,[1 3])')';
Active = diff(mean_rmse(:,[1 3])')';

% Table: All/Active vs formant number
% Percentage reduction relative to wavesurfer
%num2str([All Active; sum([All Active])])
num2str([All All./wsAll*100 Active Active./wsAll*100],4)

if 0
    % Compute RMSE mean and variance over all utterances
    mean_rmse = mean(rmse,3);
    [m n] = size(mean_rmse);
    var_rmse = zeros(m,n);
    for i = 1:m
        for j = 1:n
            var_rmse(i,j) = var(rmse(i,j,:));
        end
    end
    diff(mean_rmse')';

    % Plot RMSE per utterance
    k = 1;
    figure(7); 
    for i= 1:m
        for j= 1:n
            subplot(m,n,k);
            %hold on; 
            plot(squeeze(rmse(i,j,:)));
            axis tight
            grid
            k = k+1;
        end
    end
end
% 
% % 176 Looks good:
% E = MyTrackExp('Timit',[],15,16000,3,[],[],176,1,0.15,1,1);
% 
% % 12 Looks on par
% E = MyTrackExp('Timit',[],15,16000,3,[],[],12,1,0.15,1,1);
% 
% % 236 Looks bad
% E = MyTrackExp('Timit',[],15,16000,3,[],[],236,1,0.15,1,1);
% 
% % 428 Randomly chosen
% E = MyTrackExp('Timit',[],15,16000,3,[],[],428,1,0.15,1,1);
% 
% % 274 Randomly chosen
% E = MyTrackExp('Timit',[],15,16000,3,[],[],274,1,0.15,1,1);

keyboard;
