function [mn md mse] = mseAvg(matFName)
%function [mn md mse] = mseAvg(matFName)
% [mn md mse] = mseAvg'(rmse_3500_cj1_qt15_uc1_se0');

% Load experiement results
load(matFName);

% Summarize mean stats
mean_rmse = mean(rmse,3); 
mean_rmseAll = mean(rmseAll,3);

mean_delta = diff(mean_rmse')'; 
mean_deltaAll = diff(mean_rmseAll')';

mn = [mean_delta mean_deltaAll];

avg_mean = mean(mn,1);

% Summarize median stats
med_rmse = median(rmse,3); 
med_rmseAll = median(rmseAll,3);

med_delta = diff(med_rmse')'; 
med_deltaAll = diff(med_rmseAll')';

md = [med_delta med_deltaAll];

avg_med = mean(md,1);

mse = [avg_mean; avg_med];