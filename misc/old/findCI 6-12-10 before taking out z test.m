% Daryush Mehta
% findCI.m
% Created: May 17, 2010
% Modified: May 17, 2010
% 
% The purpose of this function is to compute a confidence interval for a given data set.
% 
% Cohen, B. H. (1996). Explaining Psychological Statistics. Pacific Grove, CA, Brooks/Cole
% Publishing Company, page 353.
% 
% function y = corrCI(x, per)
%
%   inputs:   x -- vector of data
%             per -- confidence interval percentage; e.g., per=95 for 95%
%                           confidence interval (p < 0.05, two-sided)
%
%   outputs:  low -- per% confidence interval, lower bounds
%             high -- per% confidence interval, higher bounds
%             ave -- mean estimates

function [low high ave] = findCI(x, per)

ydata = mean(x, 1);
yerror = std(x, 0, 1);

nn = length(ydata);
if nn > 40
    low = ydata + norminv(0.5-per/100/2, 0, 1)*yerror/sqrt(nn);
    high = ydata + norminv(0.5+per/100/2, 0, 1)*yerror/sqrt(nn);
else
    low = ydata + tinv(0.5-per/100/2, nn-1)*yerror/sqrt(nn);
    high = ydata + tinv(0.5+per/100/2, nn-1)*yerror/sqrt(nn);
end

ave = ydata;