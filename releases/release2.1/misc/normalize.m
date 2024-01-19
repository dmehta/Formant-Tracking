% Daryush Mehta
% normalize.m
% August 22, 2005
% The purpose of this function is to take in a vector and normalize it so
% that the signal is maximized to 1
%
% function [norm] = normalize(x)
%
%   inputs:   x -- input vector to be normalized
%
%   outputs:  norm -- norm/(max(norm))

function [norm] = normalize(x)

norm = 0.99*x./(max(abs(x)));