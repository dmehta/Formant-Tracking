function [C] = lpc2c(arCoeffs, cepOrder)

% Using AR coefficients computes a specified number of the coefficients
% of the corresponding real cepstrum
%
% INPUT:
%    arCoeffs: p x numObs -- the ar coefficients for each frame of data
%    cepOrder: Number of cepstral coefficients to compute
% OUTPUT:
%    C - LPC Cepstral coefficients

% Author: Daniel Rudoy
% Created:  3/18/2007
% Modified: 2/15/2010

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[p, N] = size(arCoeffs); % Get LPC coeffs size
C = zeros(cepOrder, N);   % Initialize LPCC matrix

% Compute cepstrum for every frame
for i = 1:N
    C(:,i) = ar2cp(arCoeffs(:,i),cepOrder);
end