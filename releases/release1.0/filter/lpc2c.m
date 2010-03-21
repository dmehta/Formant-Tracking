function [C] = lpc2c(arCoeffs, cepOrder)
%function [C] = lpc2c(arCoeffs, cepOrder)
% Using AR coefficients computes a specified number of the coefficients
% of the corresponding real cepstrum
%
% INPUT:
%    1. arCoeffs: p x numObs -- the ar coefficients for each frame of data
%    2. cepOrder: Number of cepstral coefficients to compute

% Author: Dan Spendley 3/18/2007
% Modified:  Dan Rudoy 3/18/2007

% Get LPC coeffs size
[p, numObser] = size(arCoeffs);

% Put together C matrix structure for use in loop
C = zeros(cepOrder,numObser);

% Compute cepstrum for every frame
for ii = 1:numObser
    C(:,ii) = lp2cf(arCoeffs(:,ii),cepOrder);
end

% Auxililary function implementing the necessary recursion
function [cepCoeffs] = lp2cf(arCoeffs, cepOrder)

p = length(arCoeffs);  % Number AR coefficients given
C = zeros(cepOrder,1); % Number Cepstral coefficients to produce

% Initial case in recursion
C(1) = arCoeffs(1);

% Recurse: index lower than max AR coeffs
for jj = 2:p
    a_seq = arCoeffs(jj:-1:1);
    c_seq = C(1:1:jj-1);
    aaInds = 1:1:length(c_seq);
    tmp = aaInds'.*c_seq(aaInds).*a_seq(aaInds+1)/jj;
    C(jj) = a_seq(1) + sum(tmp);
end

% Recurse: index higher than max AR coeffs
for jj = p+1:cepOrder
    a_seq = arCoeffs(p:-1:1)';
    c_seq = C(jj-p:1:jj-1);
    aaInds = 1:1:length(c_seq);
    C(jj) = sum(aaInds.*c_seq(aaInds)'.*a_seq(aaInds))/p;
end

% Set output
cepCoeffs = C;
