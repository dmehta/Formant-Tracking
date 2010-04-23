function cepCoeffs = ar2cp(arCoeffs, cepOrder)

% Using AR coefficients computes a specified number of the coefficients
% of the corresponding real cepstrum
%
% INPUT:
%    arCoeffs: p x numObs -- the ar coefficients for each frame of data
%    cepOrder: Number of cepstral coefficients to compute
% OUTPUT:
%    C - LPC Cepstral coefficients

p = length(arCoeffs);  % Number AR coefficients given
C = zeros(cepOrder,1); % Number Cepstral coefficients to produce
C(1) = arCoeffs(1);    % Initial case in recursion

% Case 1: index lower than max num AR coeffs
for j = 2:p
    a_seq = arCoeffs(j:-1:1);
    aaInds = 1:1:j-1;
    c_seq = C(aaInds);
    C(j) = a_seq(1) + sum(aaInds'.*c_seq(aaInds).*a_seq(aaInds+1))/j;
end

% Case 2: index higher than max num AR coeffs
for j = p+1:cepOrder
    a_seq = arCoeffs(p:-1:1)';
    c_seq = C(j-p:1:j-1);
    aaInds = 1:1:length(c_seq);
    C(j) = sum(aaInds.*c_seq(aaInds)'.*a_seq(aaInds))/p;
end

cepCoeffs = C; % Set output