function C = lpc2cz(arCoeffs, maCoeffs, cepOrder)

% Using ARMA coefficients computes a specified number of the coefficients
% of the corresponding real cepstrum
%
% INPUT:
%    arCoeffs: p x numObs -- the ar coefficients for each frame of data
%    maCoeffs: q x numObs -- the ma coefficients for each frame of data
%    cepOrder: Number of cepstral coefficients to compute
% 
% OUTPUT:
%    C - LPC Cepstral coefficients

% Author: Daryush
% Created:  4/23/2010
% Modified: 4/23/2010

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[p, N] = size(arCoeffs); % Get LPC coeffs size
[q, N] = size(maCoeffs); % Get LPC coeffs size

C = zeros(cepOrder, N);   % Initialize LPCC matrix

% Compute cepstrum for every frame
for i = 1:N
    C(:,i) = ar2cp(arCoeffs(:,i),cepOrder);
    C(:,i) = C(:,i) - ar2cp(maCoeffs(:,i),cepOrder); % update to include zeros
end