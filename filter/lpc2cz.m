function C = lpc2cz(arCoeffs, maCoeffs, cepOrder)

% Using ARMA coefficients computes a specified number of the coefficients
% of the corresponding real cepstrum. NB: Do not input leading AR or MA
% coefficient. ALSO, negate output of lpc() or arcov() or fb2tf().
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

[p, Np] = size(arCoeffs); % Get LPC coeffs size
[q, Nz] = size(maCoeffs); % Get LPC coeffs size

C = zeros(cepOrder, max(Np, Nz));   % Initialize LPCC matrix
C1 = C; C2 = C;

% Compute cepstrum for every frame
for i = 1:max(Np, Nz)
    if p ~= 0 % if AR coefficients input
        C1(:,i) = ar2cp(arCoeffs(:,i),cepOrder);
    end
    
    if q ~= 0 % if MA coefficients input
        C2(:,i) = ar2cp(maCoeffs(:,i),cepOrder);
    end
        
    C(:,i) = C1(:,i) - C2(:,i);
end