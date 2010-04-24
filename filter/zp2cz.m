function C = zp2cz(pol, zer, cepOrder, fs)

% Using zeros and poles, computes a specified number of the coefficients
% of the corresponding real cepstrum
%
% INPUT:
%    p: p x numObs -- the poles for each frame of data
%    z: q x numObs -- the zeros for each frame of data
%    cepOrder: Number of cepstral coefficients to compute
%    fs: sampling rate, in Hz
% 
% OUTPUT:
%    C - LPC Cepstral coefficients

% Author: Daryush
% Created:  4/23/2010
% Modified: 4/23/2010

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[p, N] = size(pol);
[q, N] = size(zer);

C = zeros(cepOrder, N);   % Initialize LPCC matrix

% Compute cepstrum for every frame
for i = 1:N
    for n = 1:cepOrder
        for ii = 1:p
            C1(:,ii) = pol(ii).^n;
        end
        
        for jj = 1:q
            C2(:,jj) = zer(jj).^n;
        end
        
        C(:,i) = 1/n*(C1+C2);
    end
end