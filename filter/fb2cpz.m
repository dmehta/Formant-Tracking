function cepCoeffs = fb2cpz(F, Fbw, Z, Zbw, cepOrder, fs)

% Convert center frequencies/bandwidths of BOTH resonances and
% anti-resonances to a specified number of cepstral coefficients. Extension
% of equation in Interspeech 2007 paper [fb2cp.m].
%
% Author: Daryush Mehta
% Created : 4/23/10
% Modified: 4/23/10

% for kk = 1:length(cepOrder)
%     cepCoeffs = fb2cp(F(kk),Fbw(kk),cepOrder,fs); % LPCC coeffs due to poles
%     cepCoeffs = cepCoeffs - fb2cp(Z(kk),Zbw(kk),cepOrder,fs); % Adjust due to effect of zeros
% end
% 
% cepCoeffs=cepCoeffs';

C_int = zeros(cepOrder,length(F));
C_int2 = zeros(cepOrder,length(Z));
cepCoeffs = zeros(1,cepOrder);

for n = 1:cepOrder
    for i = 1:length(F)
        bw_term = exp(-pi*n*(Fbw(i)/fs));
        C_int(n,i) = bw_term*cos(2*pi*n*F(i)/fs);
    end
    
    for j = 1:length(Z)
        bw_term = exp(-pi*n*(Zbw(j)/fs));
        C_int2(n,j) = bw_term*cos(2*pi*n*Z(j)/fs);
    end

    cepCoeffs(n) = (2/n)*(sum(C_int(n,:))-sum(C_int2(n,:)));
end