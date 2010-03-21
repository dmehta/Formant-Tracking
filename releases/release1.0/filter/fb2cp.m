function [cepCoeffs] = fb2cp(F, BW, numCoeffs, fs)
% function [cepCoeffs] = fb2cp(F, BW, numCoeffs, fs)
% Convert formant locations and bandwidths into a pre-specified number
% of cepstral coefficients, according to the standard mapping as given
% in the Interspeech 2007 paper 
%
% Author: Daniel Rudoy
% Created : 3/11/07
% Modified: 12/13/2007, 3/11/07 

% This code implements a vectorized version of the below:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i = 1:numCoeffs
%     for p = 1:length(F)
%         bw_term = (2/i)*exp(-pi*i*(BW(p)/fs));
%         C_int(i,p) = bw_term*cos(2*pi*i*(F(p))/fs);
%     end
%     cepCoeffs(i) = sum(C_int(i,:));
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lenF = length(F);
av = 1:numCoeffs;
bv = 1:lenF;
a = av(ones(lenF,1),:);
b = bv(ones(numCoeffs,1),:)';

% Compute cepstrum
bw = 2*exp(-pi*a.*(BW(b))/fs)./a;
C_int2 = bw.*cos(2*pi*a.*F(b)/fs);
cepCoeffs = sum(C_int2,1);