function [N, D] = fb2tf(F, Fbw, Z, Zbw, fs)

% Computes coefficients for difference equation given center frequencies
% and bandwidths of resonances and anti-resonances
%
% INPUT:
%    F:     center frequencies of the resonances, in Hz
%    Fbw:   corresponding bandwidths of the resonances, in Hz
%    Z:     center frequencies of the anti-resonances, in Hz
%    Zbw:   corresponding bandwidths of the anti-resonances, in Hz
%    fs:    sampling rate, in Hz
% 
% OUTPUT:
%    N:     coefficients of transfer function numerator, N(1) forced to 1
%    D:     coefficients of transfer function denominator, D(1) forced to 1
% 
% Reference: 
% Klatt, Dennis (1980). Software for a cascade/parallel formant synthesizer, JASA, Vol. 67, No.
% 3.
% Gold, B. & Rabiner, L. (1968). Analysis of digital and analog formant
% synthesizers. IEEE Transactions on Audio and Electroacoustics, 16, 81-94.
% Eq. (4) and two paragraphs below Eq. (4).
% 
% Author: Daryush
% Created:  4/23/2010
% Modified: 4/23/2010

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a=zeros(1, length(F));
b=zeros(1, length(F));
c=zeros(1, length(F));
N = 1;
D = 1;

if ~isempty(F)
    % resonances (formants)
    for i = 1:length(F)
        
        % based on Klatt, Eq. (2) but taking into account Matlab transfer
        % function with positive coefficients in denominator
        b(i) = -2 * cos(2*pi*F(i)/fs) * exp(-pi*Fbw(i)/fs);
        c(i) = exp(-2*pi*Fbw(i)/fs);
        a(i) = 1 + b(i) + c(i);

        N = N*a(i);
        D = conv(D,[1 b(i) c(i)]);
    end
end

if ~isempty(Z)
    % anti-resonance (zeros)
    for i = 1:length(Z)
        
        b(i) = -2 * cos(2*pi*Z(i)/fs) * exp(-pi*Zbw(i)/fs);
        c(i) = exp(-2*pi*Zbw(i)/fs);
        a(i) = 1 + b(i) + c(i);
        
        N = conv(N,[1 b(i) c(i)]);
        D = D*a(i);
    end
end

N = N/N(1);
D = D/D(1);