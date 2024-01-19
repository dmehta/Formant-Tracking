% addpath(genpath('../')); % Set paths
 
fs = 16e3; % Hz
cepOrder = 15;
 
%% set center frequencies and bandwidths
F = 500; Fbw = 100;
F=F'; Fbw = Fbw';
 
Z = 2000; Zbw = 80;
Z=Z'; Zbw = Zbw';
 
%% route 1 vs route 2
[num, denom] = fb2tf(F, Fbw, Z, Zbw, fs);
figure, freqz(num, denom, 512, fs)
C = lpc2cz(-denom(2:end)',-num(2:end)',cepOrder); % for speech to LPCC coefficients
 
% route 2
C2 = fb2cpz(F, Fbw, Z, Zbw, cepOrder, fs); % for estimation equation
 
[C C2']
sum(C-C2')