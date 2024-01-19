%% fb2p
% given center frequency and bandwidth, calculate pole (or zero) pair

function p = fb2p(f, BW, fs)

for i = 1:length(f)
    p = exp(-BW(i)/fs+1i*2*pi*b(i)/fs*[1 -1]); % conjugate pair
end