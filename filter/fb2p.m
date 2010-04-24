function p = fb2p(f, BW, fs)

% given center frequency and bandwidth, calculate pole (or zero) pair

p = zeros(2*length(f),1);
for i = 1:length(f)
    p(2*i-1:2*i) = exp(-BW(i)/fs+1i*2*pi*BW(i)/fs*[1 -1]); % conjugate pair
end