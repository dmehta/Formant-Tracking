% /n/
dur = 1; % in s
Fbeg = [257 1491]'; Fend = [257 1491]';
Fbwbeg = [32 100]'; Fbwend = [32 100]';
Zbeg = [1223]'; Zend = [1223]';
Zbwbeg = [52]'; Zbwend = [52]';

runSynth_OLA_wrapper

audio = x{1,1};
wavwritedm(audio, fs, 16, 'n_synth.wav')

% /n a/
dur = .5; % in s
Fbeg = [257 1491]'; Fend = [850 1100]';
Fbwbeg = [32 100]'; Fbwend = [80 120]';
Zbeg = []'; Zend = []';
Zbwbeg = []'; Zbwend = []';

runSynth_OLA_wrapper

audio = x{1,1};
wavwritedm(audio, fs, 16, 'na_synth.wav')

% /a/
dur = 1; % in s
Fbeg = [850 1100]'; Fend = [850 1100]';
Fbwbeg = [80 120]'; Fbwend = [80 120]';
Zbeg = []'; Zend = []';
Zbwbeg = []'; Zbwend = []';

runSynth_OLA_wrapper

audio = x{1,1};
wavwritedm(audio, fs, 16, 'a_synth.wav')

% /a n/
dur = .5; % in s
Fbeg = [850 1100]'; Fend = [257 1491]';
Fbwbeg = [80 120]'; Fbwend = [32 100]';
Zbeg = []'; Zend = []';
Zbwbeg = []'; Zbwend = []';

runSynth_OLA_wrapper

audio = x{1,1};
wavwritedm(audio, fs, 16, 'an_synth.wav')

%%
n = wavread('n_synth.wav');
na = wavread('na_poles_synth.wav');
a = wavread('a_synth.wav');
an = wavread('an_poles_synth.wav');

%%
nan = zeros(76*500, 1);
nan(1:(19+1)*500) = n;
figure, plot(nan)
nan((19)*500+1:(19)*500+(9+1)*500) = nan((19)*500+1:(19)*500+(9+1)*500) + na;
hold on
plot(nan, 'r--')
nan((28)*500+1:(28)*500+(19+1)*500) = nan((28)*500+1:(28)*500+(19+1)*500) + a;
plot(nan, 'g--')
nan((47)*500+1:(47)*500+(9+1)*500) = nan((47)*500+1:(47)*500+(9+1)*500) + an;
plot(nan, 'm.-')
nan((56)*500+1:(56)*500+(19+1)*500) = nan((56)*500+1:(56)*500+(19+1)*500) + n;
plot(nan, 'y-')

wavwritedm(nan, fs, 16, 'nan_synth.wav')