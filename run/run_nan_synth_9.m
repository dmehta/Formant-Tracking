clear all

% /n/
dur = 1; % in s
Fbeg = [257 1891]'; Fend = [257 1891]';
Fbwbeg = [32 100]'; Fbwend = [32 100]';
Zbeg = [1223]'; Zend = [1223]';
Zbwbeg = [52]'; Zbwend = [52]';

runSynth_OLA_wrapper

audio = x{1,1};
wavwritedm(audio, fs, 16, 'n_synth.wav')

trueState1 = trueState{1,1};

% /n a/
dur = .5; % in s
Fbeg = [257 1891]'; Fend = [850 1500]';
Fbwbeg = [32 100]'; Fbwend = [80 120]';
Zbeg = []'; Zend = []';
Zbwbeg = []'; Zbwend = []';

runSynth_OLA_wrapper

audio = x{1,1};
wavwritedm(audio, fs, 16, 'na_synth.wav')

trueState2 = trueState{1,1};

% /a/
dur = 1; % in s
Fbeg = [850 1500]'; Fend = [850 1500]';
Fbwbeg = [80 120]'; Fbwend = [80 120]';
Zbeg = []'; Zend = []';
Zbwbeg = []'; Zbwend = []';

runSynth_OLA_wrapper

audio = x{1,1};
wavwritedm(audio, fs, 16, 'a_synth.wav')

trueState3 = trueState{1,1};

% /a n/
dur = .5; % in s
Fbeg = [850 1500]'; Fend = [500 2100]';
Fbwbeg = [80 120]'; Fbwend = [32 50]';
Zbeg = []'; Zend = []';
Zbwbeg = []'; Zbwend = []';

runSynth_OLA_wrapper

audio = x{1,1};
wavwritedm(audio, fs, 16, 'an_synth.wav')

trueState4 = trueState{1,1};

% second /n/
dur = 1; % in s
Fbeg = [500 2100]'; Fend = [500 2100]';
Fbwbeg = [32 50]'; Fbwend = [32 50]';
Zbeg = [1800]'; Zend = [1800]';
Zbwbeg = [52]'; Zbwend = [52]';

runSynth_OLA_wrapper

audio = x{1,1};
wavwritedm(audio, fs, 16, 'n_synth2.wav')

trueState5 = trueState{1,1};

%%
n = wavread('n_synth.wav');
na = wavread('na_synth.wav');
a = wavread('a_synth.wav');
an = wavread('an_synth.wav');
n2 = wavread('n_synth2.wav');

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
nan((56)*500+1:(56)*500+(19+1)*500) = nan((56)*500+1:(56)*500+(19+1)*500) + n2;
plot(nan, 'y-')

wavwritedm(nan, fs, 16, 'nan_synth9.wav')

%% true state
trueState2b = [trueState2; zeros(2, size(trueState2,2))];
trueState2b(5:6,1:end) = NaN;
trueState3b = [trueState3; zeros(2, size(trueState3,2))];
trueState3b(5:6,1:end) = NaN;
trueState4b = [trueState4; zeros(2, size(trueState4,2))];
trueState4b(5:6,1:end) = NaN;
trueState = [trueState1 trueState2b trueState3b trueState4b trueState5];
save('trueState9', 'trueState')