% Under development for testing with synthetic vowels
function [s,fs] = genVowel()


% Set random seeds
rand('seed',134);
randn('seed',2342);

% Set options
doPlots = 1; % Plot results 1/0 = yes/no
doSound = 0; % Listen to results 1/0 = yes/no
doWrite = 0; % Write out results 1/0 = yes/no

fileName = 'Data/SynthVowel/';

% Synthesize vowel
gender = 'f';  % 'f' or 'm'
vowel  = 'e';  % Can be: a,o,e,i,u
fs     = 16000;% Sampling Rate
len    = 1.0;    % In seconds
pitchContour = [200 200];

fileName = strcat(fileName, vowel);

s = synthVowel(vowel, gender, fs, len, pitchContour);

if(doPlots)
   figure;
   plot(s);
   title('Synthesized Vowel');
end

function [s] = synthVowel(vowel, gender, fs, len, pitchContour)

playSound = 0;

[F BW] = getVTFilterParams(vowel, gender);
s = genVowelPB(pitchContour, fs, len, F, BW);

if(playSound)
   soundsc(s,fs);
end
% Generate a vowel based on the Petersen-Barney
% Formant specification
%
%
% INPUT
% Vowel :       Which vowel: /a/, /i/, /e/,/o/,/u/, /ae/
% Sex   :       Male or Female: 'm', 'f'
% Pitch :       Vector of pitch values [start finish]
% fs    :       Sampling rate
% duration:       Length in seconds
%
% Example call:
% Lookup formant locations and bandwidths for vowel/gender pairs
% [F BW] = getVTFilterParams(vowel, gender);
% vowelOut = getVowelPB([100 125],16000,1, F, BW);

function [vowelOut] = genVowelPB(pitch, fs, duration, F, BW)

% Argument checking

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Vocal Tract Filter Construction%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute filter coefficients (rational form)

a=[];b=[];c=[];
for i = 1:3
   b(i) = 2 * cos(2*pi*F(i)/fs) * exp(-pi*BW(i)/fs);
   c(i) = -exp(-2*pi*BW(i)/fs);
   a(i) = 1 - b(i) - c(i);
end

% Convolve to get coefficients
A = conv(conv([ 1 -b(1) -c(1)],[ 1 -b(2) -c(2)]),[ 1 -b(3) -c(3)]);
B = a(1)*a(2)*a(3);

% Plot frequency response of filter
figure;
freqz(B,A,512,fs) % freq response of filter
figure;
freqz(B,A,0:.01:pi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Glottal Train %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
glottalTrain = getGlottalTrain(pitch, fs, duration);

sLength = fs*duration;
sLength = sLength-mod(sLength, 512);
% Hacked to just do the stupid thing and drive with noise
glottalTrain = randn(1, sLength);
glottalTrain = glottalTrain/max(glottalTrain);

% Filter the glottal train through the vocal tract filter
vowelOut = filter(B, A, glottalTrain);    % Filter
vowelOut = vowelOut(1,1:floor(sLength));  % Truncate to length
%vowelOut = vowelOut/(max(abs(vowelOut)));       % Scale
db = 5;
% Listen
% soundsc(vowelOut,fs);

% Glottal train genesis based on Daryush's function
function[glottalTrain] = getGlottalTrain(f0, fs, duration)

sLength = fs*duration;
sLength = sLength-mod(sLength, 512);
% Generate pitch contour
pitchContour = [1 f0(1); floor(sLength) f0(2)]; % only 2 endpoints specified
pitchContour = interp1( pitchContour(:,1), pitchContour(:,2), 1:floor(sLength));
%pitchContour = synthPitchContour;

% Put in a rapidly varying sinusoid
pitchContour = sin((1:length(pitchContour))/100)+100;


% Glottal flow
t = 0:1/(1*fs):1;
cumUprime = []; cumU = []; done_flag = 0;
timeOffset = 0.001;
OQ = 60;      % Open quotient in percent
dcFlow = 0.1; % fraction of maximum noise amplitude for DC air

f0 = pitchContour(length(cumU)+1);
T0 = 1/f0;

while done_flag == 0,
    f0 = pitchContour(length(cumU)+1);
    T0 = 1/f0;

    % U is the velocity flow
    U = ( t.^2 - ( ( t.^3 ) ./ ( OQ/100*T0 ) ) );
    U = U(1:round(T0*fs));
    U(floor(( (size(U,2))*OQ/100) ):size(U,2)) = 0;
    U = U/(max(abs(U)));
    U = U + dcFlow*max(abs(U));

    cumU = [cumU U];
    if length(cumU) >= floor(sLength)
        cumU = cumU(1:floor(sLength));
        done_flag = 1;
    end
end
glottalTrain = cumU;
glottalTrain = glottalTrain(1:floor(sLength));
