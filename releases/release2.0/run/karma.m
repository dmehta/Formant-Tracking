% karma.m
%
% Author: Daniel Rudoy, Daryush Mehta
% Created:  12/13/2007
% Last modified: 08/01/2011
% 
% The purpose of this function is to obtain frequency and bandwidth values
% for formant and antiformant trajectories of a given speech waveform. A
% Kalman-based autoregressive moving average approach is employed [1-3].
% 
% References:
% [1]   D. D. Mehta, D. Rudoy, and P. J. Wolfe, "Kalman-based autoregressive
%       moving average modeling and inference for formant and antiformant
%       tracking," In review, 2011.
% [2]   D. Rudoy, D. N. Spendley, and P. J. Wolfe, "Conditionally linear
%       Gaussian models for estimating vocal tract resonances," Proceedings of
%       Interspeech, Antwerp, Belgium, 2007.
% [3]   D. Rudoy, "Nonstationary time series modeling with application to speech
%       signal processing," Doctor of Philosophy thesis, School of Engineering
%       and Applied Sciences, Harvard University, Cambridge, MA, 2010.
%       Chapter 3.
% 
% function [x_est, x_errVar, wav, params] = runWaveZ(input, varargin)
%       For inputting a wave file. See below for optional parameters. Set
%       to [] to keep as default
% function [x_est, x_errVar, wav, params] = runWaveZ(input, fs_in, varargin)
%       For inputting a vector. See below for optional parameters. Set
%       to [] to keep as default
% 
% INPUT:
%    input:                     Relative path to Microsoft WAV file (include .wav extension)
%                                   or Matlab vector containing speech waveform
%    fs_in:                     Sampling rate of input, in Hz (omit if .wav input)
%    numFormants (optional):    Number of formants to track; default: 3
%    numAntiF (optional):       Number of anti-formants to track; default: 0
%    aParams (optional):        Analysis parameters, structure with one or more of the following fields:
%                                   peCoeff     % Pre-emphasis coefficient; default: 0.7
%                                   wType       % window type; default: 'Hamming'
%                                   wLengthMS   % Length of window, in ms; default: 20
%                                   wOverlap    % Window overlap factor; default: 0.5
%                                   lpcOrder    % Number of AR coefficients; default: 12
%                                   zOrder      % Number of MA coefficients; default: 0
%                                   fs          % Downsampling frequency, in Hz; default: 7000
%    cepOrder (optional):       Number of cepstral coefficients to observe; default: 15
%    cepType (optional):        Type of cepstrum to compute; 1 = ARMA, 2 = real
%                                   default: 1
%    algFlag (optional):        Flag for algorithm to run, 1 for EKF, 2 for EKS;
%                                   default: 2
%    x0 (optional):             Column vector of initialized center frequency values
%                                   for formants then antiformants, in Hz
%                                   Length is (numFormants + numAntiformants);
%                                   default: [500 1500 2500 ... ] for formants
%                                   [1000 2000 3000 ... ] for antiformants
%    bwStates (optional):       Column vector of initialized bandwidth values
%                                   for formants then antiformants, in Hz
%                                   Length is (numFormants + numAntiformants);
%                                   default: [80 120 160 ... ] for formants
%                                   [80 120 160 ... ] for antiformants
%    formantInds (optional):        Use speech activity detection? 1 for yes (tracks will be
%                                   coasted during non-speech frames), [] for no; 
%                                   can also input masking matrix of
%                                   formant indices: 0 if masked, 1 if
%                                   unmasked; default: 1
%    trackBW (optional):        Flag to calculate bandwidths (1) or fix to
%                                   initial values in bwStates (2)
%                                   default: 1
% 
% OUTPUT:
%    x_est:                     Matrix of estimated tracks
%    x_errVar:                  Covariance matrix of estimated tracks
%    wav:                       Downsampled waveform returned as column vector
%    params:                    Parameters applied, for plotting routines, etc.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USAGE
% see karma_demo.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x_est, x_errVar, wav, params] = karma(input, varargin)

%% load default parameters and then overwrite user-defined values
[numFormants, numAntiF, aParams, cepOrder, cepType, algFlag, x0, bwStates, formantInds, trackBW] = getDefault();

if strcmpi(input(end-3:end), '.wav')
    [wav, fs_in] = wavread(input);
    
    if ~isempty(varargin)
        if ~isempty(varargin{1}), numFormants = varargin{1};
            % Set initial state of formant trackers
            % Using values suggested by Li Deng and otherwise in literature
            x0 = 500 + 1000*(0:(numFormants - 1))';
            x02 = 1000 + 1000*(0:(numAntiF - 1))';
            x0 = [x0; x02];
            
            bwStates = 80 + 40*(0:(numFormants - 1))';
            bwStates2 = 80 + 40*(0:(numAntiF - 1))';
            bwStates = [bwStates; bwStates2];
        end
    end
    if length(varargin) > 1
        if ~isempty(varargin{2}), numAntiF = varargin{2};
            % Set initial state of formant trackers
            % Using values suggested by Li Deng and otherwise in literature
            x0 = 500 + 1000*(0:(numFormants - 1))';
            x02 = 1000 + 1000*(0:(numAntiF - 1))';
            x0 = [x0; x02];
            
            bwStates = 80 + 40*(0:(numFormants - 1))';
            bwStates2 = 80 + 40*(0:(numAntiF - 1))';
            bwStates = [bwStates; bwStates2];
        end
    end
    if length(varargin) > 2
        if ~isempty(varargin{3}), 
            aParamsNew = varargin{3};
            aParamsNewC = struct2cell(aParamsNew);
            aParamsNewFields = fieldnames(aParamsNew); 
            
            aParamsFields = fieldnames(aParams);
            aParamsC = struct2cell(aParams);
            for ii = 1:length(aParamsNewC)
                indToChange = strcmpi(aParamsNewFields(ii), aParamsFields);
                aParamsC(indToChange) = aParamsNewC(ii);
            end
            aParams = cell2struct(aParamsC, aParamsFields, 1);
        end
    end
    if length(varargin) > 3
        if ~isempty(varargin{4}), cepOrder = varargin{4}; end
    end
    if length(varargin) > 4
        if ~isempty(varargin{5}), cepType = varargin{5}; end
    end
    if length(varargin) > 5
        if ~isempty(varargin{6}), algFlag = varargin{6}; end
    end
    if length(varargin) > 6
        if ~isempty(varargin{7}), x0 = varargin{7}; end
    end
    if length(varargin) > 7
        if ~isempty(varargin{8}), bwStates = varargin{8}; end
    end
    if length(varargin) > 8
        if ~isempty(varargin{9}), formantInds = varargin{9}; end
    end
    if length(varargin) > 9
        if ~isempty(varargin{10}), trackBW = varargin{10}; end    
    end
else % input is a vector
    wav = input;
    fs_in = varargin{1};

    if length(varargin) > 1
        if ~isempty(varargin{2}), numFormants = varargin{2}; 
            % Set initial state of formant trackers
            % Using values suggested by Li Deng and otherwise in literature
            x0 = 500 + 1000*(0:(numFormants - 1))';
            x02 = 1000 + 1000*(0:(numAntiF - 1))';
            x0 = [x0; x02];
            
            bwStates = 80 + 40*(0:(numFormants - 1))';
            bwStates2 = 80 + 40*(0:(numAntiF - 1))';
            bwStates = [bwStates; bwStates2];
        end            
    end
    if length(varargin) > 2
        if ~isempty(varargin{3}), numAntiF = varargin{3};
            % Set initial state of formant trackers
            % Using values suggested by Li Deng and otherwise in literature
            x0 = 500 + 1000*(0:(numFormants - 1))';
            x02 = 1000 + 1000*(0:(numAntiF - 1))';
            x0 = [x0; x02];
            
            bwStates = 80 + 40*(0:(numFormants - 1))';
            bwStates2 = 80 + 40*(0:(numAntiF - 1))';
            bwStates = [bwStates; bwStates2];
        end
    end
    if length(varargin) > 3
        if ~isempty(varargin{4}), 
            aParamsNew = varargin{4};
            aParamsNewC = struct2cell(aParamsNew);
            aParamsNewFields = fieldnames(aParamsNew); 
            
            aParamsFields = fieldnames(aParams);
            aParamsC = struct2cell(aParams);
            for ii = 1:length(aParamsNewC)
                indToChange = strcmpi(aParamsNewFields(ii), aParamsFields);
                aParamsC(indToChange) = aParamsNewC(ii);
            end
            aParams = cell2struct(aParamsC, aParamsFields, 1);
        end        
    end
    if length(varargin) > 4
        if ~isempty(varargin{5}), cepOrder = varargin{5}; end
    end
    if length(varargin) > 5
        if ~isempty(varargin{6}), cepType = varargin{6}; end
    end        
    if length(varargin) > 6
        if ~isempty(varargin{7}), algFlag = varargin{7}; end
    end        
    if length(varargin) > 7
        if ~isempty(varargin{8}), x0 = varargin{8}; end
    end        
    if length(varargin) > 8
        if ~isempty(varargin{9}), bwStates = varargin{9}; end
    end        
    if length(varargin) > 9
        if ~isempty(varargin{10}), formantInds = varargin{10}; end
    end        
    if length(varargin) > 10
        if ~isempty(varargin{11}), trackBW = varargin{11}; end    
    end        
end

% Make sure that wav is a column vector
[row, col] = size(wav);
if col > row % if row vector, make col vector
    wav = wav';
end

%% Pre-processing steps

% Analysis parameters
wType       = aParams.wType;     % Window type
wLengthMS   = aParams.wLengthMS; % Length of window (in milliseconds)
wOverlap    = aParams.wOverlap;  % Factor of overlap of window
lpcOrder    = aParams.lpcOrder;  % Number of AR Coefficients
zOrder      = aParams.zOrder;    % Number of MA Coefficients
peCoeff     = aParams.peCoeff;   % Pre-emphasis factor
fs          = aParams.fs;        % Downsampling rate

% Downsample
if(fs ~= fs_in)
    display(['Original fs = ' int2str(fs_in) ' Hz; resampling to ' num2str(fs) ' Hz'])
    wav = resample(wav,fs,fs_in,2048);
end

% Compute window length in samples, now that sampling rate is set
wLength = floor(wLengthMS/1000*fs);
wLength = wLength + (mod(wLength,2)~=0); % Force even
win = feval(wType,wLength);

% Pack the analysis parameters for later use
aParams.wLength    = wLength;   % Length of window (in samples)
aParams.win        = win;       % The actual window

%% Intraframe observation sequence
% Generate cepstral data
switch cepType
    case 1 % ARMA
        y = genLPCCz(wav, win, wOverlap, peCoeff, lpcOrder, zOrder, cepOrder);
    case 2 % real
        y = genCeps(wav, win, wOverlap, peCoeff, lpcOrder, zOrder, cepOrder);
    otherwise
        error('Unknown cepType input.')
end

%% Interframe tracking
%%% Set parameters for tracking algorithms %%%
% x_{k+1} = Fx_{k} + w_k, w_k ~ N(0, Q)
% y_{k}   = Hx_{k} + v_k, v_k ~ N(0, R)
% We need to set the parameters: F, Q and R
% H is obtained in the EKF via linearization about the state
nP = numFormants;
numObs = size(y,2);

if trackBW
    stateLen = 2*numFormants + 2*numAntiF;
    numAntiF = length(x0)-numFormants;
    x0 = [x0(1:numFormants); bwStates(1:numFormants); ...
        x0(numFormants+1:numFormants+numAntiF); bwStates(numFormants+1:numFormants+numAntiF)];
    bwStates = []; % If we are tracking bandwidths do not provide them

    % Process noise covariance matrix Q
    Q = zeros(size(x0));
    Q(1:numFormants) = 10e4;
    Q(numFormants+1:numFormants*2) = 1e4;
    Q(2*numFormants+1:2*numFormants+numAntiF) = 10e4;
    Q(2*numFormants+numAntiF+1:end) = 1e4;
    Q = diag(Q);
else
    stateLen = numFormants + numAntiF;
    bwStates = repmat(bwStates, 1, numObs);
    
    % Process noise covariance matrix Q
    Q = diag(ones(stateLen,1)*10e4);
end

% formantInds is input
if isempty(formantInds)
    display('Not using Voice Activity Detection');
    formantInds = ones(stateLen, numObs);
else
    if length(formantInds) > 1
        % use formantInds input
    else
        coastJoint = 1;    % Coast all formants jointly, or not
        quantThresh = .15; % power quantile threshold for coasting
        plotVad = 0;
        multiBand = ~coastJoint;
        [frameInds] = simpleVAD(wav, aParams, quantThresh, multiBand, [], plotVad);
        if trackBW
            formantInds = repmat(frameInds,2*numFormants+2*numAntiF,1)';
        else
            formantInds = repmat(frameInds,numFormants+numAntiF,1)';
        end
        formantInds = formantInds';
    end
end

% Process Matrix F, Correlation is not being tested here
Fmatrix = eye(stateLen);

% Measurement noise covariance matrix R
% Set/choose estimated observation noise variance, which should decrease
% as the cepstral order increases
% Using sigExp = 2 for synthetic vowels, but =1 in formant experiments
lambda = 10^(0); sigExp = 1;
R = diag(1/lambda.*ones(cepOrder,1)./(((1:cepOrder).^sigExp)'));
% R = diag(1e-0*ones(cepOrder,1));

switch algFlag
    case 1 % Run Extended Kalman Filter
        smooth = 0;
        [x_est x_errVar] = formantTrackEKSZ(y, Fmatrix, Q, R, x0, formantInds, fs, bwStates, nP, smooth);    
    case 2 % Run Extended Kalman Smoother
        smooth = 1;
        [x_est x_errVar] = formantTrackEKSZ(y, Fmatrix, Q, R, x0, formantInds, fs, bwStates, nP, smooth);
    otherwise
        error('Unknown algorithm type (algFlag) set.')
end

%% pack parameters for output
clear params
params.numFormants = numFormants;
params.numAntiF = numAntiF;
params.aParams = aParams;
params.cepOrder = cepOrder;
params.cepType = cepType;
params.algFlag = algFlag;
params.x0 = x0;
params.bwStates = bwStates;
params.formantInds = formantInds;
params.trackBW = trackBW;