function f2_evaluateDNNs(numSpeakers, mixNumber, preset, featureType, bHeadMovement, azRes)
% f2_evaluateDNNs  Sound localisation using Neural Network.
%
%USAGE  
%  f2_evaluateDNNs(preset, nHiddenLayers, nHiddenNodes, bHeadMovement)
%
%INPUT ARGUMENTS
%    numSpeakers : number of simultaneous talkers, 1...3
%      mixNumber : Mixed utterance number (for parallel runs), 1...100
%         preset : 'CLEAN' or 'MCT-DIFFUSE'
%    featureType : 'itd-ild' or 'cc-ild'
%  bHeadMovement : flag for using head movement during localiation (false)
%          azRes : Azimuth resolution in degrees, default 5
%
% Ning Ma, 29 Jan 2015
%

if nargin < 2
    error('Please specify num speakers and num mixtures');
end


if nargin < 6
    azRes = 5;
end
if nargin < 5
    bHeadMovement = false;
end
if nargin < 4
    featureType = 'ild-cc'; % 'ild-cc', 'itd-ild', or 'cc'
end
if nargin < 3
    preset = 'MCT-DIFFUSE'; % 'CLEAN';
end

% Parameters
nHiddenLayers = 2;

%% Install software 
% 
[dataRoot, twoearsRoot] = get_data_root;

modelRoot = fullfile(pwd,'..','models');

% Get to correct directory and add working directories to path
gitRoot = fileparts(fileparts(mfilename('fullpath')));

% Add TwoEars WP1 functionality
addpath(genpath([twoearsRoot, filesep, 'binaural-simulator', filesep, 'src']));

% Add TwoEars AFE functionality
addpath(genpath([twoearsRoot, filesep, 'auditory-front-end', filesep, 'src']));

% Add TwoEars tools
addpath(genpath([twoearsRoot, filesep, 'main', filesep, 'src']));

% Add local tools
addpath Tools

% Add common scripts
addpath([gitRoot, filesep, 'tools', filesep, 'common']);

% Add common scripts
addpath(genpath([gitRoot, filesep, 'tools', filesep, 'DeepLearnToolbox']));


% Reset internal states of random number generator. This allows to use
% different settings, while still obtaining the "same" random matrix with
% sound source positions.
try
    % Since MATLAB 2011a 
    rng(0);
catch
    % For older versions
    rand('seed',0); %#ok
end


%% Experimental parameters
% 
AFE_param = initialise_AFE_parameters;
nChannels = AFE_param.fb_nChannels;
strEvaluationSetting = fullfile(dataRoot, 'Evaluation_localisation_setting');
load(strEvaluationSetting);

strRootModels = fullfile(modelRoot, 'LearnedDNNs');
strRootModels = sprintf('%s_%s_%s_%ddeg_%dchannels', strRootModels, preset, featureType, azRes, nChannels);
strRootResults = fullfile(dataRoot, 'EvaluationDNN');
strRootResults = sprintf('%s_%s_%s_%ddeg_%dchannels', strRootResults, preset, featureType, azRes, nChannels);
if ~exist(strRootResults, 'dir')
    mkdir(strRootResults);
end
strResults = fullfile(strRootResults, sprintf('evaluation_%s_%dlayers', preset, nHiddenLayers));
if bHeadMovement
    strResults = strcat(strResults, '_headMovement');
end

NNs = cell(nChannels, 1);
normFactors = cell(nChannels, 1);

for c = 1:nChannels
    strModels = fullfile(strRootModels, sprintf('DNN_%s_%ddeg_%dchannels_channel%d_%dlayers', ...
        preset, azRes, nChannels, c, nHiddenLayers));
    % Load localisation module
    load(strModels);
    NNs{c} = C.NNs;
    normFactors{c} = C.normFactors;
end
C.NNs = NNs;
C.normFactors = normFactors;
C.AFE_param = AFE_param;
clear NNs normFactors;

% Sampling frequency in Herz of the noisy speech mixtures
fsHz_Ref = 16E3;


%% Measure computation time
tstart = tic;
% Counter
iter = 0;
niters = nRooms*nAzim;

fprintf('\n');

azEstEval = cell(nAzim,nRooms);

% Loop over the number of speaker conditions
jj = numSpeakers;

% Loop over number of acoustic mixtures
ii = mixNumber;

% Read audio signals
audio = readAudio(files{ii,jj},fsHz_Ref);

% Loop over number of azimuth directions
for aa = 1 : nAzim

    % Loop over number of different rooms
    for rr = 1 : nRooms

        % Report progress
        fprintf('Localization progress %.1f%%: mixture %d/%d, azimuth %d/%d, speaker %d/%d, room %d/%d (%s)\n',...
                100*iter/niters, ii, nMixtures, aa, nAzim, jj, nSCond,  rr, nRooms, rooms{rr});
            
        fprintf('Ref: %d\n', azRefEval{ii,aa,jj});

        sim = RobotSimulator(audio, fsHz_Ref, azRefEval{ii,aa,jj}, rooms{rr});
        azEstWP1 = localise_model('DNN', sim, C, nSpeakers(jj), bHeadMovement, azRefEval{ii,aa,jj});
        azEstEval{aa,rr} = azEstWP1;
        for n=1:length(azEstWP1)
            fprintf('Loc: %.1f\n', azEstWP1(n));
        end

        % Increase counter
        iter = iter + 1;

        % Workout remaining time
        avgtime = toc(tstart) / iter;
        remtime = avgtime * (niters - iter + 1);
        days = remtime/3600/24;
        if days >= 1
            fprintf('\n---- %s: estimated remaining time: %d days %s\n', datestr(now), floor(days), datestr(rem(days,1), 'HH:MM:SS'));
        else
            fprintf('\n---- %s: estimated remaining time: %s\n', datestr(now), datestr(days, 'HH:MM:SS'));
        end
    end
end

strResults = sprintf('%s_nspk%d_utt%d', strResults, numSpeakers, mixNumber);
save(strResults, 'azEstEval');

