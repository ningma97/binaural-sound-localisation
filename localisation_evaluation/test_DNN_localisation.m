function [probs, azimuths] = test_DNN_localisation(sig, fsHz)
%
%USAGE  
%  [prob_AFN_F] = test_DNN_localisation(sig, fsHz)
%
%INPUT ARGUMENTS
%    sig : binaural signals
%   fsHz : sampling rate
%
% Ning Ma, 29 Jan 2015
%

%% Parameters
%
nHiddenLayers = 2;
azRes = 5;
featureType = 'ild-cc'; % 'ild-cc', 'itd-ild', or 'cc'
preset = 'MCT-DIFFUSE'; % 'CLEAN';

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

%% Setup DNNs
% 
AFE_param = initialise_AFE_parameters;
nChannels = AFE_param.fb_nChannels;
strEvaluationSetting = fullfile(dataRoot, 'Evaluation_localisation_setting');
load(strEvaluationSetting);

strRootModels = fullfile(modelRoot, 'LearnedDNNs');
strRootModels = sprintf('%s_%s_%s_%ddeg_%dchannels', strRootModels, preset, featureType, azRes, nChannels);

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

if fsHz ~= fsHz_Ref
    sig = resample(sig, fsHz_Ref, fsHz);
end

probs = compute_posteriors_DNN(sig, fsHz_Ref, C);

azimuths = convertAzimuthsWP1ToSurrey(C.azimuths(:));
[azimuths, si] = sort(azimuths);
probs = probs(si);

if nargout < 1
    plot(azimuths, probs);
    axis([-180 180 0 1]);
    xlabel('Azimuth', 'FontSize', 12);
    ylabel('Probability', 'FontSize', 12);
end


