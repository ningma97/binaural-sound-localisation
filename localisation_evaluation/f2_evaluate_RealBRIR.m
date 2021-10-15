function f2_evaluate_RealBRIR(model, bHeadMovement, preset, featureType)
% Sound localisation using realistic BRIRs (head rotation)
%
%USAGE  
%  f2_evaluate_RealBRIR(model, bHeadMovement, preset, featureType)
%
%INPUT ARGUMENTS
%         model : 'DNN' or 'GMM'
% bHeadMovement : flag for using head movement during localiation (false)
%         preset : 'CLEAN' or 'MCT-DIFFUSE'
%    featureType : 'itd-ild' or 'cc-ild'
%
% Ning Ma, 29 Jan 2015
%

% rooms = {'auditorium3', 'spirit'};
room = 'auditorium3';

azRes = 5;
if nargin < 1
    error('Please specify model DNN or GMM');
end
if nargin < 2
    bHeadMovement = false;
end
if nargin < 3
    preset = 'MCT-DIFFUSE'; % 'CLEAN';
end
if nargin < 4
    featureType = 'ild-cc'; % 'ild-cc', 'itd-ild', or 'cc'
end


%% Install software 
% 
[~, twoearsRoot] = get_data_root;
dataRoot = '../results';

modelRoot = fullfile(pwd,'..','models');

% Get to correct directory and add working directories to path
gitRoot = fileparts(fileparts(mfilename('fullpath')));

% Add TwoEars WP1 functionality
addpath(genpath([twoearsRoot, filesep, 'binaural-simulator', filesep, 'src']));

% Add TwoEars AFE functionality
addpath(genpath([twoearsRoot, filesep, 'auditory-front-end', filesep, 'src']));

% Add TwoEars tools
addpath(genpath([twoearsRoot, filesep, 'main', filesep, 'src']));

% Add SOFA tools
addpath(genpath([twoearsRoot, filesep, 'SOFA', filesep, 'API_MO']));

% Add local tools
addpath Tools

% Add common scripts
addpath([gitRoot, filesep, 'tools', filesep, 'common']);
addpath(genpath([gitRoot, filesep, 'tools', filesep, 'DeepLearnToolbox']));
addpath([gitRoot, filesep, 'tools', filesep, 'GMM_Netlab']);

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
% String specifying different BRIRs
%rooms = {'auditorium3', 'spirit'};

% Head orientation index:
% 1:90   -> right
% 91     -> front
% 92:181 -> left
frontHeadIdx = 91;

% Define number of competing speech sources
nSpeakers = [1 2 3];

% Number of acoustic mixtures for each acoustic condition
nMixtures = 100;

% minimum source distance
minDistDeg = 10;

% Initial head position
initHeadOrientation = 0;

% Rooms 
switch lower(room)
    case 'spirit'
        brirs{1} = SOFAload(db.getFile('impulse_responses/qu_kemar_rooms/spirit/QU_KEMAR_spirit_src1_30deg.sofa'));
        brirs{2} = SOFAload(db.getFile('impulse_responses/qu_kemar_rooms/spirit/QU_KEMAR_spirit_src2_0deg.sofa'));
        brirs{3} = SOFAload(db.getFile('impulse_responses/qu_kemar_rooms/spirit/QU_KEMAR_spirit_src3_-30deg.sofa'));
        nSourcePositions = length(brirs);
        azRangeEval = zeros(nSourcePositions,1);
        for n = 1 : nSourcePositions
            asv = SOFAcalculateAPV(brirs{n});
            azRangeEval(n) = convertAzimuthsSurreyToWP1(-1*asv(frontHeadIdx));
        end

    case 'auditorium3'
        
        brirs{1} = SOFAload(db.getFile('impulse_responses/qu_kemar_rooms/auditorium3/QU_KEMAR_Auditorium3_src5_xs-0.75_ys+1.30.sofa'));
        brirs{2} = SOFAload(db.getFile('impulse_responses/qu_kemar_rooms/auditorium3/QU_KEMAR_Auditorium3_src1_xs+0.00_ys+3.97.sofa'));
        brirs{3} = SOFAload(db.getFile('impulse_responses/qu_kemar_rooms/auditorium3/QU_KEMAR_Auditorium3_src6_xs+0.75_ys+1.30.sofa'));
        brirs{4} = SOFAload(db.getFile('impulse_responses/qu_kemar_rooms/auditorium3/QU_KEMAR_Auditorium3_src2_xs+4.30_ys+3.42.sofa'));
        brirs{5} = SOFAload(db.getFile('impulse_responses/qu_kemar_rooms/auditorium3/QU_KEMAR_Auditorium3_src3_xs+2.20_ys-1.94.sofa'));
        
        %brirs{4} = SOFAload(db.getFile('impulse_responses/qu_kemar_rooms/auditorium3/QU_KEMAR_Auditorium3_src4_xs+0.00_ys+1.50.sofa'));
        
        nSourcePositions = length(brirs);
        azRangeEval = zeros(nSourcePositions,1);
        %srcDistance = zeros(nSourcePositions,1);
        for n = 1 : nSourcePositions
            asv = SOFAcalculateAPV(brirs{n});
            azRangeEval(n) = convertAzimuthsSurreyToWP1(-1*asv(frontHeadIdx));
            %srcDistance(n) = sqrt(brirs{n}.SourcePosition(1)^2+brirs{n}.SourcePosition(2)^2);
        end
        % Load scaling factors for BRIRs to correct source level difference
        load(strcat('RealBRIR_scaling_factors_',lower(room))); 

    otherwise
        error('Room is not supported!')
end

%% Model parameters
% 
AFE_param = initialise_AFE_parameters;
nChannels = AFE_param.fb_nChannels;

% Load models
switch model
    case 'DNN'
        nHiddenLayers = 2;
        strRootModels = fullfile(modelRoot, 'LearnedDNNs');
        strRootModels = sprintf('%s_%s_%s_%ddeg_%dchannels', strRootModels, preset, featureType, azRes, nChannels);
        strRootResults = fullfile(dataRoot, 'EvaluationDNN');
        strRootResults = sprintf('%s_RealBRIR_%s_%s_%ddeg_%dchannels', strRootResults, preset, featureType, azRes, nChannels);
        if ~exist(strRootResults, 'dir')
            mkdir(strRootResults);
        end
        strResults = fullfile(strRootResults, sprintf('evaluation_%s_%s_%dlayers', room, preset, nHiddenLayers));
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

    case 'GMM'
        nMix = 16;

        strRootModels = fullfile(modelRoot, 'LearnedGMMs');
        strModels = fullfile(strRootModels, sprintf('GMM_%s_%s_%ddeg_%dchannels_%dmix_Norm', preset, featureType, azRes, nChannels, nMix));
        load(strModels);

        strRootResults = fullfile(dataRoot, 'EvaluationGMM');
        strRootResults = sprintf('%s_RealBRIR_%s_%s_%ddeg_%dchannels', strRootResults, preset, featureType, azRes, nChannels);
        if ~exist(strRootResults, 'dir')
            mkdir(strRootResults);
        end
        strResults = fullfile(strRootResults, sprintf('evaluation_%s_%s_%dmix', room, preset, nMix));
        if bHeadMovement
            strResults = strcat(strResults, '_headMovement');
        end
        
    otherwise
        error('Model %s is not supported!', model);
end
% Sampling frequency in Herz of the noisy speech mixtures
fsHz_Ref = 16E3;


%% Initialise speech source parameters
% 
% 
evalFile = fullfile(dataRoot, 'Evaluation_localisation_RealBRIR_setting.mat');
if exist(evalFile, 'file')
    load(evalFile);
else

    % Root directory of GRID database
    rootGRID = fullfile(twoearsRoot, 'twoears-data/sound_databases/grid_subset');

    % Test set
    flist = 'flist_test.txt';
    if ~exist(flist, 'file')
        error('Please generate test file list %s', flist);
    end
    fid = fopen(flist);
    allFiles = textscan(fid, '%s');
    fclose(fid);
    allFiles = allFiles{1};
    nSentences = numel(allFiles);

    % Check if enough material is available
    if nSentences < nMixtures*max(nSpeakers)
       error('Not enough audio material available.')
    end


    % Number of different conditions
    nAzim  = numel(azRangeEval);
    nSCond = numel(nSpeakers);

    files = cell(nMixtures,nSCond);
    azRefEvalIdx = cell(nMixtures,nAzim,nSCond);
    azEstEval = cell(nMixtures,nAzim,nSCond);

    % Loop over the number of speaker conditions
    for jj = 1 : nSCond

        randFileIdx = randperm(nSentences, nMixtures*nSpeakers(jj));
        azIdx = zeros(nSpeakers(jj), 1);

        % Loop over number of acoustic mixtures
        for ii = 1 : nMixtures

            % Randomly select "nSpeakers" sentences
            fileIdx = randFileIdx((ii-1)*nSpeakers(jj)+1:ii*nSpeakers(jj));
            files{ii,jj} = cell(numel(fileIdx),1);
            for n = 1:numel(fileIdx)
                fn = sprintf('%s/%s.wav', rootGRID, strrep(allFiles{fileIdx(n)}, '_', '/'));
                files{ii,jj}{n} = fn;
            end

%             % Loop over number of azimuth directions
%             for aa = 1 : nAzim
%                 azIdx(1) = aa;
%                 if nSpeakers(jj) > 1
%                     azIdxAll = randperm(nSourcePositions);
%                     azIdxAll = azIdxAll(azIdxAll~=aa);
%                     azIdx(2:end) = azIdxAll(1:nSpeakers(jj)-1);
%                 end
%                 azRefEvalIdx{ii,aa,jj} = azIdx;
%             end
            
            % Loop over number of azimuth directions
            for aa = 1 : nAzim
                azIdx(1) = aa;
                azEvalIdx = 1:nAzim;
                azEvalIdx = azEvalIdx(azEvalIdx~=aa);
                % Enforce a "minDistance" spacing between all sources
                for ss = 2:nSpeakers(jj)
                    azDist = calc_azimuth_distance(azRangeEval(azIdx(ss-1)), azRangeEval(azEvalIdx));
                    azEvalIdx = azEvalIdx(azDist > minDistDeg);
                    azIdx(ss) = azEvalIdx(randperm(numel(azEvalIdx), 1));
                end
                azRefEvalIdx{ii,aa,jj} = azIdx;
            end
            
        end
    end
    
    save(evalFile, 'files', 'azRefEvalIdx', 'azEstEval', 'nAzim', 'nSCond');

end

%% Measure computation time
tstart = tic;
% Counter
iter = 0;
niters = nMixtures*nSCond*nAzim;

fprintf('\n');

% Loop over the number of speaker conditions
for jj = 1 : nSCond
    
    % Loop over number of acoustic mixtures
    for ii = 1 : nMixtures
        
        % Read audio signals
        audio = readAudio(files{ii,jj},fsHz_Ref);

        % Loop over number of azimuth directions
        for aa = 1 : nAzim

            azRef = azRangeEval(azRefEvalIdx{ii,aa,jj});
            % Report progress
            fprintf('Localization progress %.1f%%: mixture %d/%d, azimuth %d/%d, speaker %d/%d, room %s\n',...
                    100*iter/niters, ii, nMixtures, aa, nAzim, jj, nSCond, room);
            
            fprintf('Ref: %.0f\n', azRef);

            sim = RobotRealBRIR(audio, fsHz_Ref, azRefEvalIdx{ii,aa,jj}, brirs, brirScales);
            sim.rotateHead(initHeadOrientation);
            
            az = localise_model(model, sim, C, nSpeakers(jj), bHeadMovement);
            azEstWP1 = mod(az+initHeadOrientation, 360);
            azEstEval{ii,aa,jj} = azEstWP1;
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
end



save(strResults, 'azEstEval', 'azRefEvalIdx', 'azRangeEval', 'room', 'bHeadMovement', 'preset', 'featureType', 'model');

