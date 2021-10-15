function f1_prepareEvaluation
% Prepares evaluation so that all experiments use the same setting in term
% of azimuth range, data, and BRIRs
%
%
%
% Ning Ma, 29 Jan 2015
%



%% Install software 
% 
[dataRoot, twoearsRoot] = get_data_root;

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
strEvaluationSetting = fullfile(dataRoot, 'Evaluation_localisation_setting');


% *************************************************************************
% List of different rooms:
% *************************************************************************
% 'ANECHOIC' 'ROOM_A'        'ROOM_B'         'ROOM_C'       'ROOM_D'
%  
%  RT60 = 0s  RT60 = 0.32s    RT60 = 0.47s    RT60 = 0.68s   RT60 = 0.89s
%  DRR  = inf DRR  = 6.09dB   DRR  = 5.31dB   DRR  = 8.82dB  DRR  = 6.12dB
% 
% 
% String specifying different BRIRs

%rooms = {'QU_ANECHOIC' 'SURREY_ANECHOIC' 'SURREY_ROOM_A' 'SURREY_ROOM_B' 'SURREY_ROOM_C' 'SURREY_ROOM_D'};
rooms = {'QU_ANECHOIC' 'SURREY_ROOM_A' 'SURREY_ROOM_B' 'SURREY_ROOM_C' 'SURREY_ROOM_D'};

% Azimuth range that should be used for testing using the WP1 convention
%azRangeEval = [60:-10:0 350:-10:300]'; % For evaluating Surrey BRIRs
azRangeEval = convertAzimuthsSurreyToWP1(-90:5:90); % For evaluating Surrey BRIRs

% Define number of competing speech sources
nSpeakers = [1 2 3];

% Minimum distance between competing sound sources in degree
minDistDeg = 10;

% Number of acoustic mixtures for each acoustic condition
nMixtures = 50;


%% Initialise speech source parameters
% 
% 
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
nRooms = numel(rooms);
nAzim  = numel(azRangeEval);
nSCond = numel(nSpeakers);

files = cell(nMixtures,nSCond);
azRefEval = cell(nMixtures,nAzim,nSCond);

% Loop over the number of speaker conditions
for jj = 1 : nSCond

    randFileIdx = randperm(nSentences, nMixtures*nSpeakers(jj));
    srcAz = zeros(nSpeakers(jj), 1);
    
    % Loop over number of acoustic mixtures
    for ii = 1 : nMixtures

        % Randomly select "nSpeakers" sentences
        fileIdx = randFileIdx((ii-1)*nSpeakers(jj)+1:ii*nSpeakers(jj));
        files{ii,jj} = cell(numel(fileIdx),1);
        for n = 1:numel(fileIdx)
            fn = sprintf('%s/%s.wav', rootGRID, strrep(allFiles{fileIdx(n)}, '_', '/'));
            files{ii,jj}{n} = fn;
        end

        % Loop over number of azimuth directions
        for aa = 1 : nAzim

            srcAz(1) = azRangeEval(aa);
            azEval = azRangeEval;
            % Enforce a "minDistance" spacing between all sources
            for ss = 2:nSpeakers(jj)
                azDist = calc_azimuth_distance(srcAz(ss-1), azEval);
                azEval = azEval(azDist > minDistDeg);
                % Remove -90 -85 85 90 from sources for Surrey database
                azEval = azEval(azEval <= 80 | azEval >= 280);
                azIdx = randperm(numel(azEval), 1);
                srcAz(ss) = azEval(azIdx);
            end
            azRefEval{ii,aa,jj} = srcAz;
        end
    end
end

save(strEvaluationSetting);
fprintf('Evaluation setting saved in %s\n', strEvaluationSetting);



