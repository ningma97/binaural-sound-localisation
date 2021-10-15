function [azEst] = localise_model(model, sim, C, nSources, bHeadRotation, azRef)
%localise_model   Localisation module.
%
%
%INPUT ARGUMENTS
%         model : Localisation model: 'DNN', 'GMM' etc
%           sim : Robot simulator
%             C : classifier structure used for localisation
%      nSources : number of sources to be localised (default 1)
% bHeadRotation : if true, rotate head half way through the signal (default false)
%         azRef : Reference azimuths. Only used to make sure we do not 
%                 rotate the head so that no BRIRs are available
%
%OUTPUT ARGUMENTS
%     azEst : estimated azimuths
% 
%
% Ning Ma, 29 Jan 2015
% n.ma@sheffield.ac.uk
%   
%

 
% Check for proper input arguments
if nargin < 4
    bHeadRotation = false;
end

% We use the middle segment (2 x blockDur) of signal for evaluation
blockDur = 0.5; % in seconds
% Skip the beginning
sigDur = sim.remainingDuration;
skipDur = (sigDur / 2) - blockDur;
sim.getSignal(skipDur);

% If head rotation is used, divide the signal into two blocks
if ~bHeadRotation
	blockDur = blockDur * 2;
end

[sig,fsHz] = sim.getSignal(blockDur);

post1 = feval(strcat('compute_posteriors_', model), sig, fsHz, C);

if bHeadRotation
%     % Randomly select a head rotation angle but make sure head 
%     % orientation stays inside [-30 30] for the Surrey database
%     rotationAngles = [60:-5:10 -10:-5:-60];
%     while true
%         rotateAngle = rotationAngles(randi(length(rotationAngles)));
%         newHO = mod(sim.headOrientation + rotateAngle, 360);
%         if newHO <= 30 || newHO >= 330
%             break;
%         end
%     end
    rotateAngle = getHeadRotationAngle(C.azimuths, post1);

    fprintf('Rotating head by %d degrees\n', rotateAngle);
    sim.rotateHead(rotateAngle);

    % Grab the rest of the signal and compute azimuth posteriors
    [sig,fsHz] = sim.getSignal(blockDur);
    post2 = feval(strcat('compute_posteriors_', model), sig, fsHz, C);

    % Combine two blocks
    [post1, post2] = removeFrontBackConfusion(C.azimuths, post1, post2, rotateAngle);
    idxDelta = round(rotateAngle / (C.azimuths(2) - C.azimuths(1)));
    post2 = circshift(post2, idxDelta);
    post = post1 + post2;
    prob_AFN_F = post ./ sum(post);
else
    prob_AFN_F = post1;
end

% Find peaks, also consider endpoints as peak candidates
pIdx = findpeaks([0; prob_AFN_F(:); 0]);
pIdx = pIdx - 1;

% Rank peaks
[~,idx] = sort(prob_AFN_F(pIdx),'descend');

% Number of azimuth estimates
nEst = min(numel(idx),nSources);

% Apply exponential interpolation to refine peak position
%delta = interpolateParabolic(prob_AFN_F,pIdx(idx(1:nEst)));

% Azimuth estimate: Take most significant peaks
% Determine azimuth step size
az = C.azimuths(:);
%deltaAz = abs(diff(az(1:2)));
azEst = az(pIdx(idx(1:nEst))); % + deltaAz * delta;

