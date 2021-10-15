function [prob_AFN_F] = compute_posteriors_GMM(sig, fsHz, C)
%compute_posteriors_GMM   GMM-based localisation module.
%
%USAGE
%  [prob_AFN_F] = compute_posteriors_GMM(sig, fsHz, C)
%
%INPUT ARGUMENTS
%          sig : binaural input signal [nSamples x 2]
%         fsHz : sampling frequency in Hertz
%            C : classifier structure used for localisation
%
%OUTPUT ARGUMENTS
%   prob_AFN_F : integrated posteriors
% 
%
% Ning Ma, 29 Jan 2015
% n.ma@sheffield.ac.uk
%   
%

 
% AFE processing
dObj = dataObject(sig, fsHz);
mObj = manager(dObj, C.AFE_request_mix, C.AFE_param);

mObj.processSignal();

% Work out which features to use
features = strsplit(C.ftrType, '-');
bUseFeatures = false(4, 1); % ITD, ILD, CC, IC
for n = 1:length(features)
    switch lower(features{n})
        case 'itd'
            itd = dObj.itd{1}.Data(:);
            bUseFeatures(1) = true;
            nFrames = size(itd,1);
        case 'ild'
            ild = dObj.ild{1}.Data(:);
            bUseFeatures(2) = true;
            nFrames = size(ild,1);
        case 'cc'
            cc = dObj.crosscorrelation{1}.Data(:);
            % Use only -1ms to 1ms
            idx = ceil(size(cc,3)/2);
            mlag = ceil(fsHz/1000);
            cc = cc(:,:,idx-mlag:idx+mlag);
            bUseFeatures(3) = true;
            nFrames = size(cc,1);
        case 'ic'
            ic = dObj.ic{1}.Data(:);
            bUseFeatures(4) = true;
            nFrames = size(ic,1);
    end
end

nChannels = C.AFE_param.fb_nChannels;

post = zeros(nFrames, C.nAzimuths, nChannels);
for c = 1:nChannels
    
    testFeatures = [];
    if bUseFeatures(1)
        testFeatures = [testFeatures itd(:,c)];
    end
    if bUseFeatures(2)
        testFeatures = [testFeatures ild(:,c)];
    end
    if bUseFeatures(3)
        testFeatures = [testFeatures squeeze(cc(:,c,:))];
    end
    if bUseFeatures(4)
        testFeatures = [testFeatures ic(:,c)];
    end
    
    % Normalise features
    testFeatures = testFeatures - repmat(C.featNorm{c}(1,:),[size(testFeatures,1) 1]);
    testFeatures = testFeatures ./ sqrt(repmat(C.featNorm{c}(2,:),[size(testFeatures,1) 1]));

    % Determine performance of trained classifier using test material
    post(:,:,c) = classifyGMM(testFeatures,C.gmmFinal{c});
end

% % Select T-F units based on energy
% rate = (dObj.ratemap{1}.Data(:) + dObj.ratemap{2}.Data(:)) ./ 2;
% rate = log(max(rate,eps));
% energyMask = rate > (median(rate(:)*1.1));
% frameIdx = sum(energyMask,2) > 2;

% Integrate probabilities across all frequency channels
%prob_AF = exp(squeeze(nansum(log(post(frameIdx,:,:)),3)));
prob_AF = exp(squeeze(nansum(log(post),3)));

% Normalise each frame such that probabilities sum up to one
prob_AFN = prob_AF ./ repmat(sum(prob_AF,2),[1 C.nAzimuths]);
%prob_AFN = transpose(prob_AF ./ repmat(sum(prob_AF,2),[1 C.nAzimuths]));
        
% Integrate probabilities acrros all time frames
prob_AFN_F = nanmean(prob_AFN, 1);

% subplot(211); plot(C.azimuths, prob_AFN_F); title('Energetic frames')
% subplot(212); plot(C.azimuths, nanmean(prob_AFN,1)); title('All frames')

