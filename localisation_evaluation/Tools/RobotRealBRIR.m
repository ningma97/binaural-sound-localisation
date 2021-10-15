classdef RobotRealBRIR < handle
    %ROBOTMCT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        headOrientation = 0
        audio                     % Source audio [nSamples x nSources]
        nTotalSamples             % Number of audio samples
        nTotalSources             % Number of sources
        fsHz                      % Sampling rate
        azIdx                     % Source reference position index
        brirs                     % binaural room impulse response
        brirScales                % BRIR scales for correcting level diff
        isFinished = false;       % Return true if end of audio reached
        sampleIndex = 0;          % Current processing index
    end
    
    methods
        function obj = RobotRealBRIR(audio, fsHz, azIdx, brirs, brirScales)
            obj.audio = audio;
            obj.fsHz = fsHz;
            obj.azIdx = azIdx;
            obj.brirs = brirs;
            if nargin > 4
                obj.brirScales = brirScales;
            else
                obj.brirScales = ones(length(brirs), 1);
            end
            [obj.nTotalSamples,obj.nTotalSources] = size(audio);
        end
        
        function rotateHead(obj, angleDeg)
            if angleDeg ~= 0
                obj.headOrientation = mod(obj.headOrientation+angleDeg, 360);
            end
        end
        
        % Returns the duration of signal in seconds if there is less than
        % durSec seconds left
        function durSec = canGetSignal(obj, durSec)
            if obj.isFinished
                durSec = 0;
            else
                durSamples = floor(durSec * obj.fsHz);
                leftSamples = obj.nTotalSamples - obj.fsHz * 10E-3 - obj.sampleIndex;
                durSec = min(durSamples, leftSamples) / obj.fsHz;
            end
        end
        
        function durSec = remainingDuration(obj)
            if obj.isFinished
                durSec = 0;
            else
                durSec = (obj.nTotalSamples - obj.sampleIndex) / obj.fsHz;
            end
        end
        
        function [sig,fsHz] = getSignal(obj, durSec)
            fsHz = obj.fsHz;
            if obj.isFinished
                sig = [];
            else
                if ~exist('durSec', 'var')
                    durSec = obj.remainingDuration;
                end
                
                durSamples = floor(durSec * obj.fsHz);
                blockSamples = (1:durSamples) + obj.sampleIndex;
                
                % Make sure the last block does not exceed nTotalSamples
                blockSamples = blockSamples(blockSamples<=obj.nTotalSamples);
                
                obj.sampleIndex = obj.sampleIndex + durSamples;
                
                % If there is less than 10 ms we stop processing
                if obj.sampleIndex >= obj.nTotalSamples - obj.fsHz * 10E-3
                    obj.isFinished = true;
                end
                
                % Direction   : Left         Front         Right
                % WP1 azimuth : 90   60   30   0   330  300  270  
                % headIdx     : 181  151  121  91  61   31   1
                if obj.headOrientation <= 90
                    headIdx = 91 + obj.headOrientation;
                elseif obj.headOrientation >= 270
                    headIdx = 91 - (360 - obj.headOrientation);
                else
                    error('Unsupported head orientation: %d', obj.headOrientation);
                end
                
                sig = [];
                brirFs = obj.brirs{1}.Data.SamplingRate;
                % Mix sources
                for n = 1 : obj.nTotalSources
                    
                    if brirFs == obj.fsHz
                        x = obj.audio(blockSamples,n);
                    else
                        x = resample(obj.audio(blockSamples,n),brirFs,obj.fsHz);
                    end

                    if isempty(sig)
                        sig = zeros(numel(x), 2);
                    end
                    
                    brirScale = obj.brirScales(obj.azIdx(n));
                    % Left ear
                    sig(:,1) = sig(:,1) + fftfilt(brirScale.*squeeze(obj.brirs{obj.azIdx(n)}.Data.IR(headIdx,1,:))', x);
                    % Right ear
                    sig(:,2) = sig(:,2) + fftfilt(brirScale.*squeeze(obj.brirs{obj.azIdx(n)}.Data.IR(headIdx,2,:))', x);
                end
                
                if brirFs ~= obj.fsHz
                    sig = resample(sig,obj.fsHz,brirFs);
                end
                
            end
        end
        
    end
    
end

