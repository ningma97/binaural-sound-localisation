function runme
% Following these steps to run everything at a single computer node, but
% this wil take a while to complete.
%
% Following the qsub_xxx steps to run everything on a computer cluster.
%

%% Define parameters
%
% 'MCT-DIFFUSE' for multicondition training
% 'CLEAN' for clean training
preset = 'MCT-DIFFUSE';

% 'cc-ild' to use cross-correlation (CC) and ILD (34-d at 16kHz FS)
% 'itd-ild' to use ITD and ILD (2-d)

%% GMM
model = 'GMM';
featType = 'itd-ild';
f2_evaluate_RealBRIR(model, false, preset, featType);
f2_evaluate_RealBRIR(model, true, preset, featType);

%% DNN
model = 'DNN';
featType = 'ild-cc';
f2_evaluate_RealBRIR(model, false, preset, featType);
f2_evaluate_RealBRIR(model, true, preset, featType);


% 
% azErrPer = squeeze(mean(f3_extractResults('GCC')));
% fprintf('Room\tAnech\tRoom A\tRoom B\tRoom C\tRoom D\tAvg.\n');
% [nRooms, nSpks] = size(azErrPer);
% for s=1:nSpks
%     fprintf('%d-spk', s);
%     for r=1:nRooms
%         fprintf('\t%.1f', azErrPer(r,s)*100);
%     end
%     fprintf('\t%.1f\n', mean(azErrPer(:,s))*100);
% end
% fprintf('\n');
% for r=1:nRooms
%     for s=1:nSpks
%         err = azErrPer(r,s)*100;
%         if err > 99.9
%             fprintf('&100   ');
%         else
%             fprintf('&%.1f  ', err);
%         end
%     end
%     fprintf(' ');
% end
% fprintf('\\\\\n');

        
