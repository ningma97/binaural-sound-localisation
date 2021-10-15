function rotateAngle = getHeadRotationAngle(azimuths, azimuthPosteriors)
%
%


% Get the most likely source position
[~,idx] = max(azimuthPosteriors);
mlAz = azimuths(idx);

% Roate the head towards the ML azimuth
rotateAngle = 20;
signs = [-1 1];
if mlAz > 180
    rotateAngle = -rotateAngle;
elseif mlAz == 0 || mlAz == 180
    rotateAngle = signs(randi(2)) * rotateAngle;
end

% % Random head rotation
% rotationAngles = [20 -20];
% rotateAngle = rotationAngles(randi(length(rotationAngles)));


% 
% if nargin < 1
%     rotateAngle = rotationAngles(randi(length(rotationAngles)));
% else
%     while true
%         rotateAngle = rotationAngles(randi(length(rotationAngles)));
%         newHO = mod(headOrientation + rotateAngle, 360);
% 
% 
%         relAzRef = mod(azRef - newHO, 360);
%         relAzRef = convertAzimuthsWP1ToSurrey(relAzRef);
% 
%         % Make sure source azimuths relative to the new head orientation stay 
%         % within [-90 90]
%         if max(abs(relAzRef)) <= 90
%             break;
%         end
%     end
% end
