function makePointVideo

vPath = '/Volumes/GoogleDrive/My Drive/Projects/Andres sea stars/Up-side down SS41';
vName = 's04.MOV';

 % Load video info (v)
v = defineVidObject([vPath filesep vName]);

% Define frames vector
frames = v.UserData.FirstFrame:v.UserData.LastFrame;

% Load initial conditions (iC)
load([vPath filesep 'Initial conditions'])

% Get H
load([vPath filesep 'Manual tracking'])

for i = 1:length(frames)
    % Current frame
    cFrame = frames(i);
    
    % Get all data
    [allBase,allTip,currBase,currTip] = getFootData(H,cFrame);
    
    % Read input video frame
    frame = getFrame(H.vid_path,H.v,cFrame,H.imInvert,'rgb');
    
    imshow(frame);
    
    ddd=3;
end




function [allBase,allTip,currBase,currTip] = getFootData(H,cFrame)
% Returns data for visualization

    % Initialize containers
    allBase.x  = [];   allBase.y  = [];
    allTip.x   = [];   allTip.y   = [];

    % Loop thru feet, prior to current
    for i = 1:length(H.ft)-1
        
        % If there is a coordinate for current frame . . .
        if ~isnan(H.ft(i).xBase(cFrame))
            
            % Add coordinates
            allBase.x   = [allBase.x; H.ft(i).xBase(cFrame)];
            allBase.y   = [allBase.y; H.ft(i).yBase(cFrame)];
        end
        
        % If there is a coordinate for current frame . . .
        if ~isnan(H.ft(i).xTip(cFrame))
            
            % Add coordinates
            allTip.x    = [allTip.x; H.ft(i).xTip(cFrame)];
            allTip.y    = [allTip.y; H.ft(i).yTip(cFrame)];
        end
    end
    
    % Current base point
    currBase.x    = H.ft(end).xBase(cFrame);
    currBase.y    = H.ft(end).yBase(cFrame); 
    
    % Current tip point
    currTip.x    = H.ft(end).xTip(cFrame);
    currTip.y    = H.ft(end).yTip(cFrame); 

