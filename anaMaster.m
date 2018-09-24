function anaMaster(Orientation,SSnum,seqname)
% Runs analysis of kinematic data
% Direction       - 'h', 'u',  or 'v'
% SSnum             - Individual number of sea star.
% seqnum            - number of sequenc


if nargin<3
    error('You need to specific sequence info')
elseif 
    
end


% Assign defaults
if nargin<3
    %seqname = 's03';
    seqname = 's01';
    if nargin < 2
        %SS = 'SS37';
        SS = 'SS38';
        if nargin < 1
            orientation = 'Horizontal';
            %orientation = 'Upside-down';
        end
    end
end


%% Manage paths



if ~isempty(dir(['/Users/mmchenry/Documents/Matlab code']))
    
    % Path to kineBox
    %kinePath = '/Users/mmchenry/Documents/Matlab code/kineBox';
    kinePath = '/Users/mmchenry/Documents/Matlab code/kineBox_old';
    
    % Path to root dir of video (CSULB project, external drive)
    vidPath = '/Volumes/GoogleDrive/My Drive/Flow backup/Shared_flow';
    %vidPath = '/Volumes/Video/Sea stars/CSULB/Raw video';
    
    % Location of video frames
    %vidFramePath = '/Volumes/Video/Sea stars/CSULB test/Video frames';
    
    % Path to root of data
    dataPath = '/Volumes/GoogleDrive/My Drive/Projects/Andres sea stars/Kinematics';
    
% Line to assign single vids    
elseif ~isempty(dir(['C:\Program Files\MATLAB\R2016a']))
    
    %vidPath = '\\flow.local\shared\Sea stars';
    vidPath = 'C:\Users\andres\Documents\SS Assign';
    %special vid path
    %vidpath=
    % dataPath = '\\flow.local\andres\SS Assign\CSULB data'; %% by CG
    dataPath = 'C:\Users\andres\Documents\dataPath';
    
    kinePath = 'C:\Users\andres\Documents\GitPath\kineBox';
else
    error('Do not recognize computer')
    
end


%% Display data

% Camera view
camView = 'canon';

imInvert = 0;

skipPlay = 10;

% Current paths
currVidPath   = [vidPath filesep orientation filesep SS filesep camView filesep seqname '.MOV'];
currDataPath  = [dataPath  filesep orientation filesep SS filesep camView filesep seqname];

% Load video info (v)
v = defineVidObject(currVidPath,'MOV');

% Load data (H)
warning off
load([currDataPath filesep 'Manual tracking.mat'],'-mat')
warning on

% Define frames vector
%frames = v.UserData.FirstFrame:v.UserData.LastFrame;

% Load initial conditions (iC)
load([currDataPath filesep 'Initial conditions'])

% % Cue up coordinates
% xBase = nan(length(H.ft(1).xBase),length(H.ft));
% yBase = xBase;
% xTip  = xBase;
% yTip  = xBase;

% Create matrices of 
for i = 1:(length(H.ft)-2)
    xBase(:,i)  = H.ft(i).xBase;
    yBase(:,i)  = H.ft(i).yBase;
    xTip(:,i)   = H.ft(i).xTip;
    yTip(:,i)   = H.ft(i).yTip;
end

% Data frames
tmp = ~isnan(xTip);
dFrames = H.frames(sum(tmp,2)>0);
clear tmp


% Loop trhu frames
for i = min(dFrames):max(dFrames)
    
    cFrame = H.frames(i);
    
    %  frame
    im = getFrame(currVidPath,v,cFrame,imInvert,'gray');
    
    
    imshow(im,'InitialMag','fit')
    hold on
    
    for j = 1:length(H.ft)
        plot(H.ft(j).xTip(cFrame),H.ft(j).yTip(cFrame),'or')
    end
    
    title(['Frame ' num2str(cFrame)])
    
    hold off
    
    pause(0.1)
end
ttt= 3;



