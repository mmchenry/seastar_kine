function manMaster3D(Orientation,SSnum,seqnum)
% Acquisition of sea star kinematics
% Orientation - indicates the wall orientation ('h','v','u')
% SSnum       - number of seastar
% seqnum      - sequence number
%
% Code can only be run for sequences that have completed manMasterx



%% Deal with inputs

if nargin >=3
    
    if strcmp(Orientation,'h')
        Orientation = 'Horizontal';
        
    elseif strcmp(Orientation,'v')
        Orientation = 'Vertical';
        
    elseif strcmp(Orientation,'u')    
        Orientation = 'Upside-down';
        
    else
        error('Do not recorgnize the orientation requested');
    end
    
    % Make string out of seastar number
    SSnum = ['0' num2str(SSnum)];
    SSnum = ['SS' SSnum(end-1:end)];
    
    % Make string out of sequence number
    seqnum = ['0' num2str(seqnum)];
    seqnum = seqnum(end-1:end);
    
    good(1).path_c   =  [Orientation filesep SSnum filesep 'canon']; 
    good(1).path_s    =  [Orientation filesep SSnum filesep 'sony']; 
    good(1).fName        = ['s' seqnum];
    
    clear Orientation SSnum seqnum
    
else
    % Add sequences to this list for analysis ---------
    good(1).path_c  =  ['Horizontal' filesep 'SS38' filesep 'canon'];
    good(1).path_s  =  ['Horizontal' filesep 'SS38' filesep 'sony'];
    good(1).fName = 's01';
    
end
    
clear Orientation SSnum seqnum


%% Define paths and check for data

paths = givePaths;

% Vist of all video files
cList0 = catVidfiles(paths.vid,'canon');
cList1 = catVidfiles(paths.vid,'sony');

% Check for videos
if isempty(cList0)
    error(['No canon videos in ' paths.vid]);
elseif isempty(cList1)
    error(['No sony videos in ' paths.vid]);
end

% Make cList to only include good sequences
cList = struct();

% Loop thru canon sequences
for i = 1:length(cList0.path)
    
    % If path and file match the 'good' sequence . . .
    if strcmp(good.path_c,cList0.path{i}) &&  strcmp(good.fName,cList0.fName{i})

        % Transfer all data to cList
        cList.age       = cList0.age(i);
        cList.indiv     = cList0.indiv(i);
        cList.orient    = cList0.orient(i);
        cList.vidType   = cList0.vidType{i};
        cList.path_c    = cList0.path{i};
        cList.fName     = cList0.fName{i};
        cList.ext       = cList0.ext{i};
        cList.calPath   = cList0.calPath{i};
        
        break        
    end  
end

% Look for matches in sony view
for i = 1:length(cList1.path)
    
    % If match on sony path . . .
    if strcmp(good.path_s,cList1.path{i}) && strcmp(good.fName,cList1.fName{i})
        
        % Store path
        cList.path_s     = cList1.path{i};

        break
    end
end


% Check for sony video
if ~isfield(cList,'path_s')
    error('Matching video file missing for sony camera');
end

% Define path for canon data
paths.canonData = [paths.data filesep cList.path_c filesep cList.fName];

% Define path for sony data
paths.sonyData = [paths.data filesep cList.path_s filesep cList.fName];

% Path for canon video
paths.canonVid = [paths.vid filesep cList.path_c filesep cList.fName cList.ext];

% Path for sony video
paths.sonyVid = [paths.vid filesep cList.path_s filesep cList.fName cList.ext];

% Path for sony calibration
paths.sonyCalVid = [paths.vid filesep cList.path_s filesep 'scale' cList.ext];

% Path for sony dot movie
paths.sonyDotVid = [paths.vid filesep cList.path_s filesep 'dot' cList.ext];

% Path for canon dot movie
paths.canonDotVid = [paths.vid filesep cList.path_c filesep 'dot' cList.ext];

% Check canon bundled data
if isempty(dir([paths.canonData filesep 'Bundled Data.mat']))
    error(['Bundled Data.mat file missing -- you first need ' ...
        ' to complete manMasterx'])
end

% Check for videos
if isempty(dir([paths.sonyVid]))
    error('Sony video missing')
elseif isempty(dir([paths.sonyCalVid]))
    error('Sony calibration video missing')
elseif isempty(dir([paths.sonyDotVid]))
    error('Sony dot video missing')
elseif isempty(dir([paths.canonDotVid]))
    error('Canon dot video missing')
end

% Make data directory, if none exists for sony
if isempty(dir(paths.sonyData))
    mkdir(paths.sonyData);
end

% Check for canon calibration data
if isempty(dir([paths.canonData filesep 'Calibration.mat']))
    error(['Calibration.mat file missing -- you first need ' ...
        ' to run calibration in manMasterx'])
end

% Load data bundle for canon view ('S')
load([paths.canonData filesep 'Bundled Data.mat'])

% Clear variables  
clear good cList0 iSony i j cList0 cList1 good currDataPath


%% Spatial calibration

% Number of times to repeat calibration
nRepeat = 3;

% Load video info (v)
v = defineVidObject(paths.sonyCalVid);

% If analysis not done already
if isempty(dir([paths.sonyData filesep 'Calibration.mat']))
    
    % First frame
    im = getFrame(paths.sonyCalVid,v,1,0,'gray');
    
    for j = 1:nRepeat
        
        % Calibration length
        pixLen = imInteract(im,'length',3);
        
        % Prompt for real length
        answer = inputdlg({'Actual length (mm)'},'',1,{''});
        
        % Log constant
        calConst(j,1) = (str2num(answer{1})/1000)/pixLen;
        
        % Store
        cal_sony.const = mean(calConst);
        
        clear answer pixLen
    end
    
    % Save data
    save([paths.sonyData filesep 'Calibration'],'cal_sony')
      
else
    
   % Load calibration data for sony
   load([paths.sonyData filesep 'Calibration'],'cal_sony') 
end

clear nRepeat i v im  

% Load calibration data for canon
load([paths.canonData filesep 'Calibration'],'cal_canon') 


%% Common dot calibration

% If analysis not done already
if isempty(dir([paths.sonyData filesep 'Dot Calibration.mat']))
    
    % Load video info
    v_s = defineVidObject(paths.sonyDotVid);
    v_c = defineVidObject(paths.canonDotVid);
    
    % First frames
    im_c = getFrame(paths.canonDotVid,v_c,1,0,'gray');
    im_s = getFrame(paths.sonyDotVid,v_s,1,0,'gray');
    
    disp(' '); disp('Select common point in both views'); disp(' '); 
    
    % Dot position 
    [dot.sony.x,dot.sony.y]     = imInteract(im_s,'points',1);
    [dot.canon.x,dot.canon.y]   = imInteract(im_c,'points',1);
     
    % Save
    save([paths.sonyData filesep 'Dot Calibration.mat'],'dot')
    
    clear im_c im_s v_s v_c
else 
    % Load 'dot' data
    load([paths.sonyData filesep 'Dot Calibration.mat'],'dot')
end


%% Sync cameras wrt time

% If analysis not done already
if isempty(dir([paths.sonyData filesep 'cam sync.mat']))
    
    % Load video info
    v_s = defineVidObject(paths.sonyVid);
    v_c = defineVidObject(paths.canonVid);
    
    cam.info  = ' delay is for sony relative to canon';
    cam.delay = audio_sync(paths.canonVid,paths.sonyVid);
    
    % Save
    save([paths.sonyData filesep 'cam sync.mat'],'cam')
    
    clear im_c im_s v_s v_c
else 
    % Load 'cam' data
    load([paths.sonyData filesep 'cam sync.mat'],'cam')
end


%% Determine corresponding frames in two videos

% Load video info
v_s = defineVidObject(paths.sonyVid);
v_c = defineVidObject(paths.canonVid);

% Delay in frames
frame_delay = round(cam.delay * v_s.FrameRate);

if frame_delay>0
    
    % Canon frames
    frames_c = S.frames;
    
    % Sony frames
    frames_s = [v_s.UserData.FirstFrame:v_s.UserData.LastFrame]  + frame_delay;

else
    
    % Sony frames
    frames_s = [v_s.UserData.FirstFrame:v_s.UserData.LastFrame];
    
    % Canon frames
    frames_c = S.frames + frame_delay;
    
end

% Number of frames
nFrames = min([length(frames_s) length(frames_c)]);

% Trim ends
frames_s = frames_s(1:nFrames);
frames_c = frames_c(1:nFrames);

% Trim negative values
if max(frames_c<0)
    idx = frames_c>0;
    frames_s = frames_s(idx);
    frames_c = frames_c(idx);
end

% Trim negative values
if max(frames_s<0)
    idx = frames_s>0;
    frames_s = frames_s(idx);
    frames_c = frames_c(idx);
end


clear v_s v_c nFrames frame_delay


%% Gather together camera info

camIn = cam; clear cam

% Store sony data
cam.s.delay    = camIn.delay;
cam.s.xDot     = dot.sony.x;
cam.s.yDot     = dot.sony.y;
cam.s.cal      = cal_sony.const;
cam.s.frames   = frames_s;

% Store canon data
cam.c.xDot   = dot.canon.x;
cam.c.yDot   = dot.canon.y;
cam.c.cal    = cal_canon.const;
cam.c.frames = frames_c;

clear dot cal_sony cal_canon camIn frames_s frames_c


%% Measure the margins of the body from both views

% If analysis not done already
if isempty(dir([paths.sonyData filesep 'margins.mat']))
    
    % Number of measuremeents to take
    nMeas = 7;
    
    % Vector of corresponding frames
    frames_s = round(linspace(min(cam.s.frames),max(cam.s.frames),nMeas));
    frames_c = round(linspace(min(cam.c.frames),max(cam.c.frames),nMeas));
    
    % Load video info
    v_s = defineVidObject(paths.sonyVid);
    v_c = defineVidObject(paths.canonVid);
    
    for i = 1:nMeas

        % First frames
        im_c = getFrame(paths.canonVid,v_c,frames_c(i),0,'gray');
        im_s = getFrame(paths.sonyVid,v_s,frames_s(i),0,'gray');
        
        % Title text
        t_txt = ['Measurement ' num2str(i) ' of ' num2str(nMeas)];
        
        disp(' '); disp('Select L-R margins of body in both views'); disp(' ');
        
        % Dot position
        tmp_s     = imInteract(im_s,'margins',t_txt);
        tmp_c     = imInteract(im_c,'margins',t_txt);
        
        % Store
        m.frame_c(i,1) = frames_c(i);
        m.frame_s(i,1) = frames_s(i);
        m.c(i,:)       = [min(tmp_c) max(tmp_c)];
        m.s(i,:)       = [min(tmp_s) max(tmp_s)];
    end
    
    % Save
    save([paths.sonyData filesep 'margins.mat'],'m')
    
    clear im_c im_s v_s v_c
else 
    % Load margin ('m') data
    load([paths.sonyData filesep 'margins.mat'],'m')
end


%% Analyze changes in relative size of body

% Linear fit for relative body size in sony view, relative to canon
%c = polyfit(m.frame_s,hypot(m.s(:,1),m.s(:,2))./hypot(m.c(:,1),m.c(:,2)),1);

% Margin fits for sony
cam.s.mFrames = m.frame_s;
cam.s.mLcoef  = polyfit(m.frame_s,m.s(:,1),1);
cam.s.mRcoef  = polyfit(m.frame_s,m.s(:,2),1);

% Margin fits for canon
cam.c.mFrames = m.frame_c;
cam.c.mLcoef  = polyfit(m.frame_c,m.c(:,1),1);
cam.c.mRcoef  = polyfit(m.frame_c,m.c(:,2),1);

% cam.s.mL = m.s(1,1);
% cam.s.mR = m.s(1,2);
% 
% % Margin fits for canon
% cam.c.mL = m.c(1,1);
% cam.c.mR = m.c(1,2);

clear c m


%% Eyespot selection

% Load video info
v_s = defineVidObject(paths.sonyVid);
v_c = defineVidObject(paths.canonVid);

% Path for data saving
savePath = [paths.sonyData filesep 'eye_data.mat'];

% Load or create H
if ~isempty(dir(savePath))
    load(savePath,'H')
else
    H = [];
end

% Get coordinates via interactive mode
videoGUI3D(paths.sonyVid,v_s,paths.canonVid,v_c,S,cam,'arms',savePath,H);




