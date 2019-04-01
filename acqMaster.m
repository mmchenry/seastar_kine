function acqMaster
% Acquisition of sea star kinematics

%% Code execution

% Create movies of the Centroid movies for review
do.MakeCentroidMovie = 1;

% Make movie to evaluate centroid and rotation tracking
do.MakeRotationMovies = 1;

% Make movie to evaluate centroid and rotation tracking
do.MakeFootMovie = 1;

% Re-run the rotation anlysis from the beginning 
reRunRotation = 0;

% Visualize frames to survey all steps of the analysis
do.anaSurvey = 0;


%% Verify deleting data

if reRunRotation==1
    bName = questdlg('You are about to delete exiting rotation data!',...
        'Warning!!','Proceed','Cancel','Proceed');
    if ~strcmp(bName,'Proceed')
        return
    else
        delete([currDataPath filesep 'Body.mat'])
    end
end


%% General parameters

% Invert image
imInvert = 1;

% Number of points to define the region of interest
numroipts = 200;

% Runs code on a single sequence, defined by seq_path
run_single = 1;


%% Manage paths (need to modify for new PC)
   
if ~isempty(dir(['/Users/mmchenry/Documents/Matlab code']))

    % Path to root dir of video (CSULB project, external drive)
    vidPath = '/Users/mmchenry/Documents/Video/Chip sea stars/prelim video';
    
    % Location of video frames
    %vidFramePath = '/Users/mmchenry/Documents/Video/Chip sea stars/prelim video/video frames';
    
    % Path to root of data
    dataPath = '/Users/mmchenry/Documents/Projects/Chip sea stars/prelim data';

elseif isfolder('C:\Users\tpo\Documents\seastar_kine')
    
    % Path to root dir of video (CSULB project, external drive)
    vidPath = 'C:\Users\tpo\Documents\Video\Chip sea stars\prelim video';

    % Path to root of data
    dataPath = 'C:\Users\tpo\Documents\Chip sea star data\prelim data';
    
else
    
    error('Do not recognize computer')
    
end


%% Catalog sequence videos

%cList = catVidfiles(vidPath);

cList.fName = 'S004_S001_T007';
cList.ext   = 'MOV';
cList.movtype = 'mov';
cList.path = '';

% Paths for current sequence
currDataPath = [dataPath filesep cList.path filesep cList.fName];
currVidPath  = [vidPath filesep cList.path filesep cList.fName '.' cList.ext];

% Check video path 
if ~isfile(currVidPath)
    error(['Video file does not exist at ' currVidPath]);
end

% Load video info (v)
v = defineVidObject(currVidPath);


%% Select duration of analysis

if ~isfile([currDataPath filesep 'clipInfo.mat'])
     
    % Run selecting duration
    clipInfo = selectDurations(currDataPath,currVidPath);
    
    % Save data
    save([currDataPath filesep 'clipInfo'],'clipInfo')   
    
else
    
    % Load frame intervals ('clipInfo')
    load([currDataPath filesep 'clipInfo'])
    
end


%% Interactive mode: Select initial conditions

if ~isfile([currDataPath filesep 'Initial conditions.mat']) && ...
        ~isnan(clipInfo.startFrame)
    
    disp(['Initial conditions: ' currVidPath])
    disp('')
    
    % First image
    im = getFrame(currVidPath,v,v.UserData.FirstFrame,imInvert,'gray');
    
    % Initial position
    disp(' ')
    disp('Select animal to be tracked')
    [x,y] = imInteract(im,'points',1);
    
    % Threshold
    disp(' ')
    disp('Select threshold')
    tVal = imInteract(im,'threshold');
    
    % Radius
    disp(' ')
    disp('Select roi radius')
    r = imInteract(im,'radius',x,y);
    
    % Store data
    iC.x       = x;
    iC.y       = y;
    iC.tVal    = tVal;
    iC.r       = r;
    iC.useMean = 0;
    
    % Save data
    save([currDataPath filesep 'Initial conditions'],'iC')
    
    clear im useMean im0Mean im0NoMean x y tVal r iC
    disp(' ')
    
else
    
    % Load initial conditions (iC)
    load([currDataPath filesep 'Initial conditions'])
    
end


%% Track centroid coordinates

if ~isfile([currDataPath filesep 'Centroid.mat'])
    
    disp(['Tracking centroid: ' currVidPath])
    disp('')

    % Frames
    frames = clipInfo.startFrame:clipInfo.endFrame;
    
    % Number of frames to visualize with surveyData
    numVis = 16;
    
    % Region of interest for first frame
    roi0 = giveROI('define','circular',numroipts,iC.r,iC.x,iC.y);
    
    % Run tracker code for centroid
    Centroid = tracker(currVidPath,v,imInvert,'threshold translation',...
        roi0,iC.tVal,frames);

   % Define roi for each frame
    for i = 1:length(Centroid.x)
       Centroid.roi(i) =  giveROI('define','circular',numroipts,...
           iC.r,Centroid.x(i),Centroid.y(i));     
    end
    
    % Save data
    save([currDataPath filesep 'Centroid'],'Centroid')
    
    % Visualize a bunch of frames to check results
    surveyData(currVidPath,v,imInvert,'Centroid tracking',Centroid,iC,numVis);
     
    disp(' ')
    
    clear  roi0 frames Centroid
    
else
    
     % Load centroid data (Centroid)
    load([currDataPath filesep 'Centroid.mat'])
    
end



%% Generate centroid movie for review

if do.MakeCentroidMovie && ...
   ~isfile([currDataPath filesep 'centroid movie.mat'])
    
    % Visualize steps in analysis
    imVis = 1;

    % Name of movie file
    movFile = 'Centroid tracking';

    % Frames
    frames = clipInfo.startFrame:clipInfo.endFrame;
    
    % Region of interest for first frame
    roi0 = giveROI('define','circular',numroipts,iC.r,iC.x,iC.y);
    
    % Create coordinate transformation structure
    %S = defineSystem2d('roi',roi0,Centroid);
    
    % Define roi for each frame
    for i = 1:length(Centroid.x)
       Centroid.roi(i) =  giveROI('define','circular',numroipts,...
           iC.r,Centroid.x(i),Centroid.y(i));     
    end
    
    % Update status
    disp(' '); disp(['Making Centroid Movie: ' currVidPath]); disp(' ')
    
    % Make movie
    aniData(currVidPath,v,currDataPath,movFile,imInvert,...
                'Centroid tracking',Centroid,imVis,iC);
    
    clear iC Centroid clipInfo S M mov roi0 frames imVis
end


%% Track rotation
    
if ~isfile([currDataPath filesep 'Body.mat'])
    
    % Visualize steps
    visSteps = 0;
    
    % Downsample frames for image registration (speed up processing)
    dSample = 1;

    % Numver of frames to visualize
    numVis = 16;
    
    % Run tracker for predator
    trackRotation(currVidPath,v,currDataPath,'advanced',...
        dSample,visSteps,reRunRotation,imInvert);
    
    % Load data ('Body')
    load([currDataPath filesep 'Body.mat'])
    
    % Apply post-processing
    Body = rotationPostProcess(Body);
    
    % Visualize a bunch of frames to check results
    surveyData(currVidPath,v,imInvert,'Centroid & Rotation',Body,iC,numVis);
    
    clear dSample visSteps 
    
else
    
    % Load body kinematics (Body)
    load([currDataPath filesep 'Body.mat'])
    
end


%% Apply post-processing (Body data)

% Apply post-processing
Body = rotationPostProcess(Body);

% Save
%save([currDataPath filesep 'Body, post.mat'],'Body')


%% Review rotation: make movies

if do.MakeRotationMovies
    
    imVis = 1;

    % File name of movie to be created
    fName = 'Centroid and rotation';
    
    disp(' ')
    disp(['Making Pred Rotation Movie: ' currVidPath])
    disp(' ')
    
    % Make movie
    aniData(currVidPath,v,currDataPath,fName,imInvert,...
        'Centroid & Rotation',Body,imVis);
    
    clear M mov Rotation S iC currDataPath currVidPath
    
end


%% Track feet

% Downsample
dSample = 0;

% Visualize steps
visSteps = 1;

% Number of frames to visualize
numVis = 16;

% Duration used for mean image (s)
meanDur = 10;

% Duration use for stream image (s)
streakDur = 0.5;

% Number of frames used for each mean image
numMean = 10;

% Parameters for blobs describing the tube feet
blobParam.tVal = 0.3;
blobParam.areaMin = 100;
blobParam.areaMax = 2000; 
blobParam.AR_max  = 6;

% Mean image duration in frames
meanDr_fr = round(meanDur * v.FrameRate);

% Streak image duration in frames
strakDr_fr = round(streakDur * v.FrameRate);

% Set intervals for mean frames
interval_fr = [min(Body.frames):meanDr_fr:max(Body.frames)]';

% Extend duration of last interval
interval_fr = [interval_fr(1:end-1); max(Body.frames)];

%time = (Body.frames-Body.frame(1))./v.FrameRate;

% Create mean images over regular interval of movie -----------
if ~isfile([currDataPath filesep 'mean_roi.mat'])
    
    % Loop thru intervals
    for i = 1:(length(interval_fr)-1)
       
        % Update status
        disp(['INTERVAL ' num2str(i) ' of ' num2str(length(interval_fr)-1)])
        
        % Interval of mean image
        roiM(i).frStart = interval_fr(i);
        roiM(i).frEnd   = interval_fr(i+1);
        
        % Range of frames
        dFrame = round(range([roiM(i).frStart  roiM(i).frEnd])./numMean);
        
        % Frame numbers
        roiM(i).frames = [roiM(i).frStart:dFrame:roiM(i).frEnd]';
        
        % Calc mean image
        [imMean,imStd] = motionImage(currVidPath,v,'mean roi','none',imInvert,...
                        Body,dSample,roiM(i).frames);
                    
        % Store
        roiM(i).im    = imMean;
        roiM(i).imStd = imStd;
        
        clear imMean imStd                   
    end

    % Save data
    save([currDataPath filesep 'mean_roi'],'roiM');
    
% Otherwise
else
    
    % Load mean images (roiM)
    load([currDataPath filesep 'mean_roi.mat'])
end


% Loop thru frames, create mask that excludes stationary objects ---------
if ~isfile([currDataPath filesep 'blobs.mat'])
        
    % Run blob analysis
    B = anaBlobs(currVidPath,v,'G&L props',Body,blobParam,imInvert,...
        dSample,roiM,visSteps);
    
    % Save blob data
    save([currDataPath filesep 'blobs'],'B');
      
end


% Produce data for motion images ------------------------------------------
if ~isfile([currDataPath filesep 'motion image data.mat'])

     % Load blob data (B)
    load([currDataPath filesep 'blobs']) 
    
    % Produce image stack
    imStack = motionImage(currVidPath,v,'mask static',Body.frames,B,Body.frames,imInvert);
 
    % Save images to be analyzed
    save([currDataPath filesep 'motion image data'],'-v7.3','imStack');
end


% Loop thru frames, track feet --------------------------------------------
if ~isfile([currDataPath filesep 'foot blobs.mat'])
    
    % Load blob data (B)
    load([currDataPath filesep 'blobs'])  
    
    % Load imStack
    load([currDataPath filesep 'motion image data'])
    
    B_ft = anaBlobs(currVidPath,v,'filter motion',B,strakDr_fr,Body,blobParam,...
        visSteps,imInvert,imStack);
    
    % Save blob data
    save([currDataPath filesep 'foot blobs'],'B_ft');
    
    % Visualize a bunch of frames to check results
    surveyData(currVidPath,v,0,'Feet',Body,B_ft,numVis);
    
else
    
    % Load data of feet (B_ft)
    load([currDataPath filesep 'foot blobs'])
    
end

clear dSample blobParam visSteps meanDr_fr interval_fr streakDur numMean numVis



%% Make movie of feet

if do.MakeFootMovie 
    
    imVis = 0;
    
    imInvert = 0;

    % File name of movie to be created
    fName = 'Foot tracking, global';
    
    % Load video info (v)                added by CG
    v = defineVidObject(currVidPath);
    
    disp(' ')
    disp(['Making Foot tracking Movie: ' fName])
    disp(' ')
    
    % Make movie
%     aniData(currVidPath,v,currDataPath,fName,imInvert,...
%         'Feet',Body,imVis,B_ft);

    aniData(currVidPath,v,currDataPath,fName,imInvert,...
        'Global feet',Body,imVis,B_ft);
end


%% Visualize frames from all steps of the analysis

if do.anaSurvey
    
    % Number of frames to visualize
    numVis = 16;
    
    % Centroid tracking
    surveyData(currVidPath,v,imInvert,'Centroid tracking',Centroid,iC,numVis);
    
    % Rotation
    surveyData(currVidPath,v,imInvert,'Centroid & Rotation',Body,iC,numVis);
    
    % Feet
    surveyData(currVidPath,v,0,'Feet',Body,B_ft,numVis);
    
    clear numVis
    
end





function Body = rotationPostProcess(Body)

% Loop thru frames
for i = 1:length(Body.Rotation.tform)
    
    % Scaling factor for downsampling
    Body.Rotation.roi(i).imFactor = 300./Body.Rotation.roi(i).rect(3);
    
    % Compensate tform for downsampling
    Body.Rotation.tform(i).T(3,1:2) = Body.Rotation.tform(i).T(3,1:2) ./ ...
                                      Body.Rotation.roi(i).imFactor;
    
    % Current region of interest
    roi = Body.Rotation.roi(i);
    
    % Current tform
    tform = Body.Rotation.tform(i);
    
    % And inverse
    Body.Rotation.tformInv(i) = invert(tform);  
    
    % Get corrected body center point in roi
    [xCntr,yCntr] = transformPointsForward(invert(tform),roi.r,roi.r);
    
    % Get corrected origin of roi
    [xOr,yOr] = transformPointsForward(invert(tform),0,0);
    
    % Other referece point to look at rotation
    [xOther,yOther] = transformPointsForward(invert(tform),roi.r,2*roi.r);
    
    % Angular rotation up to this point
    Body.Rotation.rot_ang(i,1)  = atan2(tform.T(1,2),tform.T(1,1)) * 180/pi;
    
    % Corrected body center point in global FOR
    Body.xCntr(i,1)  = Body.x(i)-roi.r+xCntr;
    Body.yCntr(i,1)  = Body.y(i)-roi.r+yCntr;
    Body.xOther(i,1) = Body.x(i)-roi.r+xOther;
    Body.yOther(i,1) = Body.y(i)-roi.r+yOther;
    
    % Adjust roi coordinates
    Body.Rotation.roi(i).rect(1) = Body.x(i)-roi.r+round(xOr);
    Body.Rotation.roi(i).rect(2) = Body.y(i)-roi.r+round(yOr);
    Body.Rotation.roi(i).xCntr   = Body.xCntr(i,1);
    Body.Rotation.roi(i).yCntr   = Body.yCntr(i,1);
    Body.Rotation.roi(i).xPerimG = Body.Rotation.roi(i).xPerimG + xOr;
    Body.Rotation.roi(i).yPerimG = Body.Rotation.roi(i).yPerimG + yOr; 
    
    
end


ttt = 3;
    

function clipInfo = selectDurations(currDataPath,currVidPath)
% Interactively prompts to select a duration for analysis

% Loop thru sequences

% Check for path in data dir
if isempty(dir(currDataPath))
    % Make directory, if not there
    mkdir(currDataPath);
end

if isempty(dir([currDataPath filesep 'clipInfo.mat']))
    
    getDefaults = 1;
    
    % Make figure
    f = figure;
    
    % Loop for multiple tries at frame numbers
    while true
        
        if getDefaults
            % Current path
%             currVidPath = [vidPath filesep cList.path{i} filesep ...
%                 cList.fName{i} cList.ext{i}];
            
            % Current video object
            v = defineVidObject(currVidPath,'JPG');
            
            firstFrame = v.UserData.FirstFrame;
            lastFrame  = v.UserData.LastFrame;
        end
        
        im1 = getFrame(currVidPath,v,firstFrame);
        im2 = getFrame(currVidPath,v,lastFrame);
        
        % Plot candidate frames
        subplot(1,2,1)
        imshow(im1,'InitialMag','fit')
        title('First frame')
        
        subplot(1,2,2)
        imshow(im2,'InitialMag','fit')
        title('Last frame')
        
        % Ask for approval on quality
        an = questdlg('Is this video worth analyzing?','','Yes',...
            'No','Cancel','Yes');
        
        if strcmp(an,'Yes')
            % Do nothing
            
        elseif strcmp(an,'No')
            
            firstFrame = nan;
            lastFrame  = nan;
            
            break
            
        else
            return
        end
        
        
        % Ask for approval on frames
        an = questdlg('Is the sea star in both frames?','','Yes',...
            'No','Cancel','Yes');
        
        if strcmp(an,'Yes')
            
            close(f)
            break
            
        elseif strcmp(an,'No')
            
            % Prompt for frame intervals
            [firstFrame,lastFrame] = getFrameNum(firstFrame,lastFrame,...
                im1,im2);
            getDefaults = 0;
            
        else
            return
        end
    end
    
    clipInfo.startFrame = firstFrame;
    clipInfo.endFrame   = lastFrame;
    
    %save([currDataPath filesep 'clipInfo'],'clipInfo')
end


function [firstFrame,lastFrame] = getFrameNum(firstFrame,lastFrame,im1,im2)

prompt={'Start frame num:','Last frame num:'};
name='Choose clip duration';
numlines=1;
defaultanswer={num2str(firstFrame),num2str(lastFrame)};

answer = inputdlg(prompt,name,numlines,defaultanswer);

if isempty(answer)
    return
end

firstFrame   = str2num(answer{1});
lastFrame    = str2num(answer{2});


function cList = catVidfiles(vidPath)


a1 = dir(vidPath);

n = 1;

novid = 1;

% Loop thru for adult and juvenile directories
for i = 1:length(a1)
    
    if a1(i).isdir && strcmp(a1(i).name,'Adults')      
        currAge = 'a';       
    elseif a1(i).isdir && strcmp(a1(i).name,'Juveniles')       
        currAge = 'j';        
    else       
        currAge = [];        
    end
    
   
    
    if ~isempty(currAge)
        
        a2 = dir([vidPath filesep a1(i).name]);
        
        % Loop thru oriention directories
        for j = 1:length(a2)
            
            if a2(j).isdir && strcmp(a2(j).name,'Horizontal')
                
                currOrient = 'h';
                
            elseif a2(j).isdir && strcmp(a2(j).name,'Vertical')
                
                currOrient = 'v';
                
            elseif a2(j).isdir && ...
                    (strcmp(a2(j).name,'Upside-down') || strcmp(a2(j).name,'UpsideDown'))
                
                currOrient = 'u';
                
            else
                currOrient = [];
            end
            
            
            if ~isempty(currOrient)
                
                a3 = dir([vidPath filesep a1(i).name filesep a2(j).name]);
                
                % Loop trhu individual directories
                for k = 1:length(a3)
                    
                    if a3(k).isdir && length(a3(k).name)>2 &&...
                            strcmp(a3(k).name(1:2),'SS')
                        
                        indivNum = str2num(a3(k).name(3:4));
                        
                        % Directory contents for individual
                        a4 = dir([vidPath filesep a1(i).name filesep ...
                            a2(j).name filesep a3(k).name]);

                        m = 1;cal.Type=[];
                        
                        % Loop thru sequences for the indiviual
                        for l = 1:length(a4)
                            
                            % If video is a calibration
                            if length(a4(l).name) > 10 && ...
                                    (strcmp(a4(l).name(1:11),'calibration') || ...
                                     strcmp(a4(l).name(1:11),'Calibration'))
                                     
                                if a4(l).isdir
                                    cal.Type = 'image';
                                    
                                elseif ~a4(l).isdir && strcmp(a4(l).name(end-2:end),'MP4')
                                    cal.Type = 'mp4';
                                    
                                elseif ~a4(l).isdir && strcmp(a4(l).name(end-2:end),'MOV')
                                    cal.Type = 'mov';
                                    
                                else
                                    error('no match for calibration type');
                                end
                                
                                % Store calibration path
                                cal.Path = [a1(i).name filesep ...
                                    a2(j).name filesep a3(k).name ...
                                    filesep a4(l).name];
                                
                                % If video is an MOV file . . .
                            elseif ~a4(l).isdir && length(a4(l).name)>3 && ...
                                    strcmp(a4(l).name(end-2:end),'MOV')
                                
                                movPath{m} = [a1(i).name filesep ...
                                    a2(j).name filesep a3(k).name ...
                                    filesep a4(l).name];
                                movType{m} = 'mov';
                                
                                m = m + 1;
                                
                            elseif ~a4(l).isdir && length(a4(l).name)>3 && ...
                                    strcmp(a4(l).name(end-2:end),'MP4')
                                
                                movPath{m} = [a1(i).name filesep ...
                                    a2(j).name filesep a3(k).name ...
                                    filesep a4(l).name];
                                movType{m} = 'mp4';
                                
                                m = m + 1;    
                                % If video is an image sequence
                            elseif a4(l).isdir && length(a4(l).name)>3 && ...
                                    strcmp(a4(l).name(1:4),'time')
                                
                                movType{m} = 'image seq';
                                movPath{m} = [a1(i).name filesep ...
                                    a2(j).name filesep a3(k).name ...
                                    filesep a4(l).name];
                            end
                        end

                        % Loop trhu sequences
                        for f = 1:length(movType)                            
                            
                            % Extract file parts
                            [pathstr,name,ext] = fileparts(movPath{f});
                            
                            
                            % Store info on inidvidual
                            cList.age(n,1)    = currAge;
                            cList.indiv(n,1)  = indivNum;
                            cList.orient(n,1) = currOrient;
                            
                            cList.vidType{n} = movType{f};
                            cList.path{n}    = pathstr;
                            cList.fName{n}   = name;
                            cList.ext{n}     = ext;
                            
                            novid = 0;
                            
                            if isempty(cal.Type)
                                cList.calPath{n} = [];
                                %disp(' ')
                                warning(['No calibration for: ' movPath{f}]);
                            else
                                cList.calPath{n} = cal.Path;
                            end
                            
                            % Advance sequence index
                            n = n + 1;                            
                        end                     
                    end
                end              
            end
        end     
    end   
end

if novid==1
    error('No video files found')
end