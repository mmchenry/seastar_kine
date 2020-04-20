function acqMaster(dataPath,vidPath,action)
% Acquisition of sea star kinematics, only from a bottom view of the tube
% feet on glass


%% Code execution defaults

% Re-runs acquisition
do.rerunAcq = 0;

% Re-runs analysis (not including acquisition)
do.rerunAna = 0;

% Create movies of the Centroid movies for review
do.MakeCentroidMovie = 0;

% Make movie to evaluate centroid and rotation tracking
do.MakeRotationMovies = 0;

% Make movie of foot tracking
do.MakeFootMovie = 0;

% Make movie of foot tracking for a presentation
do.MakeFootMoviePretty = 0;

% Make video for analysis by DeepLabCut
d.MakeDeepMovie = 0;

% Make movie to evaluate centroid and rotation tracking, after
% post-processing
do.MakeFootMoviePost = 0;

% Re-run the rotation anlysis from the beginning 
reRunRotation = 0;

% Visualize frames to survey all steps of the analysis
do.anaSurvey = 0;

% Visualize steps of analysis executed
visSteps = 0;

% Surveys the data to show tracking over a fixed number of video stills
visFootTracking = 0;


%% 

% Run MakeDeepMovie, if requested
if strcmp(action,'run acq')
    
    % Do nothing
    
elseif strcmp(action,'rerun acq')
    
    % Re-runs acquisition
    do.rerunAcq = 1;
    
elseif strcmp(action,'deep movie')
    
    d.MakeDeepMovie = 1;
    
else
    error(['Do not recognize ' action])
end


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

% Number of frames to visualize with surveyData
numVis = 16;

% Duration used for mean image (s)
meanDur = 10;

% Duration use for stream image (s)
streakDur = 0.5;

% Number of frames used for each mean image
numMean = 10;

% Parameters for blobs describing the tube feet
blobParam.tVal    = 0.3;
blobParam.areaMin = 50;
blobParam.areaMax = 2000; 
blobParam.AR_max  = 6;

maxSize = 350;


%% Manage paths 
   
paths = givePaths;


%% Catalog sequence videos

%cList = catVidfiles(paths.vid);

% %cList.fName = 'S004_S001_T007';
% if nargin<1
%     %cList.fName = 'SS001_S001_T013';
%     %cList.fName = 'STUDIO0_S009_S001_T020';
%     cList.fName = ['weights' filesep 'STUDIO0_S009_S001_T009'];
% else
%     
%     cList.fName = fileName;
% end
% cList.ext   = 'MOV';
% cList.movtype = 'mov';
% cList.path = '';

% Paths for current sequence
currDataPath = [paths.data filesep dataPath];
currVidPath  = [paths.vid filesep vidPath];

% Check video path 
if ~isfile(currVidPath)
    error(['Video file does not exist at ' currVidPath]);
end

% Check data path 
if ~isfolder(currDataPath)
    error(['Data folder does not exist at ' currDataPath]);
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
    disp('Select center of mouth')
    [x,y] = imInteract(im,'points',1);
    
    % Arm tip position
    disp(' ')
    disp('Select arm tips, starting at 12:00 and going clockwise')
    [xArms,yArms] = imInteract(im,'points',5);
    
    % Check
    if length(xArms)~=5
        error('You need to select 5 arms');
    end
    
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
    iC.xArms   = xArms;
    iC.yArms   = yArms;
    iC.tVal    = tVal;
    iC.r       = r;
    iC.useMean = 0;
    
    
    % First image
    im = getFrame(currVidPath,v,v.UserData.FirstFrame,imInvert,'gray');
    
    % Binary
    bw = ~imbinarize(im,iC.tVal);
    bw = imfill(bw,'holes');
    bw = bwselect(bw,iC.x,iC.y);
    
    % Survey blobs
    props = regionprops(bw,'Area');
    
    if length(props)>1
        error('Too many blobs')
    else
        iC.area = props.Area;
    end
    
    % Save data
    save([currDataPath filesep 'Initial conditions'],'iC')
    
    clear im useMean im0Mean im0NoMean x y tVal r
    disp(' ')
    
else
    
    % Load initial conditions (iC)
    load([currDataPath filesep 'Initial conditions'])
    
end


%% Track centroid coordinates

if do.rerunAcq || ~isfile([currDataPath filesep 'Centroid.mat'])
    
    disp(['Tracking centroid: ' currVidPath])
    disp('')

    % Frames
    frames = v.UserData.FirstFrame:v.UserData.LastFrame;
    
    % Region of interest for first frame
    roi0 = giveROI('define','circular',numroipts,iC.r,iC.x,iC.y);
    
    % Run tracker code for centroid
    Centroid = trackCenter(currVidPath,v,imInvert,'threshold translation',...
        roi0,iC,frames);

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
                'Centroid tracking',Centroid,visSteps,iC);
    
    clear iC Centroid clipInfo S M mov roi0 frames 
end


%% Track rotation
    
if do.rerunAcq || ~isfile([currDataPath filesep 'Body.mat'])
    
    % Downsample frames for image registration (speed up processing)
    dSample = 1;
    
    % Run tracker for predator
    trackRotation(currVidPath,v,currDataPath,'advanced',...
        dSample,visSteps,reRunRotation,imInvert,maxSize);
    
    % Load data ('Body')
    load([currDataPath filesep 'Body.mat'])
    
    % Apply post-processing
    %Body = rotationPostProcess(Body,v,iC);
    Body = postProcess('post rotation',Body,v,iC,maxSize);
    
    % Save with post-processing
    save([currDataPath filesep 'Body.mat'],'-v7.3','Body')
    
    % Visualize a bunch of frames to check results
    surveyData(currVidPath,v,imInvert,'Centroid & Rotation',Body,iC,numVis);
    
    clear dSample 
    
else
    
    % Load body kinematics (Body)
    load([currDataPath filesep 'Body.mat'])
    
    % Make sure post-processing done
    if ~isfield(Body,'xCntr')
        % Apply post-processing
        Body = postProcess('post rotation',Body,v,iC,maxSize);
        
        % Save with post-processing
        save([currDataPath filesep 'Body.mat'],'Body')
    end
end


%% Review rotation: make movies

if do.MakeRotationMovies
    % File name of movie to be created
    fName = 'Centroid and rotation';
    
    disp(' ')
    disp(['Making Pred Rotation Movie: ' currVidPath])
    disp(' ')
    
    % Make movie
    aniData(currVidPath,v,currDataPath,fName,imInvert,...
        'Centroid & Rotation',Body,visSteps);
    
    clear M mov Rotation S iC currDataPath currVidPath   
end


%% Foot tracking step 1: Create series of mean images in local FOR

% Downsample
dSample = 0;

% Mean image duration in frames
meanDr_fr = round(meanDur * v.FrameRate);

% Streak image duration in frames
strakDr_fr = round(streakDur * v.FrameRate);

% Set intervals for mean frames
interval_fr = [min(Body.frames):meanDr_fr:max(Body.frames)]';

% Extend duration of last interval
interval_fr = [interval_fr(1:end-1); max(Body.frames)];

% Path for mean images
mPath = [currDataPath filesep 'mean_images'];

% Create mean images over regular interval of movie -----------
if do.rerunAcq || ~isfolder([currDataPath filesep 'mean_images'])
    
    % Make directory
    if ~isfolder(mPath)
        mkdir(mPath);
    end
    
    % Loop thru intervals
    for i = 1:(length(interval_fr)-1)
       
        % Update status
        disp(['INTERVAL ' num2str(i) ' of ' num2str(length(interval_fr)-1)])
        
        % Interval of mean image
        roiM.frStart = interval_fr(i);
        roiM.frEnd   = interval_fr(i+1);
        
        % Range of frames
        dFrame = round(range([roiM.frStart  roiM.frEnd])./numMean);
        
        % Frame numbers
        roiM.frames = [roiM.frStart:dFrame:roiM.frEnd]';
        
        % Calc mean image
        [imMean,imStd] = motionImage(currVidPath,v,'mean roi','none',imInvert,...
                        Body,dSample,roiM.frames,iC);
                    
        % Store
        roiM.im    = imMean;
        roiM.imStd = imStd;
        
        % Strings for starting and ending frames
        frStartStr = ['00000' num2str(roiM.frStart)];
        frStartStr = frStartStr(end-5:end);
        frEndStr   = ['00000' num2str(roiM.frEnd)];
        frEndStr   = frEndStr(end-5:end);
        
        % Current filename
        fName = ['mean_image_' frStartStr '_' frEndStr];
        
        % Save data
        save([mPath filesep fName],'roiM');
        
        clear imMean imStd fName frStartStr frEndStr roiM                   
    end

    clear dSample meanDr_fr interval_fr streakDur numMean mPath
    
end


%% Foot tracking step 2: mask for feet applied to global FOR
% Creates local mask that excludes stationary objects 

if do.rerunAcq || ~isfolder([currDataPath filesep 'blobs']) 
      
    % Downsample
    dSample = 0;
    
    % Path for mean images
    mPath = [currDataPath filesep 'mean_images'];

    % Load initial conditions (iC)
    load([currDataPath filesep 'Initial conditions'])
    
    % Path for saving blobs
    savePath = [currDataPath filesep 'blobs'];
    
    % Run blob analysis (generates B, saved in 'blobs' files)
    anaBlobs(currVidPath,v,'G&L props',Body,blobParam,imInvert,...
        dSample,mPath,visSteps,iC,savePath);   
end

% Produce data for motion images im imStack
if do.rerunAcq || ~isfolder([currDataPath filesep 'mask_static'])
    
    % Path for motion images
    motionPath = [currDataPath filesep 'mask_static'];
    
    % Path to blob data
    blobPath = [currDataPath filesep 'blobs'];
    
    % Load initial conditions (iC)
    load([currDataPath filesep 'Initial conditions'])
    
    tStart = tic;
    
    % Produce image stack (Creates series of blur images in mask_static)
    motionImageStack(currVidPath,v,'mask static',Body,blobPath, ...
                   imInvert,iC,motionPath,strakDr_fr);
    disp(' ');
    disp(['Time (min) = ' num2str(toc(tStart)/60)])
 
end


%% Foot tracking step 3: Track feet

% Loop thru frames, track feet --------------------------------------------
if 1 %do.rerunAcq || ~isfolder([currDataPath filesep 'foot_blobs'])

    % Stores linked feet in B_ft, in 'foot_blobs' folder
    anaBlobs(currVidPath,v,'filter motion',currDataPath,strakDr_fr,Body,blobParam,...
        visSteps,imInvert);
    
    % Load initial conditions (iC)
    load([currDataPath filesep 'Initial conditions'])

    % Visualize a bunch of frames to check results
    surveyData(currVidPath,v,0,'Feet',currDataPath,Body,numVis,iC);
    
    clear B B_ft imStack
       
end


%% Foot tracking step 4: Post-processing

% Distance threshold for including feet
dist_thresh = 5;

% Arm number assignment ------------------------------------------
if do.rerunAna || ~isfile([currDataPath filesep 'post- arms.mat'])
    
    % Path to data of feet (B_ft)
    Bpath = [currDataPath filesep 'foot_blobs'];
    
    % Find arm numbers and matching global data
    B2 = postProcess('package data',Body,iC,Bpath);
    
    % Find arm numbers and matching global data
    B2 = postProcess('find offset',Body,iC,Bpath,B2);

    % Add arms in global FOR to Body
    [Body,B2] = postProcess('add arms',Body,iC,B2);
    
    % Save Body and B2
    save([currDataPath filesep 'Body, post.mat'],'-v7.3','Body')
    save([currDataPath filesep 'post- arms'],'-v7.3','B2')
    
    % Add arms and trajectory coordinate systems to body
    Body = postProcess('Traj body system',Body);

    % Save
    save([currDataPath filesep 'Body, post.mat'],'-v7.3','Body')
    
    % Visual check with a random frame
    %visArms(currVidPath,v,Body,B2,200)
end

% Connect blobs across frames --------------
if do.rerunAna || ~isfile([currDataPath filesep 'post- foot refined.mat'])
    
    % Load B2
    if ~exist('B2','var')
        load([currDataPath filesep 'post- arms'])
    end
    
    % Load Body
    load([currDataPath filesep 'Body, post.mat'])
    
    % Update status
    disp('postProcess: Connecting feet across frames . . .')
    
    % Package into F structure
    F = postProcess('connect',Body,iC,B2,dist_thresh);
    
    % Add arm and trajectory coordinate systems to feet
    F = postProcess('Traj foot system',Body,F);
    
    % Remove erroneous feet
    F = postProcess('Remove bad feet',Body,F);

    % Save data
    save([currDataPath filesep 'post- foot refined'],'F')
    
    clear F 
end

% Visualize result
% surveyData(currVidPath,v,0,'Individual feet',Body,B_ft,numVis);


%% Define frames used

% Get listing of frame numbers
[a,frNums] = fileList([currDataPath filesep 'foot_blobs'],'foot_blobs');

% Index for frames
iFrames = Body.frames>=frNums(1) & Body.frames<=frNums(end);

% Save data
save([currDataPath filesep 'Frames used'],'iFrames')


%% Make movie for analyzing via DeepLabCut

if d.MakeDeepMovie
    
    ttt= 3;
    
    
end


%% Make movie of individual feet, after post-processing

if do.MakeFootMoviePost
    
    % Turn of image inversion
    imInvert = 0;
    
     % Load B2
    if ~exist('B2','var')
        load([currDataPath filesep 'post- arms'])
    end
    
    % Load Body
    load([currDataPath filesep 'Body, post.mat'])
    
    % Load F
    load([currDataPath filesep 'post- foot refined'])
    
%     % Get index of frames used
%     for i = 1:length(B2)
%         if isempty(B2(i).frIdx)
%             iFrames(i) = 0==1;
%         else
%             iFrames(i) = 1==1;
%         end
%     end
%     
%     % Save data
%     save([currDataPath filesep 'Frames used'],'iFrames')
%     
    clear B2

    % File name of movie to be created
    fName = 'Foot tracking, indivdual white';
    
    % Load video info (v)                added by CG
    v = defineVidObject(currVidPath);
    
    disp(' ')
    disp(['Making Foot tracking Movie: ' fName])
    disp(' ')

    aniData(currVidPath,v,currDataPath,fName,imInvert,...
        'Individual feet',Body,visSteps,F,iFrames,iC);
%     aniData(currVidPath,v,currDataPath,fName,imInvert,...
%         'Global feet',Body,visSteps,B_ft);

    
    clear B2 F

end


%% Make movie of feet for presentation

if do.MakeFootMoviePretty

    % Turn of image inversion
    imInvert = 0;

    % Load iFrames
    load([currDataPath filesep 'Frames used'])
    
    % Load F data
    load([currDataPath filesep 'post- foot refined']);
    
    % Load Body
    load([currDataPath filesep 'Body, post.mat'])
    
    % File name of movie to be created
    fName = 'Foot tracking, pretty';
    
    % Load video info (v)                added by CG
    v = defineVidObject(currVidPath);
    
    % Update status
    disp(' ')
    disp(['Making Foot tracking Movie: ' fName])
    disp(' ')

    % Create animation
    aniData(currVidPath,v,currDataPath,fName,imInvert,...
        'Individual feet, pretty',Body,visSteps,F,iFrames,iC);
end


%% Visualize a sampling of video frames

if visFootTracking
    
    % Numebr of frames to visualize
    nFrames = 20;

    % Load iFrames
    load([currDataPath filesep 'Frames used'])
    
    % Load F data
    load([currDataPath filesep 'post- foot refined']);
    
    % Load Body
    load([currDataPath filesep 'Body, post.mat'])
    
    % File name of movie to be created
    fName = 'Foot tracking, pretty';
    
    % Load video info (v)                added by CG
    v = defineVidObject(currVidPath);
    
    % Load video info (v)                added by CG
%     v = defineVidObject(currVidPath);
    
    % Create animation
    visTracking(currVidPath,v,currDataPath,imInvert,...
                'Basic',Body,visSteps,F,iFrames,iC,nFrames);
    
end


%% Visualize frames from all steps of the analysis

if do.anaSurvey
    
    if ~exist('Body','var')
         % Load Body
         load([currDataPath filesep 'Body, post.mat'])
    end
    
    if ~exist('B_ft','var')
        % Load data of feet (B_ft)
        load([currDataPath filesep 'foot blobs'])
    end
    
    % Number of frames to visualize
    numVis = 16;
    
    % Centroid tracking
    surveyData(currVidPath,v,imInvert,'Centroid tracking',Centroid,iC,numVis);
    
    % Rotation
    surveyData(currVidPath,v,imInvert,'Centroid & Rotation',Body,iC,numVis);
    
    % Feet
    surveyData(currVidPath,v,0,'Feet local',Body,B_ft,numVis);
    
end



function visArms(currVidPath,v,Body,B2,cFrame)

% Colormap of line colors
cmap(1,:) = [0         0.4470    0.7410];
cmap(2,:) = [0.8500    0.3250    0.0980];
cmap(3,:) = [0.9290    0.6940    0.1250];
cmap(4,:) = [0.4940    0.1840    0.5560];
cmap(5,:) = [0.4660    0.6740    0.1880];
cmap(6,:) = [0.3010    0.7450    0.9330];
cmap(7,:) = [0.6350    0.0780    0.1840];

figure

% Frame number and index
iFrame = find(Body.frames==cFrame,1,'first');

% Current whole frame
im = getFrame(currVidPath,v,cFrame,0,'gray');

imshow(im,'InitialMag','fit');
hold on

% Loop thru potential feet
for i = 1:length(B2(iFrame).L)
    
    xC = B2(iFrame).G(i).Centroid(1);
    yC = B2(iFrame).G(i).Centroid(2);
    
    aNum  = B2(iFrame).L(i).armNum;
    
    scatter(xC,yC,'SizeData',200,'MarkerEdgeColor',...
        'w','MarkerEdgeAlpha',0.5);
    hold on
    
    if ~isnan(aNum)
        
        scatter(B2(iFrame).G(i).Centroid(1),...
            B2(iFrame).G(i).Centroid(2),...
            'SizeData',200,'MarkerEdgeColor',...
            cmap(aNum,:));
        title(num2str(aNum))
    end
    ttt=3;
end

hold off


% function Body = rotationPostProcess(Body,v,iC)
% 
% % Define time
% Body.t = Body.frames'./v.FrameRate;
% 
% % Define local coordinates for arms 
% Body.xArmL = iC.xArms' - Body.x(1) + Body.Rotation.roi(1).r;
% Body.yArmL = iC.yArms' - Body.y(1) + Body.Rotation.roi(1).r;
% 
% % Loop thru frames
% for i = 1:length(Body.Rotation.tform)
%     
%     % Scaling factor for downsampling
%     Body.Rotation.roi(i).imFactor = 300./Body.Rotation.roi(i).rect(3);
%     
%     % Compensate tform for downsampling
%     Body.Rotation.tform(i).T(3,1:2) = Body.Rotation.tform(i).T(3,1:2) ./ ...
%                                       Body.Rotation.roi(i).imFactor;
%     
%     % Current region of interest
%     roi = Body.Rotation.roi(i);
%     
%     % Current tform
%     tform = Body.Rotation.tform(i);
%     
%     % And inverse
%     Body.Rotation.tformInv(i) = invert(tform);  
%    
%     % Get corrected body center point in roi
%     [xCntr,yCntr] = transformPointsForward(invert(tform),roi.r,roi.r);
%     
%     % Get corrected origin of roi
%     [xOr,yOr] = transformPointsForward(invert(tform),0,0);
%     
%     % Other referece point to look at rotation
%     [xOther,yOther] = transformPointsForward(invert(tform),roi.r,2*roi.r);
%     
%     % Angular rotation up to this point
%     Body.Rotation.rot_ang(i,1)  = atan2(tform.T(1,2),tform.T(1,1)) * 180/pi;
%     
%     % Corrected body center point in global FOR
%     Body.xCntr(i,1)  = Body.x(i)-roi.r + xCntr;
%     Body.yCntr(i,1)  = Body.y(i)-roi.r + yCntr;
%     Body.xOther(i,1) = Body.x(i)-roi.r + xOther;
%     Body.yOther(i,1) = Body.y(i)-roi.r + yOther;
% %     Body.xArmL(i,:)  = Body.x(i)-roi.r + xArm';
% %     Body.yArmL(i,:)  = Body.y(i)-roi.r + yArm';
%     
%     % Adjust roi coordinates
%     Body.Rotation.roi(i).rect(1) = Body.x(i)-roi.r+round(xOr);
%     Body.Rotation.roi(i).rect(2) = Body.y(i)-roi.r+round(yOr);
%     Body.Rotation.roi(i).xCntr   = Body.xCntr(i,1);
%     Body.Rotation.roi(i).yCntr   = Body.yCntr(i,1);
%     Body.Rotation.roi(i).xPerimG = Body.Rotation.roi(i).xPerimG + xOr;
%     Body.Rotation.roi(i).yPerimG = Body.Rotation.roi(i).yPerimG + yOr; 
%     
%     
% end


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
%             currVidPath = [paths.vid filesep cList.path{i} filesep ...
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
            
            return
            
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




function [a,frNums] = fileList(fPath,fPrefix)

frNums = [];

% File listing
a = dir([fPath filesep fPrefix '*']);

% Loop trhu files
for i = 1:length(a)
    
    % Index of separator
    iSep = find(a(i).name=='_',1,'last');
    
    % Get frame number
    a(i).frNum = str2num(a(i).name((iSep+1):end-4));
    
    frNums = [frNums; a(i).frNum];
end
