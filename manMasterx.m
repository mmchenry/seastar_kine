function manMasterx(Orientation,SSnum,seqnum)
% Acquisition of sea star kinematics
% Orientation - indicates the wall orientation ('h','v','u')
% SSnum       - number of seastar
% seqnum      - sequence number


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
    
    good(1).path  =  [Orientation filesep SSnum filesep 'canon']; 
    good(1).fName = ['s' seqnum];
    
    clear Orientation SSnum seqnum
    
else
    % Add sequences to this list for analysis ---------
    %good(i).path  =  ['Upside-down' filesep 'SS41' filesep 'canon'];
    %good(i).path  =  ['Horizontal' filesep 'SS37' filesep 'canon'];
    good(1).path  =  ['Horizontal' filesep 'SS38' filesep 'canon'];
    %good(i).fName = 's03';
    %good(i).fName = 's04';
    good(1).fName = 's01';
    
end
    

%% Code execution

% Copies over portions of movies as image sequences
do.choose_dur = 0;

% Interactively select initial conditions
do.initialConditions = 0;
  
% Whether to track the body centroid
do.Centroids = 0;

% Track body rotation
do.bodyRotation = 0;

% Manual tracking of tube feet
do.manTracking = 1;

% Manual tracking of arm tips
do.armTracking = 1;

% Put together manual tracking, centroid, and rotation data
do.bundleData = 1;

% Make a movie of the data overlaid onto a video 
do.makeDataMovie = 0;


%% General parameters

% Invert image
imInvert = 1;

% Number of points to define the region of interest
numroipts = 200;

% Runs code on a single sequence, defined by seq_path
run_single = 1;


%% Manage paths 

paths = givePaths;
   

%dataPath = '/Volumes/GoogleDrive/My Drive/Projects/Andres sea stars/Up-side down SS41';
%kinePath = '/Volumes/GoogleDrive/My Drive/Projects/Andres sea stars/Up-side down SS41';
camName = 'canon';

%addpath(kinePath)


%% List of sequences to analyze

% Vist of all video filese
cList0 = catVidfiles(paths.vid,camName);

% Define path relative to root

% Initialize index
i = 1;

% --------------

% Make cList to only include good sequences
cList = struct();
k = 1;
for k = 1:length(good)
   for j = 1:length(cList0.path)
      if strcmp(good(k).path,cList0.path{j}) &&  strcmp(good(k).fName,cList0.fName{j})
          
          % Transfer all data to cList
          cList.age(k,1)      = cList0.age(j);
          cList.indiv(k,1)    = cList0.indiv(j);
          cList.orient(k,1)   = cList0.orient(j);
          cList.vidType{k,1}  = cList0.vidType{j};
          cList.path{k,1}     = cList0.path{j};
          cList.fName{k,1}    = cList0.fName{j};
          cList.ext{k,1}      = cList0.ext{j};
          cList.calPath{k,1}  = cList0.calPath{j};
          
          k = k + 1;
          
          break
      end
   end
   
   if k==1
       error('No matching videos found in list')
   elseif length(good)~=length(cList)
       warning('At lest one of the sequences not found in the video path');

   end
   
end
    
   

clear good cList0


%% Choose durations (do.choose_dur)

if do.choose_dur
    
    selectDurations(cList,paths,do)
    
end


%% Interactive mode: Select initial conditions (do.initialConditions)

if do.initialConditions
    
    % Loop thru sequences
    for i = 1:length(cList.vidType)
        
        currDataPath = [paths.data filesep cList.path{i} filesep cList.fName{i}];
        currVidPath  = [paths.vid filesep cList.path{i} filesep cList.fName{i} cList.ext{i}];

        
        % Load video info (v)
        v = defineVidObject(currVidPath);
        
        % If analysis requested and not done already, meets yes_ana
        % criteria and the clipinfo has already been determined
        if ~isempty(dir([currDataPath filesep 'clipInfo.mat'])) && ...
                isempty(dir([currDataPath filesep 'Initial conditions.mat'])) 
                          
            % Load frame intervals ('clipInfo')
            load([currDataPath filesep 'clipInfo'])
            
            if ~isnan(clipInfo.startFrame)
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
            end
        end
        clear v currDataPath currVidPath 
    end    
    disp(' ')  
end


%% Track centroid coordinates (do.Centroids)

if do.Centroids
    
    % Loop thru sequences
    for i = 1:length(cList.vidType)
        
        currDataPath = [paths.data filesep cList.path{i} filesep cList.fName{i}];
        currVidPath  = [paths.vid filesep cList.path{i} filesep cList.fName{i} cList.ext{i}];
        
        % If analysis requested and not done already . . .
        if ~isempty(dir([currDataPath filesep 'Initial conditions.mat'])) && ...
                isempty(dir([currDataPath filesep 'Centroid.mat']))
            
            disp(['Tracking centroid: ' currVidPath])
            disp('')
            
            % Load initial conditions (iC)
            load([currDataPath filesep 'Initial conditions'])
            
            % Load frame intervals ('clipInfo')
            load([currDataPath filesep 'clipInfo'])
            
            % Load video info (v)
            v = defineVidObject(currVidPath);
            
            % Frames
            frames = clipInfo.startFrame:clipInfo.endFrame;
            
            % Region of interest for first frame
            roi0 = giveROI('define','circular',numroipts,iC.r,iC.x,iC.y);
            
            % Run tracker code for centroid
            Centroid = tracker(currVidPath,v,imInvert,'threshold translation',...
                roi0,iC.tVal,frames);
            
            % Save data
            save([currDataPath filesep 'Centroid'],'Centroid')
        end
        disp(' ')
        
        clear v currDataPath currVidPath roi0 frames iC clipInfo Centroid
    end    
end


%% Track rotation (do.bodyRotation)

if do.bodyRotation
    
    % Loop thru catalog
    for i = 1:length(cList.vidType)
        
        % Current paths
        currDataPath = [paths.data filesep cList.path{i} filesep cList.fName{i}];
        currVidPath  = [paths.vid filesep cList.path{i} filesep cList.fName{i} cList.ext{i}];
             
        if ~isempty(dir([currDataPath filesep 'Initial conditions.mat'])) % && ...
        %    strcmp(anaStatus.centroid,'approved')
            
            % Load initial conditions (iC)
            load([currDataPath filesep 'Initial conditions'])
 
            % Load centroid data (Centroid)
            load([currDataPath filesep 'Centroid'])
            
            % Load mean image data (imMean)
            %load([currDataPath filesep 'meanImageData']);  
            
            % Load video info (v)
            v = defineVidObject(currVidPath);
            
            % Path for saving temporary file
            savePath = [currDataPath filesep 'tmp_Rotation.mat'];
            
            % Run tracker for predator
            if 1 %isempty(dir([currDataPath filesep 'Rotation.mat']))
                disp(' ')
                disp(['-- Predator rotation :' currDataPath])
                
                Rotation = rotationTracker(currVidPath,v,imInvert,'advanced rotation',...
                    iC.r,Centroid,Centroid.frames,savePath);
                
                % Save rotation data
                if length(Rotation.rot_ang)
                    % Save data
                    save([currDataPath filesep 'Rotation'],'Rotation')
                    
                    % Delete temporary file
                    delete(savePath)
                else
                   warning('Rotation tracking incomplete.')  
                end
            end      
                
        end
      
        clear v currDataPath currVidPath anaStatus Rotation
    end
end


%% Manual tracking mode (do.manTracking)

if do.manTracking
    
    % Loop thru sequences
    for i = 1:length(cList.vidType)
        
        % Current paths
        currDataPath = [paths.data filesep cList.path{i} filesep cList.fName{i}];
        currVidPath  = [paths.vid filesep cList.path{i} filesep cList.fName{i} cList.ext{i}];
        
        % Load frame intervals ('clipInfo')
        load([currDataPath filesep 'clipInfo'])
        
        % Load video info (v)
        v = defineVidObject(currVidPath);
        
        % Define frames vector
        %frames = v.UserData.FirstFrame:v.UserData.LastFrame;
        % Frames
        frames = clipInfo.startFrame:clipInfo.endFrame;
        
        % Load initial conditions (iC)
        load([currDataPath filesep 'Initial conditions'])
        %
        %         % Region of interest for first frame
        %         roi0 = giveROI('define','circular',numroipts,iC.r,iC.x,iC.y);
        %
        %         % Load trasnformation data (S)
        %         load([currDataPath filesep 'Transformation.mat'])
        
        savePath = [currDataPath filesep 'Manual tracking.mat'];
        
        if isempty(dir([currDataPath filesep 'Manual tracking.mat']))
            H = [];
        else
            load(savePath,'-mat');
        end
        
        % Get coordinates via interactive mode
        H = videoGUI(currVidPath,v,frames,0,'simple',iC.r,[0 1 0],H,savePath);
        
        % Save data
        v = H.v;
        H = rmfield(H,'v');
        save(savePath,'H')
        
        clear roi0 currDataPath currVidPath Rotation S
    end
end


%% Track tips of arms (do.armTracking)

if do.armTracking
    
    % Number of frames to analyze
    numFrames = 5;
    
    % Number of arms on sea star
    numArms = 6;
    
    % Loop thru sequences
    for i = 1:length(cList.vidType)
        

        % Current paths
        currDataPath = [paths.data filesep cList.path{i} filesep cList.fName{i}];
        currVidPath  = [paths.vid filesep cList.path{i} filesep cList.fName{i} cList.ext{i}];
        
        % Path for  data
        savePath = [currDataPath filesep 'Arm coords.mat'];
        
        if isempty(dir(savePath))
            
            % Load frame intervals ('clipInfo')
            load([currDataPath filesep 'clipInfo'])
            
            % Load video info (v)
            v = defineVidObject(currVidPath);
            
            % Define frames vector
            %frames = v.UserData.FirstFrame:v.UserData.LastFrame;
            
            % Frames
            frames = round(linspace(clipInfo.startFrame,clipInfo.endFrame,numFrames));
            
            % Load initial conditions (iC)
            load([currDataPath filesep 'Initial conditions'])

            % Get coordinates via interactive mode
            A = videoGUIarms(currVidPath,v,frames,0,'simple',iC.r,[0 1 0],savePath,numArms);
            
            % Default 
%             dataPass = 1;
%             
%             % Check for complete data
%             if length(A.arm)~=numArms    
%                 dataPass = 0;
%             else
%                 % Loop thru arms
%                 for j = 1:numArms
%                     % Any nans?
%                     if max(isnan(A.arm(j).x))==1
%                         dataPass = 0;
%                         break
%                     end
%                 end
%             end
%             
%             % Save data
%             if dataPass            
%                 % Save data
%                 save(savePath,'A')
%             else
%                 warning('Arm data not saved because data incomplete')
%             end
            
            clear roi0 currDataPath currVidPath Rotation S
        end
    end
end


%% Bundle all data

if do.bundleData
    
    % Loop thru sequences
    for i = 1:length(cList.vidType)
        
        disp(['Bundling data for ' cList.path{i} filesep cList.fName{i}])
        
        % Current paths
        currDataPath = [paths.data filesep cList.path{i} filesep cList.fName{i}];
        currVidPath  = [paths.vid filesep cList.path{i} filesep cList.fName{i} cList.ext{i}];
        
        % Load initial conditions (iC)
        load([currDataPath filesep 'Initial conditions'])
        
        % Load centroid data (Centroid)
        load([currDataPath filesep 'Centroid.mat'])
        
        % Load rotation data (Rotation)
        load([currDataPath filesep 'Rotation.mat'])
        
        % Load manual tracking data (H)
        load([currDataPath filesep 'Manual tracking.mat'],'-mat')
        
        % Load arm data (A)
        load([currDataPath filesep 'Arm coords.mat'],'A')
        
        % Bundle the data into S
        S = bundleData(iC,Centroid,Rotation,H,A);
        
        % Transfer sequence data
        S.orient    = cList.orient;
        S.indiv     = cList.indiv(i);
        S.paths.vid   = cList.path{i};
        S.vidName   = cList.fName{i};
        S.vidExt    = cList.ext{i};
        S.calPath   = cList.calPath{i};
        
        % Save bundled data
        save([currDataPath filesep 'Bundled Data'],'S')
        
        % Clear variables
        clear iC Centroid Rotation H A
    end   
end


%% Make movie of data

if do.makeDataMovie
    
    % Output video path
    outVidPath   = [paths.vid filesep 'Data movies'];
    
    % Whether to make the figure visible
    imVis = 1;
    
    % Loop thru sequences
    for i = 1:length(cList.vidType)
          
        % Current paths
        currDataPath = [paths.data filesep cList.path{i} filesep cList.fName{i}];
        currVidPath  = [paths.vid filesep cList.path{i} filesep cList.fName{i} cList.ext{i}];
        
        % Load data bundle ('S')
        load([currDataPath filesep 'Bundled Data.mat'])
        
        % Load video info (v)
        v = defineVidObject(currVidPath);
        
        % Make filename for video 
        ssnum = ['0' num2str(S.indiv)]; 
        fName = [S.orient 'SS' ssnum(end-1:end) '_' S.vidName];
        
        % Short version
        %makeDataMovie(currVidPath,v,S,[outVidPath filesep fName],imVis,'two view skip')
        
        % Long version
        makeDataMovie(currVidPath,v,S,[outVidPath filesep fName],imVis,'two view')
        
    end
    
end




function selectDurations(cList,paths,do)
% Interactively prompts to select a duration for analysis

% Loop thru sequences
for i = 1:length(cList.vidType)

    if 1

        currDataPath = [paths.data filesep cList.path{i} filesep cList.fName{i}];

        % Check for path in data dir
        if isempty(dir(currDataPath))
            % Make directory, if not there
            mkdir(currDataPath);
        end

        if isempty(dir([currDataPath filesep 'clipInfo.mat']))

            disp(['selectDurations: ' currDataPath]);
            disp('')
            
            getDefaults = 1;

            % Make figure
            f = figure;

            % Loop for multiple tries at frame numbers
            while true

                if getDefaults
                    % Current path
                    currVidPath = [paths.vid filesep cList.path{i} filesep ...
                        cList.fName{i} cList.ext{i}];

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

            save([currDataPath filesep 'clipInfo'],'clipInfo')
        end
    end
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




function ptsL = G2L(cOrigin,theta,ptsG)
% Coordinate transformation from global to local coordinates
% cOrgin    - Coordinate origin (3x1)
% theta     - Azimuth angle of local system wrt global
% ptsG      - Coordinates in the global FOR (3xn)

if (size(cOrigin,2) > size(cOrigin,1)) || ...
   ((length(ptsG)>3) && (size(ptsG,1) > size(ptsG,2)))  
    error('All points should be arrange in column vectors');
end

% Output nans, if nans are input
if sum(~isnan(ptsG))==0
    ptsL = ptsG;
    
% Otherwise . . .
else
    
    % Remove nans
    idx = ~isnan(ptsG(1,:));
    ptsG = ptsG(:,idx);
    
    % Add z-dimension
    cOrigin = [cOrigin; 0];
    ptsG    = [ptsG; zeros(1,size(ptsG,2))];
    
    % Axes defined
    xaxis = [cos(theta) sin(theta) 0]';
    yaxis = [cos(theta+pi/2) sin(theta+pi/2) 0]';
    zaxis = [0 0 1]';
    
    % Rotation matrix
    S = [xaxis yaxis zaxis];
    
    % Local coordinate
    ptsL = global2localcoord(ptsG,'rr',cOrigin,S);
    
    % Remove z-dimension
    ptsL = ptsL(1:2,:);
end

function ptsG = L2G(cOrigin,theta,ptsL)
% Coordinate transformation from global to local coordinates
% cOrgin    - Coordinate origin (3x1)
% theta     - Azimuth angle of local system wrt global
% ptsG      - Coordinates in the global FOR (3xn)

if (size(cOrigin,2) > size(cOrigin,1)) || ...
   ((length(ptsG)>3) && (size(ptsG,1) > size(ptsG,2)))  
    error('All points should be arrange in column vectors');
end

% Output nans, if nans are input
if sum(~isnan(ptsG))==0
    ptsL = ptsG;
    
% Otherwise . . .
else
    
    % Remove nans
    idx = ~isnan(ptsG(1,:));
    ptsL = ptsL(:,idx);
    
    % Add z-dimension
    cOrigin = [cOrigin; 0];
    ptsL    = [ptsL; zeros(1,size(ptsL,2))];
    
    % Axes defined
    xaxis = [cos(theta) sin(theta) 0]';
    yaxis = [cos(theta+pi/2) sin(theta+pi/2) 0]';
    zaxis = [0 0 1]';
    
    % Rotation matrix
    S = [xaxis yaxis zaxis];
    
    % Local coordinate
    ptsG = local2globalcoord(ptsL,'rr',cOrigin,S);
    
    % Remove z-dimension
    ptsG = ptsG(1:2,:);
end
