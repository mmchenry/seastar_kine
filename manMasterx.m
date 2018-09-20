function manMasterx
% Acquisition of sea star kinematics

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


%% General parameters

% Invert image
imInvert = 1;

% Number of points to define the region of interest
numroipts = 200;

% Runs code on a single sequence, defined by seq_path
run_single = 1;


%% Manage paths (need to modify for new PC)
   
if ~isempty(dir(['/Users/mmchenry/Documents/Matlab code']))
    
    % Path to kineBox
    %kinePath = '/Users/mmchenry/Documents/Matlab code/kineBox';
    kinePath = '/Users/mmchenry/Documents/Matlab code/kineBox_old';
    
    % Path to root dir of video (CSULB project, external drive)
    %vidPath = '/Volumes/Video/Sea stars/CSULB test/Raw video';
    %vidPath = '/Volumes/Video/Sea stars/CSULB/Raw video';
    vidPath = '/Volumes/GoogleDrive/My Drive/Flow backup/Shared_flow';
    %vidPath = '/Users/mmchenry/Documents/Video/Sea stars'
    
    % Location of video frames
    vidFramePath = '/Volumes/Video/Sea stars/CSULB test/Video frames';
    
    % Path to root of data
    %dataPath = '/Users/mmchenry/Documents/Projects/Seastars/CSULB data';
  % remember to undo %%  to run all vids
  dataPath = '/Volumes/GoogleDrive/My Drive/Projects/Andres sea stars/Kinematics';
  
 elseif ~isempty(dir(['C:\Program Files\MATLAB\R2016a']))
    
    %vidPath = '\\flow.local\shared\Sea stars';
    %vidPath = 'C:\Users\andres\Documents\Sea stars';
    %special vid path
    %vidpath=
    % dataPath = '\\flow.local\andres\Sea stars\CSULB data';
    
    % kinePath = 'C:\Users\andres\Documents\GitPath\kineBox';
    
% Line to assign single vids    
elseif ~isempty(dir(['C:\Program Files\MATLAB\R2016a']))
    
    %vidPath = '\\flow.local\shared\Sea stars';
    vidPath = 'C:\Users\andres\Documents\SS Assign';
    %special vid path
    %vidpath=
    % dataPath = '\\flow.local\andres\SS Assign\CSULB data'; %% by CG
    dataPath = 'C:\Users\andres\Documents\dataPath'
    
    kinePath = 'C:\Users\andres\Documents\GitPath\kineBox';
   
    %%%Hunk of code added to Andres' Mac%%%
elseif ~isempty(dir(['/Users/andrescarrillo/seastar_kine']))
    
    vidPath = '/Volumes/AChd2TB';
    %vidPath = '/Volumes/Samsung_T3';
    %special vid path
    %vidpath=
    % dataPath = '\\flow.local\andres\SS Assign\CSULB data'; %% by CG
    dataPath = '/Users/andrescarrillo/seastar_kine/dataPath';
    
    kinePath = '/Users/andrescarrillo/seastar_kine';
    
    camName = 'canon';
else
    error('Do not recognize computer')
    
end

%dataPath = '/Volumes/GoogleDrive/My Drive/Projects/Andres sea stars/Up-side down SS41';
%kinePath = '/Volumes/GoogleDrive/My Drive/Projects/Andres sea stars/Up-side down SS41';
camName = 'canon';

%addpath(kinePath)


%% List of sequences to analyze

% Vist of all video filese
cList0 = catVidfiles(vidPath,camName);

% Define path relative to root

% Initialize index
i = 1;

% Add sequences to this list for analysis ---------
%good(i).path  =  ['Horizontal' filesep 'SS38' filesep 'canon']; 
good(i).path  =  ['Horizontal' filesep 'SS37' filesep 'canon']; 
good(i).fName = 's03';
i = i + 1;

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


%% Produce video stills

if do.choose_dur
    
    selectDurations(cList,dataPath,vidPath,do)
    
end


%% Interactive mode: Select initial conditions

if do.initialConditions
    
    % Loop thru sequences
    for i = 1:length(cList.vidType)
        
        currDataPath = [dataPath filesep cList.path{i} filesep cList.fName{i}];
        currVidPath  = [vidPath filesep cList.path{i} filesep cList.fName{i} cList.ext{i}];

        
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


%% Track centroid coordinates

if do.Centroids
    
    % Loop thru sequences
    for i = 1:length(cList.vidType)
        
        currDataPath = [dataPath filesep cList.path{i} filesep cList.fName{i}];
        currVidPath  = [vidPath filesep cList.path{i} filesep cList.fName{i} cList.ext{i}];
        
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


%% Track rotation

if do.bodyRotation
    
    % Loop thru catalog
    for i = 1:length(cList.vidType)
        
        % Current paths
        currDataPath = [dataPath filesep cList.path{i} filesep cList.fName{i}];
        currVidPath  = [vidPath filesep cList.path{i} filesep cList.fName{i} cList.ext{i}];
        
        % Check approval status
        if ~isempty(dir([currDataPath filesep 'anaStatus.mat']))
            load([currDataPath filesep 'anaStatus.mat'])
        else
            warning(['You have not yet approved centroid tracking for' ...
                currDataPath])
        end
              
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
            
            % Run tracker for predator
            if isempty(dir([currDataPath filesep 'Rotation.mat']))
                disp(' ')
                disp(['-- Predator rotation :' currDataPath])
                
                Rotation = tracker(currVidPath,v,imInvert,'advanced rotation',...
                    iC.r,Centroid,Centroid.frames);
                
                % Save rotation data
                save([currDataPath filesep 'Rotation'],'Rotation')
            end      
                
        end
      
        clear  v currDataPath currVidPath anaStatus Rotation
    end
end


%% Bundle Centroid and Rotation data

if 0
    % Loop thru sequences
    for i = 1:length(cList.vidType)
    
        % Current paths
        currDataPath = [dataPath filesep cList.path{i} filesep cList.fName{i}];
        currVidPath  = [vidPath filesep cList.path{i} filesep cList.fName{i} cList.ext{i}];
    
        % If centroid data there and centroid data not yet approved
         if isempty(dir([currDataPath filesep 'Transformation.mat']))
        
            % Load initial conditions (iC)
            load([currDataPath filesep 'Initial conditions'])
            
            % Region of interest for first frame
            roi0 = giveROI('define','circular',numroipts,iC.r,iC.x,iC.y);
        
            % Load centroid data (Centroid)
            load([currDataPath filesep 'Centroid.mat'])
        
            % Load rotation data (Rotation)
            load([currDataPath filesep 'Rotation.mat'])
        
            % Create coordinate transformation structure
            S = defineSystem2d('roi',roi0,Centroid,Rotation);
        
            % Save transformation data
            save([currDataPath filesep 'Transformation'],'S')
         end  
        clear roi0 currDataPath currVidPath Rotation S
    end
end


%% Manual tracking mode (do.manTracking)

% Loop thru sequences
for i = 1:length(cList.vidType)
    
    % Current paths
    currDataPath = [dataPath filesep cList.path{i} filesep cList.fName{i}];
    currVidPath  = [vidPath filesep cList.path{i} filesep cList.fName{i} cList.ext{i}];
    
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
    
    if isempty(dir([currDataPath filesep 'Manual tracking.mat']))
        H = [];
    else
        load([currDataPath filesep 'Manual tracking.mat']);
    end
    
    % Get coordinates via interactive mode
    H = videoGUI(currVidPath,v,frames,0,'simple',iC.r,[0 1 0],H);
    
    % Save data
    save([currDataPath filesep 'Manual tracking.mat'],'H');

        
    clear roi0 currDataPath currVidPath Rotation S
end





function selectDurations(cList,dataPath,vidPath,do)
% Interactively prompts to select a duration for analysis

% Loop thru sequences
for i = 1:length(cList.vidType)

    if 1

        currDataPath = [dataPath filesep cList.path{i} filesep cList.fName{i}];

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
                    currVidPath = [vidPath filesep cList.path{i} filesep ...
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


function cList = catVidfiles(vidPath,camName)
% Catalogs all video files in directory tree


a1 = dir(vidPath);


n = 1;
novid = 1;

% Loop thru oriention directories
for i = 1:length(a1)
    
    if a1(i).isdir && strcmp(a1(i).name,'Horizontal')
        
        currOrient = 'h';
        
    elseif a1(i).isdir && strcmp(a1(i).name,'Vertical')
        
        currOrient = 'v';
        
    elseif a1(i).isdir && ...
            (strcmp(a1(i).name,'Upside-down') || strcmp(a1(i).name,'UpsideDown'))
        
        currOrient = 'u';
        
    else
        currOrient = [];
    end
    
    
    if ~isempty(currOrient)
        
        a2 = dir([vidPath filesep a1(i).name filesep 'SS*']);

        % Loop thru juvenile directories
        for j = 1:length(a2)
            currAge = 'j';
            
            a3 = dir([vidPath filesep a1(i).name filesep a2(j).name filesep ...
                camName filesep '*.MOV']);
         
            indivNum = str2num(a2(j).name(3:4));
            
            cal.Type = [];
            
            m = 1;
            
            % Catalog, only if movies present
            if length(a3)>0
                
                movPath = [];
                movType = [];
                
                % Loop trhu individual directories
                for k = 1:length(a3)
                    
                    %                 if a3(k).isdir && length(a3(k).name)>2 &&...
                    %                         strcmp(a3(k).name(1:2),'SS')
                    %
                    %
                    %                     % Directory contents for individual
                    %                     a4 = dir([vidPath filesep a1(j).name filesep ...
                    %                         a2(j).name filesep a3(k).name]);
                    %
                    %                     m = 1;cal.Type=[];
                    
                    % Loop thru sequences for the indiviual
                    
                    % Extract file parts
                    %[pathstr,name,ext] = fileparts(movPath{f});
                    
                    % If video is a calibration
                    if length(a3(k).name) > 7 && ...
                            (strcmp(a3(k).name,'scale.MOV') || ...
                            strcmp(a3(k).name,'Scale.MOV') || ...
                            strcmp(a3(k).name,'Calibration.MOV') ||...
                            strcmp(a3(k).name,'calibration.MOV'))
                        
                        if strcmp(a3(k).name(end-2:end),'MOV')
                            cal.Type = 'mov';
                            
                        else
                            error('no match for calibration type');
                        end
                        
                        % Store calibration path
                        cal.Path = [a1(i).name filesep ...
                            a2(j).name filesep a3(k).name];
                        
                        % If video is an MOV file . . .
                    elseif strcmp(a3(k).name(1:2),'s0') && strcmp(a3(k).name(end-2:end),'MOV')
                        
                        movPath{m} = [a1(i).name filesep ...
                            a2(j).name filesep camName filesep a3(k).name ];
                        movType{m} = 'mov';
                        
                        m = m + 1;
                    end
                end

                if ~isempty(movPath)
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
                
                
            else
                warning(['No videos in: ' vidPath filesep a1(i).name ...
                         filesep a2(j).name filesep camName ])
            end
        end
    end
end





if novid==1
    warning('No video files found')
end
