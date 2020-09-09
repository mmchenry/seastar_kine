function visMaster(dataPath,vidPath,action)
% Visualizations to review acquisition of sea star kinematics

%% Check input 'action'

if ~strcmp(action,'centroid movie') && ...
   ~strcmp(action,'rotation movie') && ...
   ~strcmp(action,'feet local movie') && ...
   ~strcmp(action,'feet pretty movie') && ...
   ~strcmp(action,'feet frames') && ...
   ~strcmp(action,'feet local frames') && ...
   ~strcmp(action,'centroid frames') && ...
   ~strcmp(action,'rotation frames')

    error(['Do not recognize ' action])
end


%% General parameters

% Invert image
imInvert = 1;

numroipts = 500;

% Number of frames to visualize with surveyData
numVis = 16;

% Whether to visualize individual frames, as they are being generated
visSteps = 0;


%% Manage paths and load data
  
paths = givePaths;

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

% Load frame intervals ('clipInfo')
load([currDataPath filesep 'clipInfo'])

% Load initial conditions (iC)
load([currDataPath filesep 'Initial conditions'])


%% Generate centroid movie ('centroid movie') 

if strcmp(action,'centroid movie') 

    % Load centroid data (Centroid)
    load([currDataPath filesep 'Centroid.mat'])
    
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


%% Generate rotation movie ('rotation movie')

if strcmp(action,'rotation movie') 
    
    % File name of movie to be created
    fName = 'Centroid and rotation';
    
    disp(' ')
    disp(['Making Pred Rotation Movie: ' currVidPath])
    disp(' ')
    
    % Make movie
    aniData(currVidPath,v,currDataPath,fName,imInvert,...
        'Centroid & Rotation',Body,visSteps,iC);
    
    clear M mov Rotation S fName  
end


%% Make movie of individual feet ('feet local movie')

if strcmp(action,'feet local movie') 
    
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
    
    % Load frames used (iFrames)
    load([currDataPath filesep 'Frames used'])
    
    clear B2

    % File name of movie to be created
    fName = 'Foot tracking, indivdual white';
    
    % Load video info (v)                added by CG
    v = defineVidObject(currVidPath);
    
    disp(' ')
    disp(['Making Foot tracking Movie: ' fName])
    disp(' ')
    
    % Generate movie
    aniData(currVidPath,v,currDataPath,fName,imInvert,...
        'Individual feet, local',Body,visSteps,F,iFrames,iC);
%     aniData(currVidPath,v,currDataPath,fName,imInvert,...
%         'Global feet',Body,visSteps,B_ft);
   
    clear B2 F fName
end


%% Make movie of feet for presentation ('feet pretty movie')

if strcmp(action,'feet pretty movie') 

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


%% Visualize video stills of feet ('centroid frames')

if strcmp(action,'centroid frames') 
    
    % Visualize a bunch of frames to check results of rotation tracking
     surveyData(currVidPath,v,imInvert,'Centroid tracking',Centroid,iC,...
        numVis,v.Name);
    
end


%% Visualize video stills of feet ('rotation frames')

if strcmp(action,'rotation frames') 
    
    if ~exist('Body','var')
         % Load Body
         load([currDataPath filesep 'Body, post.mat'])
    end
    
    if ~exist('F','var')
         % Load Body
         load([currDataPath filesep 'post- foot refined'])
    end
    
    % Visualize a bunch of frames to check results of rotation tracking
    surveyData(currVidPath,v,imInvert,'Centroid & Rotation',Body,iC,numVis);
    
end


%% Visualize video stills of feet ('feet frames')

if strcmp(action,'feet frames') 
    
    % Visualize a bunch of frames to check results of foot tracking
    surveyData(currVidPath,v,imInvert,'Feet',currDataPath,Body,...
        numVis,iC,v.Name);
end


%% Visualize video stills of feet ('feet local frames')

if strcmp(action,'feet local frames') 
    
    % Visualize a bunch of frames to check results of foot tracking
    surveyData(currVidPath,v,imInvert,'Feet local',currDataPath,Body,...
        numVis,iC,v.Name);
end



%% Visualize a sampling of video frames

%TODO: Test this old code
if 0 
    
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





return

%% Generate centroid movie data (do.MakeCentroidMovie)

if do.MakeCentroidMovie
        
    % If centroid data there and centroid data not yet approved
    if isempty(dir([currDataPath filesep 'centroid movie.mat'])) && ...
        ~isempty(dir([currDataPath filesep 'Centroid.mat']))
        
        % Load initial conditions (iC)
        load([currDataPath filesep 'Initial conditions'])
        
        % Load centroid data (Centroid)
        load([currDataPath filesep 'Centroid.mat'])
        
        % Load frame intervals ('clipInfo')
        load([currDataPath filesep 'clipInfo'])
        
        % Load video info (v)
        v = defineVidObject(currVidPath);
        
        % Frames
        frames = clipInfo.startFrame:clipInfo.endFrame;
        
        % Region of interest for first frame
        roi0 = giveROI('define','circular',numroipts,iC.r,iC.x,iC.y);
        
        % Create coordinate transformation structure
        S = defineSystem2d('roi',roi0,Centroid);
        
        disp(' ')
        disp(['Making Centroid Movie: ' currVidPath])
        disp(' ')
        
        % Make movie
        M = aniData(currVidPath,v,imInvert,'Centroid tracking',S,0);
        
        mov.dataPath      = currDataPath;
        mov.currVidPath   = currVidPath;
        mov.M             = M;
        
        % Save movie data
        save([currDataPath filesep 'centroid movie'],'mov','-v7.3')
        
    elseif ~isempty(dir([currDataPath filesep 'centroid movie.mat']))
        
        disp('Centroid movie data already generated.')
        
    elseif isempty(dir([currDataPath filesep 'Centroid.mat']))
        
        error(['You have not yet run the Centroid analysis on '...
               currVidPath])       
    end
    
    clear v iC Centroid clipInfo S M mov roi0 frames
end


%% Review recent centroid tracking (do.playCentroid)

if do.playCentroid
    
    % If centroid data there and centroid data not yet approved
    if ~isempty(dir([currDataPath filesep 'centroid movie.mat']))
        
        % Load centroid movie (mov)
        load([currDataPath filesep 'centroid movie'])
        
        % Make figure window
        f = figure('units','normalized');
        
        while true
            
            movie(f,mov.M,1,20,[0 0 1 1])
            
            b = questdlg('What to you want to do?','','Replay','Close','Delete','Replay');
            
            if strcmp(b,'Delete')
                
                % Prompt for approval
                b2 = questdlg('Delete the movie data?','','Yes','No','Yes');
                
                if strcmp(b2,'Yes')
                    % Delete movie data
                    delete([currDataPath filesep 'centroid movie.mat']);
                end
                
                break
                
            elseif strcmp(b,'Close')
                
                close(f)
                
                break
                
            elseif isempty(b)
                return
            end
        end
        
    else
        error('You have not yet generated the centroid movie data')
    end
    clear v  b b2 anaStatus f mov
    
end


%% Make rotation movie data (do.MakeRotationMovie)

if do.MakeRotationMovie
    
    % If centroid data there and centroid data not yet approved
    if ~isempty(dir([currDataPath filesep 'Rotation.mat'])) && ...
            ~isempty(dir([currDataPath filesep 'Transformation.mat'])) && ...
            isempty(dir([currDataPath filesep 'Rotation movie.mat']))
        
        % Load initial conditions (iC)
        load([currDataPath filesep 'Initial conditions'])
        
        % Load tranfromation strcuture (S)
        load([currDataPath filesep 'Transformation'])
        
        % Load video info (v)                added by CG
        v = defineVidObject(currVidPath);
        
        disp(' ')
        disp(['Making Pred Rotation Movie: ' currVidPath])
        disp(' ')
        
        movieFileName = [cList.fName{i} ' Rotation'];
        
        % Make movie
        aniData(currVidPath,v,currDataPath,movieFileName,imInvert,...
            'Centroid & Rotation',S,0);
        
%         mov.dataPath      = currDataPath;
%         mov.currVidPath   = currVidPath;
%         mov.M             = M;
%         
%         % Save movie data
%         save([currDataPath filesep 'Rotation movie'],'mov','-v7.3')
        
    elseif isempty(dir([currDataPath filesep 'Rotation.mat'])) || ...
            isempty(dir([currDataPath filesep 'Transformation.mat']))
        error('You first need to run manMaster or acqMaster to generate rotation data')
        
    elseif ~isempty(dir([currDataPath filesep 'Rotation movie.mat']))
        
        disp('Rotation movie data already generated')
        
    end
    clear M mov Rotation S iC 
    
end


%% Play tracking movie (do.playRotation)

if do.playRotation 
    
    % If centroid data there and centroid data not yet approved
    if ~isempty(dir([currDataPath filesep 'Rotation movie.mat']))
        
        % Load centroid movie (mov)
        load([currDataPath filesep 'Rotation movie'])
        
        % Make figure window
        f = figure('units','normalized');
        
        % Status
        disp(' ')
        disp(['Approving tracking : ' currVidPath])
        disp(' ');
        
        while true
            
            % Play movie
            movie(f,mov.M,1,20,[0 0 1 1])
            
            % Prompt for approval
            b = questdlg('Which action?','','Replay','Close','Delete','Replay');
            
            if strcmp(b,'Delete')
                
                % Prompt for approval
                b2 = questdlg('Delete the movie data?','','Yes','No','Yes');
                
                if strcmp(b2,'Yes')
                    % Delete movie data
                    delete([currDataPath filesep 'Rotation movie.mat']);
                end

                break
                
            elseif strcmp(b,'Close')
                
                close(f)               
                break
                
            elseif isempty(b)
                return
            end
        end
        close(f)
    end
    clear f M mov Rotation b b2 anaStatus
end




function selectDurations(cList,dataPath,vidPath,do)
% Interactively prompts to select a duration for analysis

% Loop thru sequences
for i = 1:length(cList.vidType)


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