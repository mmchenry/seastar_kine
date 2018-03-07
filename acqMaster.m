function acqMaster
% Acquisition of sea star kinematics

%% Code execution

do.makeMeanImage = 0;

% Copies over portions of movies as image sequences
do.choose_dur = 0;

% Interactively select initial conditions
do.initialConditions = 0;

% Course body tracking
do.courseTrack  = 0;

% Refined body tracking
do.refinedTrack = 0;

% % Whether to track the body centroid
% do.Centroids = 1;
% 
% % Review results of centroid tracking
% do.CentroidReview = 1; 
% 
% % Create movies of the Centroid movies for review
% do.MakeCentroidMovies = 1;
% 
% % Track body rotation
% do.bodyRotation = 1;
% 
% do.MakeRotationMovies = 1;
% 
% do.RotationReview = 1;

% Designate which categories to analyze
do.ana_adult  = 0;
do.ana_juv    = 1;
do.ana_horiz  = 0;
do.ana_vert   = 1;
do.ana_upside = 0;


%% General parameters

% Invert image
imInvert = 1;

% Number of points to define the region of interest
numroipts = 200;

% Number of frames to analyze between saving
numForSave = 10;


%% Manage paths (need to modify for new PC)

if ~isempty(dir(['/Users/mmchenry/Documents/Matlab code']))
    
    % Path to kineBox
    kinePath = '/Users/mmchenry/Documents/Matlab code/kineBox';
    
    % Path to root dir of video (CSULB project, external drive)
    %vidPath = '/Volumes/Video/Sea stars/CSULB test/Raw video';
    %vidPath = '/Users/mmchenry/Documents/Video';
    vidPath = '/Volumes/Video/Sea stars/CSULB/Raw video';
    
    % Location of video frames
    %vidFramePath = '/Volumes/Video/Sea stars/CSULB test/Video frames';
    
    % Path to root of data
    dataPath = '/Users/mmchenry/Documents/Projects/Seastars/CSULB data';
    
    
    
elseif ~isempty(dir(['C:\Program Files\MATLAB\R2016a']))
    
    %vidPath = '\\flow.local\shared\Sea stars';
    vidPath = 'C:\Users\andres\Documents\Sea stars';
    
    dataPath = '\\flow.local\andres\Sea stars\CSULB data';
    
    kinePath = 'C:\Users\andres\Documents\GitPath\kineBox';
    
else
    error('Do not recognize computer')
    
end

% Check paths
if isempty(dir(dataPath))
    error('dataPath missing');
    
elseif isempty(dir(kinePath))
    error('kinePath missing');
    
elseif isempty(dir(vidPath))
    error('vidPath missing');
end


%% Initialize

% Add kinePath
addpath(kinePath)


%% Catalog sequence videos

cList = catVidfiles(vidPath);


%% Find mean image

if do.makeMeanImage
    
    % Loop thru catalog
    for i = 1:length(cList.vidType)
        
        % Current paths
        currDataPath = [dataPath filesep cList.path{i} filesep cList.fName{i}];
        currVidPath  = [vidPath filesep cList.path{i} filesep cList.fName{i} cList.ext{i}];
        
        % Load video info (v)
        v = defineVidObject(currVidPath);
        
        if isempty(dir([currDataPath filesep 'meanImageData.mat']))
            
            % Make mean images
            imMean = motionImage(currVidPath,v,'mean bright','none',imInvert);
            
            % Save mean image data
            save([currDataPath filesep 'meanImageData'],'imMean');           
        end
        
        clear v currDataPath currVidPath imMean
    end     
end


%% Interactive mode: Intialize

if  do.initialConditions
    
    % Loop thru sequences
    for i = 1:length(cList.vidType)
        
        % Current directories
        currDataPath = [dataPath filesep cList.path{i} filesep cList.fName{i}];
        currVidPath  = [vidPath filesep cList.path{i} filesep cList.fName{i} cList.ext{i}];
                
        % Load video info (v)
        v = defineVidObject(currVidPath);
        
        % Load imMean
        load([currDataPath filesep 'meanImageData'])
        
        % Define initial conditions, make data files
        initializeTracking(currVidPath,v,currDataPath,...
                           'single body, with arms',imMean,imInvert)       
    end
    
    clear currDataPath currVidPath v i imMean
end


%% Course body tracking

if do.courseTrack
    
    % Whether to visualize steps in analysis
    visSteps = 0;
    
    % Loop thru sequences
    for i = 1:length(cList.vidType)
        
        % Current directories
        currDataPath = [dataPath filesep cList.path{i} filesep cList.fName{i}];
        currVidPath  = [vidPath filesep cList.path{i} filesep cList.fName{i} cList.ext{i}];

        % Load initial conditions (iC)
        load([currDataPath filesep 'Initial conditions'])
        
        % Load imMean
        load([currDataPath filesep 'meanImageData'])
        
        % Load video info (v)
        v = defineVidObject(currVidPath);
        
        % Define initial conditions, make data files   
         trackBody(currVidPath,v,currDataPath,'Body course.mat',...
                  'course',imMean,iC.frIntervalCourse,iC.x,iC.y,iC.r,visSteps)  
    end
    
    clear currDataPath currVidPath v i imMean iC
end


%% Arm tracking

if 0 % do.refinedTrack
    
    visSteps = 1;
    
    % Loop thru sequences
    for i = 1:length(cList.vidType)
        
        % Current directories
        currDataPath = [dataPath filesep cList.path{i} filesep cList.fName{i}];
        currVidPath  = [vidPath filesep cList.path{i} filesep cList.fName{i} cList.ext{i}];
                
        % Load video info (v)
        v = defineVidObject(currVidPath);
        
        if ~isempty(dir([currDataPath filesep 'Body course.mat']))
            
            % Load course body data (Body)
            load([currDataPath filesep 'Body course.mat'])
            
            % Load imMean
            load([currDataPath filesep 'meanImageData'])
            
            % Load initial conditions (iC)
            load([currDataPath filesep 'Initial conditions'])
            
            iC.frIntervalFine = 5;
            
            % Define initial conditions, make data files
            trackBody(currVidPath,v,currDataPath,'Body arm.mat',...
                'arms',imMean,Body,iC.frIntervalFine,'none',visSteps)           
        end
    end
    
    clear currDataPath currVidPath v i Body iC imMean visSteps
end



%% Visualize arm tracking

if 0 % do.refinedTrack
    
    visSteps = 1;
    
    % Loop thru sequences
    for i = 1:length(cList.vidType)
        
        % Current directories
        currDataPath = [dataPath filesep cList.path{i} filesep cList.fName{i}];
        currVidPath  = [vidPath filesep cList.path{i} filesep cList.fName{i} cList.ext{i}];
                
        % Load video info (v)
        v = defineVidObject(currVidPath);
        
        if ~isempty(dir([currDataPath filesep 'Body arm.mat']))
            
            % Load course body data (Body)
            load([currDataPath filesep 'Body arm.mat'])
            
             % Load initial conditions (iC)
            load([currDataPath filesep 'Initial conditions'])

            % Make movie
            aniData(currVidPath,v,currDataPath,'Body arm.mat',iC.invert,...
                'Body arm',Body,visSteps,iC);
        end
    end
    
    clear currDataPath currVidPath v i visSteps
end



return

%% Refined body tracking

if 0 % do.refinedTrack
    
    % Loop thru sequences
    for i = 1:length(cList.vidType)
        
        % Current directories
        currDataPath = [dataPath filesep cList.path{i} filesep cList.fName{i}];
        currVidPath  = [vidPath filesep cList.path{i} filesep cList.fName{i} cList.ext{i}];
                
        % Load video info (v)
        v = defineVidObject(currVidPath);
        
        if ~isempty(dir([currDataPath filesep 'Body course.mat']))
            
            % Load course body data (Body)
            load([currDataPath filesep 'Body course.mat'])
            
            % Save data file
            save([currDataPath filesep 'Body fine.mat'],'Body')
            
            % Define initial conditions, make data files
            trackBody(currVidPath,v,currDataPath,'Body fine.mat',...
                'refined',[],[],[],1,Body)
        end
    end
    
    clear currDataPath currVidPath v i
end


%% Make body-tracking movie

if 1 %do.MakeRotationMovies
    
    % Loop thru sequences
    for i = 1:length(cList.vidType)
        
        % Current paths
        currDataPath = [dataPath filesep cList.path{i} filesep cList.fName{i}];
        currVidPath  = [vidPath filesep cList.path{i} filesep cList.fName{i} cList.ext{i}];
        
%         % Check approval status
%         if ~isempty(dir([currDataPath filesep 'anaStatus.mat']))
%             load([currDataPath filesep 'anaStatus.mat'])
%         else
%             warning(['You need to approve centroid tracking for' ...
%                 currDataPath])
%         end
%         
%         if ~isfield(anaStatus,'Rotation')
%             anaStatus.Rotation = '';
%         end
        
        % If centroid data there and centroid data not yet approved
        if ~isempty(dir([currDataPath filesep 'Body fine.mat'])) 
            
%             % Load initial conditions (iC)
%             load([currDataPath filesep 'Initial conditions'])
%             
%             % Load Body structure (Body)
            load([currDataPath filesep 'Body fine.mat'])
%             
    % Load video info (v)
        v = defineVidObject(currVidPath);
%             
            disp(' ')
            disp(['Making Pred Rotation Movie: ' currVidPath])
            disp(' ')
            
            % Make movie
            aniData(currVidPath,v,currDataPath,'Body fine.mat',imInvert,'Body',Body);
            
%             mov.dataPath      = currDataPath;
%             mov.currVidPath   = currVidPath;
%             mov.M             = M;
%             
            % Save movie data
            %save([currDataPath filesep 'Rotation movie'],'mov','-v7.3')
            
            clear v currDataPath currVidPath
        end
    end
end


%% Interactive mode: Select initial conditions

% if do.initialConditions
%     
%     % Loop thru sequences
%     for i = 1:length(cList.vidType)
%         
%         currDataPath = [dataPath filesep cList.path{i} filesep cList.fName{i}];
%         currVidPath  = [vidPath filesep cList.path{i} filesep cList.fName{i} cList.ext{i}];
% 
%         
%         % Load video info (v)
%         v = defineVidObject(currVidPath);
%         
%         % If analysis requested and not done already, meets yes_ana
%         % criteria and the clipinfo has already been determined
%         if ~isempty(dir([currDataPath filesep 'clipInfo.mat'])) && ...
%                 yes_ana(cList.age(i),cList.orient(i),do) && ...
%                 isempty(dir([currDataPath filesep 'Initial conditions.mat'])) 
%                           
%             % Load frame intervals ('clipInfo')
%             load([currDataPath filesep 'clipInfo'])
%             
%             if ~isnan(clipInfo.startFrame)
%                 disp(['Initial conditions: ' currVidPath])
%                 disp('')
%                 
%                 % First image
%                 im = getFrame(currVidPath,v,v.UserData.FirstFrame,imInvert,'gray');
%                 
%                 % Initial position
%                 disp(' ')
%                 disp('Select animal to be tracked')
%                 [x,y] = imInteract(im,'points',1);
%                 
%                 % Threshold
%                 disp(' ')
%                 disp('Select threshold')
%                 tVal = imInteract(im,'threshold');
%                 
%                 % Radius
%                 disp(' ')
%                 disp('Select roi radius')
%                 r = imInteract(im,'radius',x,y);
%                 
%                 % Store data
%                 iC.x       = x;
%                 iC.y       = y;
%                 iC.tVal    = tVal;
%                 iC.r       = r;
%                 iC.useMean = 0;
%                 
%                 % Save data
%                 save([currDataPath filesep 'Initial conditions'],'iC')
%                 
%                 clear im useMean im0Mean im0NoMean x y tVal r
%             end
%         end
%     end    
%     disp(' ')  
% end


%% Initialize data structures for analysis

% % Loop thru sequences
% for i = 1:length(cList.vidType)
%     
%     % Load initial conditions (iC)
%     load([currDataPath filesep 'Initial conditions'])
%     
%     % Load frame intervals ('clipInfo')
%     load([currDataPath filesep 'clipInfo'])
%     
%     if isempty(dir([currDataPath filesep 'Centroid.mat']))
%         
%         % Make empty Centroid structure
%         Centroid.frames   = [clipInfo.startFrame:clipInfo.endFrame]';
%         Centroid.x        = nan(size(Centroid.frames));
%         Centroid.y        = nan(size(Centroid.frames));
%         Centroid.y_flip   = nan(size(Centroid.frames));
%         Centroid.xPerim{1}  = nan;
%         Centroid.yPerim{1}  = nan;
%         
%         % Save data
%         save([currDataPath filesep 'Centroid'],'Centroid')
%     end
%     
% %     if isempty(dir([currDataPath filesep 'Rotation.mat']))
% %         
% %         % Make empty Centroid structure
% %         Rotation.frames   = [clipInfo.startFrame:clipInfo.endFrame]';
% %         Rotation.x        = nan(size(Rotation.frames));
% %         Rotation.y        = nan(size(Rotation.frames));
% %         Rotation.y_flip   = nan(size(Rotation.frames));
% %         
% %         % Save data
% %         save([currDataPath filesep 'Rotation'],'Rotation')
% %     end
% end


%% Track centroid coordinates

% if do.Centroids
%     
%     % Loop thru sequences
%     for i = 1:length(cList.vidType)
%         
%         currDataPath = [dataPath filesep cList.path{i} filesep cList.fName{i}];
%         currVidPath  = [vidPath filesep cList.path{i} filesep cList.fName{i} cList.ext{i}];
%         
%         % Check approval status, if file present
%         if ~isempty(dir([currDataPath filesep 'anaStatus.mat']))
%             load([currDataPath filesep 'anaStatus.mat'])
%         else
%             anaStatus.centroid = 'not approved';
%         end
%         
%         if ~strcmp(anaStatus.centroid,'approved')
%         
%         % If analysis requested and not done already . . .
% %         if ~isempty(dir([currDataPath filesep 'Initial conditions.mat'])) && ...
% %                 yes_ana(cList.age(i),cList.orient(i),do) 
%             
% %             disp(['Tracking centroid: ' currVidPath])
% %             disp('')
%             
%             % Load initial conditions (iC)
%             load([currDataPath filesep 'Initial conditions'])
%             
%             % Load frame intervals ('clipInfo')
%             load([currDataPath filesep 'clipInfo'])
%               
%             % Insert number of frames to save
%             iC.numForSave = numForSave;
%             
%             % Frames
%             frames = clipInfo.startFrame:clipInfo.endFrame;
%             
%             % Region of interest for first frame
%             roi0 = giveROI('define','circular',numroipts,iC.r,iC.x,iC.y);
%             
%             % Centroid data file
%             cDataFile = [currDataPath filesep 'Centroid.mat'];
%             
%             % Run tracker code for centroid
%             tracker(currVidPath,v,currDataPath,imInvert,...
%                            'threshold translation',roi0,iC.tVal,frames);
%             
%             % Save data
%             %save([currDataPath filesep 'Centroid'],'Centroid')
%             
%             clear cDataFile
% %         end
%         disp(' ')
%     end    
% end


%% Generate centroid movies for review

if do.MakeCentroidMovies
    
    % Loop thru sequences
    for i = 1:length(cList.vidType)
        
        % Current paths
        currDataPath = [dataPath filesep cList.path{i} filesep cList.fName{i}];
        currVidPath  = [vidPath filesep cList.path{i} filesep cList.fName{i} cList.ext{i}];
        
        % Check status
        if ~isempty(dir([currDataPath filesep 'anaStatus.mat']))
            load([currDataPath filesep 'anaStatus.mat'])
        else
            anaStatus.centroid = '';
        end
        
        % If centroid data there and centroid data not yet approved
        if isempty(dir([currDataPath filesep 'centroid movie.mat'])) && ...
           ~isempty(dir([currDataPath filesep 'Centroid.mat'])) && ...
            isempty(anaStatus.centroid)
            
            % Load initial conditions (iC)
            load([currDataPath filesep 'Initial conditions'])
            
            % Load centroid data (Centroid)
            load([currDataPath filesep 'Centroid.mat'])
            
            % Load frame intervals ('clipInfo')
            load([currDataPath filesep 'clipInfo'])
            
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
            aniData(currVidPath,v,currDataPath,imInvert,'Centroid tracking',S,1);
            
%             mov.dataPath      = currDataPath;
%             mov.currVidPath   = currVidPath;
%             mov.M             = M;
%             
%             % Save movie data
%             save([currDataPath filesep 'centroid movie'],'mov','-v7.3')
        end
    end
end


%% Review recent centroid tracking

if do.CentroidReview
    
    % Loop thru sequences
    for i = 1:length(cList.vidType)
        
        % Current paths
        currDataPath = [dataPath filesep cList.path{i} filesep cList.fName{i}];
        currVidPath  = [vidPath filesep cList.path{i} filesep cList.fName{i} cList.ext{i}];
        
        % Check status
        if ~isempty(dir([currDataPath filesep 'anaStatus.mat']))
            load([currDataPath filesep 'anaStatus.mat'])
        else
            anaStatus.centroid = '';
        end
        
        % If centroid data there and centroid data not yet approved
        if ~isempty(dir([currDataPath filesep 'centroid movie.mat'])) && ...
            isempty(anaStatus.centroid)
            
            % Load centroid movie (mov)
            load([currDataPath filesep 'centroid movie'])

            % Make figure window
            f = figure('units','normalized');
            
            % Status
            disp(' ')
            disp(['Approving tracking : ' currVidPath])
            disp(' ');
            
            while true
                
                movie(f,mov.M,1,20,[0 0 1 1])
                
                b = questdlg('Does the tracking look good?','','Yes','No','Replay','Yes');
                
                if strcmp(b,'Yes')
                    
                    % Update status
                    anaStatus.centroid = 'approved';
                    
                    % Prompt for approval
                     b2 = questdlg('Delete the movie data?','','Yes','No','Yes');
                    
                     if strcmp(b2,'Yes')
                         % Delete movie data
                         delete([currDataPath filesep 'centroid movie.mat']);
                     end
                    
                    % Save status file
                    save([currDataPath filesep 'anaStatus.mat'],'anaStatus')
                    
                    break
                    
                elseif strcmp(b,'No')
                    
                    % Update status
                    anaStatus.centroid = 'not approved';
                    
                    % Delete movie data
                    delete([currDataPath filesep 'centroid movie.mat']);
                    
                    % Save status file
                    save([currDataPath filesep 'anaStatus.mat'],'anaStatus')
                    
                    break
                    
                elseif isempty(b)
                    return
                end
            end   
            close(f)
        end     
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
              
        if ~isempty(dir([currDataPath filesep 'Initial conditions.mat'])) && ...
            strcmp(anaStatus.centroid,'approved')
            
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
    end
end


%% Bundle Centroid and Rotation data

% Loop thru sequences
for i = 1:length(cList.vidType)
    
    % Current paths
    currDataPath = [dataPath filesep cList.path{i} filesep cList.fName{i}];
    currVidPath  = [vidPath filesep cList.path{i} filesep cList.fName{i} cList.ext{i}];
    
    % If centroid data there and centroid data not yet approved
    if isempty(dir([currDataPath filesep 'Transformation.mat']))
        
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
end


%% Review rotation: make movies

if do.MakeRotationMovies
    
    % Loop thru sequences
    for i = 1:length(cList.vidType)
        
        % Current paths
        currDataPath = [dataPath filesep cList.path{i} filesep cList.fName{i}];
        currVidPath  = [vidPath filesep cList.path{i} filesep cList.fName{i} cList.ext{i}];
        
        % Check approval status
        if ~isempty(dir([currDataPath filesep 'anaStatus.mat']))
            load([currDataPath filesep 'anaStatus.mat'])
        else
            warning(['You need to approve centroid tracking for' ...
                currDataPath])
        end
        
        if ~isfield(anaStatus,'Rotation')
            anaStatus.Rotation = '';
        end
        
        % If centroid data there and centroid data not yet approved
        if ~isempty(dir([currDataPath filesep 'Rotation.mat'])) && ...
           ~isempty(dir([currDataPath filesep 'Transformation.mat'])) && ...
            isempty(dir([currDataPath filesep 'Rotation movie.mat'])) && ...
            isempty(anaStatus.Rotation)
            
            % Load initial conditions (iC)
            load([currDataPath filesep 'Initial conditions'])
            
            % Load tranfromation strcuture (S)
            load([currDataPath filesep 'Transformation'])
            
            
            disp(' ')
            disp(['Making Pred Rotation Movie: ' currVidPath])
            disp(' ')
            
            % Make movie
            M = aniData(currVidPath,v,currDataPath,imInvert,'Centroid & Rotation',S,0);
            
            mov.dataPath      = currDataPath;
            mov.currVidPath   = currVidPath;
            mov.M             = M;
            
            % Save movie data
            save([currDataPath filesep 'Rotation movie'],'mov','-v7.3')
        end
    end
end


%% Approve tracking movies

if do.RotationReview
    
    % Loop thru sequences
    for i = 1:length(cList.vidType)
        
        % Current paths
        currDataPath = [dataPath filesep cList.path{i} filesep cList.fName{i}];
        currVidPath  = [vidPath filesep cList.path{i} filesep cList.fName{i} cList.ext{i}];
        
        % Check approval status
        if ~isempty(dir([currDataPath filesep 'anaStatus.mat']))
            load([currDataPath filesep 'anaStatus.mat'])
        else
            error(['You first need to approve centroid tracking for' ...
                currDataPath])
        end
        
        % If centroid data there and centroid data not yet approved
        if ~isempty(dir([currDataPath filesep 'Rotation movie.mat'])) && ...
                ~isfield(anaStatus,'Rotation_pd')
            
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
                b = questdlg('Does the tracking look good?','','Yes','No','Replay','Yes');
                
                if strcmp(b,'Yes')
                    
                    % Update status
                    anaStatus.Rotation_pd = 'approved';
                    
                    % Prompt for approval
                     b2 = questdlg('Delete the movie data?','','Yes','No','Yes');
                    
                     if strcmp(b2,'Yes')
                         % Delete movie data
                         delete([currDataPath filesep 'Rotation movie.mat']);
                     end
                    
                    % Save status file
                    save([currDataPath filesep 'anaStatus.mat'],'anaStatus')
                    
                    break
                    
                elseif strcmp(b,'No')
                    
                    % Update status
                    anaStatus.Rotation_pd = 'not approved';
                    
                    % Delete movie data
                    delete([currDataPath filesep 'Rotation movie.mat']);
                    
                    % Save status file
                    save([currDataPath filesep 'anaStatus.mat'],'anaStatus')
                    
                    break
                    
                elseif isempty(b)
                    return
                end
            end   
            close(f)
        end     
    end   
end



%%
%curr_vidPath = uigetdir(vidPath,'Select video directory')


function yes = yes_ana(age,wallorient,do)

yes = ((do.ana_adult   &&  strcmp(age,'a')) || ...
        (do.ana_juv    &&  strcmp(age,'j'))) && ...
       ((do.ana_horiz  &&  strcmp(wallorient,'h')) || ...
        (do.ana_vert   &&  strcmp(wallorient,'v')) || ...
        (do.ana_upside &&  strcmp(wallorient,'u')));
    


% function [firstFrame,lastFrame] = getFrameNum(firstFrame,lastFrame,im1,im2)
% 
% prompt     = {'Start frame num:','Last frame num:'};
% name= 'Choose clip duration';
% numlines=1;
% defaultanswer={num2str(firstFrame),num2str(lastFrame)};
% 
% answer = inputdlg(prompt,name,numlines,defaultanswer);
% 
% if isempty(answer)
%     return
% end
% 
% firstFrame   = str2num(answer{1});
% lastFrame    = str2num(answer{2});


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
                                m = m + 1;
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