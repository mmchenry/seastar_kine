function anaSingle(fileName)
% Runs analysis of data on indvidual sequence 


%% Code execution

% Analyze foot kinematics from manual tracking
do.anaManTrack = 1;

% Draw trajectory
do.drawTraj = 0;

% Draw all points in trajectory FOR
do.drawAllT = 0;

% Draw rose plots for power stroke
do.rosePlot = 0;

% Gait diagram
do.gaitDiagram = 0;

% Create animation of foot contacts
do.footContactMovie = 0;


%% General parameters

% Color for body of sea star
clr{1} = 0.5.*[1 1 1];


plotclr = plotclrs;

% Paths 
paths = givePaths;


%% Catalog sequence videos

%cList = catVidfiles(paths.vid);

% cList.fName = 'S004_S001_T007';
if nargin<1
%      cList.fName = 'SS001_S001_T013';
% cList.fName = 'S005_S001_T011';
    %cList.fName = 'STUDIO0_S009_S001_T020';
%      cList.fName = ['weights' filesep 'STUDIO0_S009_S001_T009'];
%      cList.fName = ['floats' filesep 'STUDIO0_S009_S001_T020'];
    cList.fName = ['']
else
    
    cList.fName = fileName;
end

cList.ext   = 'MOV';
cList.movtype = 'mov';
cList.path = '';

% % Go into batch mode, if no filename given
% if nargin<1
%     
%     % Get list of sequences
%     cList = dir([dataPath filesep 'S0*']);
%     
%     % Loop thru cList, running anaSingle
%     for i = 1:length(cList)
%         disp(['Batchmode: running ' cList(i).name]); disp(' ')
%         anaSingle(cList(i).name)
%     end
%     
%     % Stop code
%     return
%     
% else  
%     cList.fName = fileName;
% end
% 
% cList.ext   = 'MOV';
% cList.movtype = 'mov';
% cList.path = '';

% Paths for current sequence
currDataPath = [paths.data filesep cList.path filesep cList.fName];
% currDataPath = [dataPath filesep cList.path filesep cList.fName];
%currVidPath  = [vidPath filesep cList.path filesep cList.fName '.' cList.ext];

% % Check video path 
% if ~isfile(currVidPath)
%     error(['Video file does not exist at ' currVidPath]);
% end
% 
% % Load video info (v)
% v = defineVidObject(currVidPath);


%% Load data for sequence

% Load initial conditions (iC)
load([currDataPath filesep 'Initial conditions'])

% Load body kinematics (Body)
load([currDataPath filesep 'Body, post.mat'])

% Load F structure
%load([currDataPath filesep 'post- foot data'])
load([currDataPath filesep 'post- foot refined'])


%% Sequence-specific information

if strcmp(cList.fName,['floats' filesep 'STUDIO0_S009_S001_T020'])
    
    %TODO: These data are made up -- need to input real values
    
    P.indivNum = 1;
    
    % Start and end of bouncing gait (s)
    P.tBounceStart = 17;
    P.tBounceEnd   = 30;
    
     % Start and end of non-bouncing walking (s)
    P.tWalkStart = 0;
    P.tWalkEnd   = 15;
    
    % Body mass (kg)
    P.mass = 10.27e-3;
    
    % Calibration constant (m/pix)
    P.cal = 1;
    
    % Total number of tube feet
    P.numFeet = 20*5;
    
elseif strcmp(cList.fName,['weights' filesep 'STUDIO0_S009_S001_T009'])
    
    %TODO: These data are made up -- need to input real values
    
    P.indivNum = 1;
    
    % Start and end of bouncing gait (s)
    P.tBounceStart = 22;
    P.tBounceEnd   = 70;
    
     % Start and end of non-bouncing walking (s)
    P.tWalkStart = 0;
    P.tWalkEnd   = 20;
    
    % Body mass (kg)
    P.mass = 10.27e-3;
    
    % Calibration constant (m/pix)
    P.cal = 1;
    
    % Total number of tube feet
    P.numFeet = 20*5;

    
elseif strcmp(cList.fName,'S005_S001_T011')
    
    P.indivNum = 6;
    
    % Start and end of bouncing gait (s)
    P.tBounceStart = 22;
    P.tBounceEnd   = 32;
    
     % Start and end of non-bouncing walking (s)
    P.tWalkStart = 0;
    P.tWalkEnd   = 10;
    
    % Body mass (kg)
    P.mass = 10.27e-3;
    
    % Calibration constant (m/pix)
    P.cal = 5.3238e-05;
    
    % Total number of tube feet
    P.numFeet = 20*5;
    
elseif strcmp(cList.fName,'S008_S001_T019')
    
    P.indivNum = 7;
    
    % Start and end of bouncing gait (s)
    P.tBounceStart = 22;
    P.tBounceEnd   = 65;
    
     % Start and end of non-bouncing walking (s)
    P.tWalkStart = 0;
    P.tWalkEnd   = 12;
    
    % Body mass (kg)
    P.mass = 12.98e-3;
    
    % Calibration constant (m/pix)
    P.cal = 5.3238e-05;
    
    % Total number of tube feet
    P.numFeet = 20*5;
    
    
elseif strcmp(cList.fName,'S008_S001_T020')
    
    P.indivNum = 7;
    
    % Start and end of bouncing gait (s)
    P.tBounceStart = 20;
    P.tBounceEnd   = 58;
    
     % Start and end of non-bouncing walking (s)
    P.tWalkStart = 0;
    P.tWalkEnd   = 12;
    
    % Body mass (kg)
    P.mass = 12.98e-3;
    
    % Calibration constant (m/pix)
    P.cal = 5.3238e-05;
    
    % Total number of tube feet
    P.numFeet = 20*5; 
    
elseif strcmp(cList.fName,'S008_S001_T037')
    
    P.indivNum = 1;
    
    % Start and end of bouncing gait (s)
    P.tBounceStart = 40;
    P.tBounceEnd   = 55;
    
     % Start and end of non-bouncing walking (s)
    P.tWalkStart = 0;
    P.tWalkEnd   = 15;
    
    % Body mass (kg)
    P.mass = 15.32e-3;
    
    % Calibration constant (m/pix)
    P.cal = 5.3238e-05;
    
    % Total number of tube feet
    P.numFeet = 20*5; 
    
    
elseif strcmp(cList.fName,'S007_S001_T051')
    
    P.indivNum = 5;
    
    % Start and end of bouncing gait (s)
    P.tBounceStart = 12;
    P.tBounceEnd   = 40;
    
     % Start and end of non-bouncing walking (s)
    P.tWalkStart = 0;
    P.tWalkEnd   = 7;
    
    % Body mass (kg)
    P.mass = 16.32e-3;
    
    % Calibration constant (m/pix)
    P.cal = 5.3238e-05;
    
    % Total number of tube feet
    P.numFeet = 20*5; 
    
    
elseif strcmp(cList.fName,'S006_S001_T015')
    
    P.indivNum = 3;
    
    % Start and end of bouncing gait (s)
    P.tBounceStart = 15;
    P.tBounceEnd   = 25;
    
     % Start and end of non-bouncing walking (s)
    P.tWalkStart = nan;
    P.tWalkEnd   = nan;
    
    % Body mass (kg)
    P.mass = 13.8e-3;
    
    % Calibration constant (m/pix)
    P.cal = 5.3238e-05;
    
    % Total number of tube feet
    P.numFeet = 20*5; 

end



%% Analyze manually-tracked foot kinematics

if do.anaManTrack 

ttt = 3;

end %do.anaManTrack 

%% Visualize trajectory

if do.drawTraj

    % Times for drawing body
    tDraw = [Body.t(10) Body.t(end-10)];
    
    figure
    
    for i = 1:length(tDraw)
        
        % Time difference btwen current and all times
        tDiff = abs(tDraw(i)-Body.t);
        
        % Index is closest time to desired
        idx = find(tDiff==min(tDiff),1,'first');
        
        % Points for whole body
        starPtsT = makeStar(Body.xArmT(idx,:)', Body.yArmT(idx,:)', Body.xT(idx),Body.yT(idx));
        
        % Draw body at interval
        h = fill(starPtsT(:,1),starPtsT(:,2),clr{1},'FaceAlpha',0.3,...
            'EdgeColor','none');
        
        hold on
        
        % Overlay the tube feet
        for j = 1:length(F)
            
            % Index for current frame   
            idxF = F(j).frames==Body.frames(idx);
                
            xVal = F(j).xT(idxF);
            yVal = F(j).yT(idxF);
            
            scatter(xVal,yVal,'SizeData',100,'MarkerEdgeColor','none', ...
                     'MarkerFaceColor',0.8.*[1 1 1]);
            clear idxF xVal yVal
        end

        clear starPntsT idx tDiff
    end
    
    % Plot center points
    h = scatter(Body.xT,Body.yT,'SizeData',10,'MarkerEdgeColor','none', ...
        'MarkerFaceColor',0.5.*[1 1 1]);
    
    line([min(Body.xT) max(Body.xT)],[0 0],'Color','k')
    axis equal
    set(gca,'YColor','none','XColor','none')
    
end


%% Animate foot contact

if do.footContactMovie

    visSteps = 0;
    
    % Write vidoe file to disk
    makeFile = 1;
    
    % Filename
    filename = 'foot traj';
    
     % Figure color and positon
    fColor = 1.*[1 1 1];
    fPos   = [ 1 530  2002  995];
    
    % Color for body of sea star
    clr{1} = 0.1.*[1 1 1];

    % Initialize file
    if makeFile
        
        % Set up output video file
        vOut = VideoWriter([currDataPath filesep filename '.mp4'],'MPEG-4');
        vOut.Quality = 50;
        
        open(vOut)
        
        set(0,'DefaultFigureWindowStyle','normal')
        
        % Make figure
        f = figure('Position',fPos,'Color',fColor);
        
        if ~visSteps
            set(f,'Visible','off')
        end
        
    else
        f = figure;
    end
    
    % Pull default colormap
    cmap = colormap(f,'cool');
    cval = linspace(0,1,size(cmap,1));
    
    % Frame rate
    frRate = 1./mean(diff(Body.t));

    % Index for starting frame
    iStart = Body.frames==F(1).frames(1);
    
    % Start time
    tStart = Body.t(iStart);

    % Index of body data to include
    iBod = Body.frames>=F(1).frames(1);
    
    % Points for whole body
    starPtsT = makeStar(Body.xArmTL(iStart,:)', Body.yArmTL(iStart,:)',0,0);
    
    % Draw body
    h = fill(starPtsT(:,1),starPtsT(:,2),clr{1},'FaceAlpha',0.3,...
             'EdgeColor','none');   
    hold on
    axis equal
    
    % Frames
    frames = Body.frames(iBod);
    
    % Time vector
    time = frames ./ frRate - tStart;
    
    % Start index
    j = 1;
    
    % Loop thru feet, saving properties
    for i = 1:length(F)
        
        % If there are coodrinates (excluding feet that start at first frame)
        if max(~isnan(F(i).armNum)) %&& F(i).frames(1)>F(1).frames(1)
            
            % Mean position along arm
%             armPos(j,1) = mean(F(i).xA);
            
            % Arm number
%             armNum(j,1)     = F(i).armNum(1);
            
            % Index of current foot
            ftIdx(j,1) = i;
            
            % Start frame
            frStart(j,1) = F(i).frames(1);
            
            % End frame
            frEnd(j,1) = F(i).frames(end);
            
            % Angle of angle data
            minAng(j,1) = min(F(i).ftAng);
            maxAng(j,1) = max(F(i).ftAng);
            
            % Advance index
            j = j + 1;
        end
    end
    
    % Calc normalized foot angle
    for i = 1:length(frStart)
        angNorm{i} = (F(ftIdx(i)).ftAng-minAng(i))./(maxAng(i)-minAng(i));
    end

    % Loop trhu time
    for i = 1:length(frames)
        
        % Current frame
        cFrame = frames(i);
        
        % Current time
        cTime = time(i);
        
         % Index for feet that include this time
        cIdx = find(cFrame>=frStart & cFrame<=frEnd);
        
        % Step trhu feet
        for j = 1:length(cIdx)
            
            % Number of current foot from F structure
            cFt = ftIdx(cIdx(j));
            
            % Index for duration up to present time
            iDur = F(cFt).frames<=cFrame;
            
            % Time values over duration
            tDur = F(cFt).frames(iDur)./frRate - tStart;
            
            % Current coordinates
            cX = F(cFt).xTL(iDur);
            cY = F(cFt).yTL(iDur);
            
            % Interpolate for angular values at ends of patch
            angVals = interp1(F(cFt).frames./frRate-tStart,angNorm{cIdx(j)},tDur);
            
            % Interpolate color
            col(:,1) = interp1(cval,cmap(:,1),angVals);
            col(:,2) = interp1(cval,cmap(:,2),angVals);
            col(:,3) = interp1(cval,cmap(:,3),angVals);

            % Plot
            h(j) = scatter(cX,cY,[],col,'filled','SizeData',100, ...
                    'MarkerEdgeColor','none');
            
            clear col cX cY angVals tDur iDur cFt
        end
        
        set(gca,'visible','off')
        hTxt = text(-500,-500,['Time = ' num2str(cTime,'%15.2f') ' s']);
        set(hTxt,'FontSize',24,'Color',0.8.*[1 1 1])
        
            disp(['Making movie: done ' num2str(i) ' of ' ...
                num2str(length(time))]);
            
        if visSteps
            pause(0.01)
        end
        
        if makeFile
            imFrame = getframe(f);
            writeVideo(vOut,imFrame);
        end
        
        if exist('h')
            delete(h)
        end
        
        if exist('hTxt')
            delete(hTxt)
        end
        clear h imFrame
        
    end
    
    
    
    if makeFile
        close(vOut)

        set(0,'DefaultFigureWindowStyle','docked')
    end


end


%% Plot 

if do.drawAllT
    
    figure;
    
    % Number of arms
    numArms = max(F(1).armNum);
    
    % Index for starting frame in body data
    iStart = find(Body.frames==F(1).frames(1),1,'first');
    
    % Points for whole body
    starPtsT = makeStar(Body.xArmT(iStart,:)', Body.yArmT(iStart,:)', ...
        Body.xT(iStart),Body.yT(iStart));
    
    % Draw body at interval
    h = fill(starPtsT(:,1),starPtsT(:,2),clr{1},'FaceAlpha',0.3,...
        'EdgeColor','none');
    hold on
    
    % Start index
    k = 1;
    
    % Loop thru blobs
    for i = 1:length(F)
        
        % Step thru frames
        for j = 1:length(F(i).frames)
            
            % Index in body
            iBod = find(Body.frames==F(i).frames(j),1,'first');
            
            % Coordinates
            xVal = F(i).xT(j) - Body.xT(iBod);
            yVal = F(i).yT(j) - Body.yT(iBod);
            
            % Plot points
            h = scatter(xVal,yVal,'MarkerEdgeColor','none',...
                'MarkerFaceColor',plotclr(k,:),'MarkerFaceAlpha',0.5);
            
        end
        
        if k==size(plotclr,1)
            k = 1;
        else
            k = k + 1;
        end
        
        % Status update
        disp(['Plotting . . . (' num2str(i) ' of ' num2str(length(F)) ')'])
    end 
end


%% Rose plots of power stroke direction

if do.rosePlot
    
% Time interval to evaluate  (default is [0 inf] 
tInterval = [0 inf];    
    
% Number of arms
numArms = max(F(1).armNum);

minNum = 50;

% Frame rate
frRate = 1./mean(diff(Body.t));

% Index for starting frame in body data
iStart = find(Body.frames==F(1).frames(1),1,'first');

% Start time
tStart = Body.t(iStart);

figure;

% Points for whole body
starPtsT = makeStar(Body.xArmT(iStart,:)', Body.yArmT(iStart,:)', ...
                    Body.xT(iStart),Body.yT(iStart));

 % Draw body at interval
    h = fill(starPtsT(:,1),starPtsT(:,2),clr{1},'FaceAlpha',0.3,...
        'EdgeColor','none');
  axis equal
  
% f1 = figure;
f2 = figure;

% Loop thru arms
for i = 1:numArms
    
    ftAng = [];
    
    % Loop thru feet
    for j = 1:length(F)
       
        if ~max(isnan(F(j).armNum)) && max(F(j).armNum==i) && ...
                length(F(j).frames)>minNum
            
            % Time vector for current foot
            currTime = F(j).frames./frRate - tStart;
            
            % Index for current period
            cIdx = currTime>=tInterval(1) & currTime<=tInterval(2);
            
            % Index in F data
            iFStart = find(cIdx,1,'first');
            iFEnd   = find(cIdx,1,'last');
            
            % Index in Body
            iStart   = find(Body.frames==F(j).frames(iFStart),1,'first');
            iEnd     = find(Body.frames==F(j).frames(iFEnd),1,'first');
            
            if ~isempty(iFStart) && ~isempty(iFEnd) && ~isempty(iStart) && ...
                    ~isempty(iEnd)
                % Start point
                pStart = [F(j).xT(iFStart)-Body.xT(iStart)  ...
                    F(j).yT(iFStart)-Body.yT(iStart)];
                
                % End point
                pEnd = [F(j).xT(iFEnd)-Body.xT(iEnd)  ...
                    F(j).yT(iFEnd)-Body.yT(iEnd)];
                
                % Add angular value
                ftAng = [ftAng; atan2(pEnd(2)-pStart(2),pEnd(1)-pStart(1))];
                
            end
            % Plot individual trajectories
%             figure(f1)
%             subplot(2,3,i)
%             plot([pStart(1) pEnd(1)],[pStart(2) pEnd(2)],'-k',...
%                 pStart(1),pStart(2),'ok')
%             hold on
            
        end
        
    end

    figure(f2)
    subplot(2,3,i)
    %h = polarhistogram(ftAng,60);
    [a,rTick] = circ_plot_mjm(ftAng);
    %title(['Arm ' num2str(i)])
  
end

end


%% Gait diagram

if do.gaitDiagram
    
    % Number of pairs of tube feet along arm
    numPairs = 20;
    
    %gaitPlot(Body,F,'????')
    
    gaitAngPlot(Body,F,cList.fName,numPairs)
    
end


%% Analyze metrics of walking and bouncing

% Min number of feet to include in polar order parameter
minFeet = 3;

% Frame rate
frRate = 1./mean(diff(Body.t));

% Extract stats on each contact phase
fS = contactStats(F,frRate);

% Polar-order parameter
[OP,nFt] = polarParam(Body,F,fS,frRate,minFeet);

% Indicies for walking and bouncing from Body data
iBodWalk   = Body.t>=P.tWalkStart & Body.t<=P.tWalkEnd;
iBodBounce = Body.t>=P.tBounceStart & Body.t<=P.tBounceEnd;

% Frame rate
frRate = 1./mean(diff(Body.t));

% Containers
wDur = [];
bDur = [];

% Loop thru feet
for i = 1:length(F)
    
    tCurr = F(i).frames./frRate;
    
    % Save time,  if duration is within walk period
    if (min(tCurr)>= P.tWalkStart) && (max(tCurr)<= P.tWalkEnd)
        
        wDur = [wDur; range(tCurr)];
        
    % Save time,  if duration is within bounce period
    elseif (min(tCurr)>= P.tBounceStart) && (max(tCurr)<= P.tBounceEnd)
        
        bDur = [bDur; range(tCurr)];
        
    end
end

% Body diameter (pix)
bDiam = max([nanmax(nanmax(Body.xArmTL,[],2)-nanmin(Body.xArmTL,[],2)), ...
             nanmax(nanmax(Body.yArmTL,[],2)-nanmin(Body.yArmTL,[],2))]);

% Calculate speed (blody lengths)
spd = hypot(diff(smooth(Body.xCntr)),diff(smooth(Body.yCntr)))./diff(Body.t);
spd = [spd(1); spd]./bDiam;

% Speed measurement (m/s)
% spd = P.cal.*hypot(diff(smooth(Body.xCntr)),...
%                    diff(smooth(Body.yCntr)))...
%                 ./diff(Body.t);
% spd = [spd(1); spd];

% Displacement 


% Store power stroke duration
D.meanPwrDur.w = nanmean(wDur);
D.meanPwrDur.b = nanmean(bDur);

% Store speed
D.meanSpeed.w  = nanmean(spd(iBodWalk));
D.meanSpeed.b  = nanmean(spd(iBodBounce));
D.maxSpeed.w   = nanmax(spd(iBodWalk));
D.maxSpeed.b   = nanmax(spd(iBodBounce));

% Store polar parameter
D.meanOP.w     = nanmean(OP(iBodWalk));
D.meanOP.b     = nanmean(OP(iBodBounce));

% Store proportion of feet in contact at each time
D.nFt.w       = mean(nFt(iBodWalk))   / P.numFeet;
D.nFt.b       = mean(nFt(iBodBounce)) / P.numFeet;

clear iBodWalk iBodBounce OP

% Indicies in foot data for walking and bouncing
iWalk   = fS.tEnd>=P.tWalkStart   & fS.tEnd<=P.tWalkEnd;
iBounce = fS.tEnd>=P.tBounceStart & fS.tEnd<=P.tBounceEnd;

% Store durations
D.meanDur.w = nanmean((fS.tEnd(iWalk)   - fS.tStart(iWalk)));
D.meanDur.b = nanmean((fS.tEnd(iBounce) - fS.tStart(iBounce)));

% Save data
save([currDataPath filesep 'summary stats'],'D');






function fS = contactStats(F,frRate)
% Extract stats on each contact phase

% Starting time of the video
%tMin = F(1).frames(1)./frRate;

j = 1;

% Loop thru feet, saving properties
for i = 1:length(F)
    
     % If there are coodrinates (excluding feet that start at first frame)
     if max(~isnan(F(i).armNum)) %&& F(i).frames(1)>F(1).frames(1)

         % Mean position along arm
         fS.armPos(j,1) = mean(F(i).xA);
        
         % Arm number
         fS.armNum(j,1)     = F(i).armNum(1);
         
         % Index of current foot
         fS.ftIdx(j,1) = i;
         
         % Start frame
         fS.frStart(j,1) = F(i).frames(1);
         
         % End frame
         fS.frEnd(j,1) = F(i).frames(end);   
         
         % Range of angular diaplacement
         fS.minAng(j,1) = min(F(i).ftAng);
         fS.maxAng(j,1) = max(F(i).ftAng);
         
         % Body displacement
         fS.displ(j,1) = range(F(i).xT);
         
         % Advance index
         j = j + 1;     
     end 
end

% Start and end times
fS.tStart  = fS.frStart./frRate;
fS.tEnd    = fS.frEnd./frRate;

% Calc normalized foot angle
for i = 1:length(fS.armNum)
    fS.angNorm{i} = (F(fS.ftIdx(i)).ftAng-fS.minAng(i))./(fS.maxAng(i)-fS.minAng(i));
end


function [OP,nFt] = polarParam(Body,F,fS,frRate,minFeet)


% Width of body along x-axis
normLen = 1*range(Body.xArmTL(1,:));

% Loop thru time values
for i = 1:length(Body.t)

    % Current frame
    cFrame = Body.frames(i);
    
    % Current time
    cTime = Body.t(i);
    
    % Indicies of foot data include current frame
    iFeet = cFrame>=fS.frStart & cFrame<=fS.frEnd;
    
    % Indicies from F structure
    iF = fS.ftIdx(iFeet);  
    
    % initialize order parameter to zero 
    z = 0;  
    
    % Number of feet in contact
    nFt(i,1) = length(iF);
    
    % If enough feet in contact
    if ~isempty(iF) && length(iF)>minFeet
        
        % Loop thru each matching foot
        for j = 1:length(iF)
            
            % Time values for current foot
            tVal = F(iF(j)).frames./frRate;
            
%             cFtAng(j,1) = 2*pi*interp1(tVal,cAngVal,cTime);
            
            % Current normalized power stroke
            cAngVal = 2.*pi.*(F(iF(j)).ftAng-min(F(iF(j)).ftAng))./...
                      range(F(iF(j)).ftAng);
            
            % Find current angle for this foot
            cFtAng(j,1) = interp1(tVal,cAngVal,cTime);
            
            % Add to order parameter
            z = z + exp(1i*cFtAng(j));
            
            cDispl(j,1) = 1;
            
        end    
        
        % Normalize x by number of feet
        z = z/length(iF);
        
        % compute absolute value which is the polar order parameter
        OP(i,1) = abs(z); 
        
    else
        OP(i,1) = nan;
    end
    
    clear iF iFeet cIdx cFrame cDispl
end





function cmap = plotclrs

% Colormap of line colors
cmap(1,:) = [0         0.4470    0.7410];
cmap(2,:) = [0.8500    0.3250    0.0980];
cmap(3,:) = [0.9290    0.6940    0.1250];
cmap(4,:) = [0.4940    0.1840    0.5560];
cmap(5,:) = [0.4660    0.6740    0.1880];
cmap(6,:) = [0.3010    0.7450    0.9330];
cmap(7,:) = [0.6350    0.0780    0.1840];


function [starPts,ALen] = makeStar(xA,yA,xC,yC)
% Create periperal shape of seastar body

nArms = length(xA);

% xA = [xA;xA(1)];
% yA = [yA;yA(1)];

nPts = 500;

% Radial position of arms, sorted
[Atheta,idx] = sort(wrapTo2Pi(atan2(yA-yC,xA-xC)));

% Arm lengths
ALen = mean(hypot(xA(idx)-xC,yA(idx)-yC));

% Diameter of cental disc
D = mean(ALen)*0.4;

% Tack on starting arm
Atheta = unwrap([Atheta; Atheta(1)]);
%ALen   = [ALen; ALen(1)];

tVal0 = linspace(0,pi,round(nPts/nArms))';
tVal1 = linspace(pi,2*pi,round(nPts/nArms))';

xP = []; yP = [];

% Loop thru arms
for i = 1:nArms
    
    % Radial positions of current and next arms
    %Atheta0 = atan2(yA(i)-yC,xA(i)-xC);
    %Atheta1 = atan2(yA(i+1)-yC,xA(i+1)-xC);
    
    % Arm lengths of current and next arms
    %ALen0 = hypot(xA(i)-xC,yA(i)-yC);
    %ALen1 = hypot(xA(i+1)-xC,yA(i+1)-yC);
    
    % Radial position values wrapTo2Pi
    
    tmp1 = min([Atheta(i),mean([Atheta(i) Atheta(i+1)])]);
    tmp2 = max([Atheta(i),mean([Atheta(i) Atheta(i+1)])]);    
    theta0 = linspace(tmp1,tmp2,round(nPts/nArms))';
    
    
    tmp1 = min([Atheta(i+1),mean([Atheta(i) Atheta(i+1)])]);
    tmp2 = max([Atheta(i+1),mean([Atheta(i) Atheta(i+1)])]); 
    theta1 = linspace(tmp1,tmp2,round(nPts/nArms))';
    
    rTmp0 = tVal0.^0.7;
    rTmp1 = tVal1.^0.7;

    rVal0 = (ALen-D).*(cos(rTmp0./max(rTmp0)*pi)+1)/2+D;
    rVal1 = rVal0(end:-1:1);
    
    
%     rVal0 = (ALen-D).*(cos(tVal0)+1)/2+D;
%     rVal1 = (ALen-D).*(cos(tVal1)+1)/2+D;
    
    [x0,y0] = pol2cart(theta0,rVal0);
    [x1,y1] = pol2cart(theta1,rVal1);
    
    xP = [xP; x0; x1];
    yP = [yP; y0; y1];
    
    if 0
        
       subplot(2,2,1)
       plot(tVal0,rVal0,'-',tVal1,rVal1,'-')
       %plot(theta0,rVal0,'-',theta1,rVal1,'-')
       hold on
       
       subplot(2,2,3)
       plot(tVal0,theta0,'-',tVal1,theta1,'-')
       hold on
       
       subplot(2,2,[2 4])
       plot(x0,y0,'-',x1,y1,'-')
       axis equal
       hold on
       
       ttt=3;
    end
end

starPts = [xP+xC yP+yC];


function gaitPlot(Body,F,ttext)

% Number of arms
numArms = size(Body.xArmT,2);

% Height of all bars for one arm
hBar = 1;

% Color of the seastar body
sClr = 0.8.*[1 1 1];

% Make figure window
f = figure('Units','normalized');

% Transparency of bars
fAlpha = 0.5;

% Pull default colormap
cmap = colormap(f,'parula');
cval = linspace(0,1.2,size(cmap,1));

% Initialize index
j = 1;

% Frame rate
frRate = 1./mean(diff(Body.t));

% Starting time of the video
%tMin = F(1).frames(1)./frRate;

% Loop thru feet, saving properties
for i = 1:length(F)
    
     % If there are coodrinates . . .
     if max(~isnan(F(i).armNum)) 

         % Mean position along arm
         armPos(j,1) = mean(F(i).xA);
        
         % Arm number
         armNum(j,1)     = F(i).armNum(1);
         
         % Start frame
         frStart(j,1) = F(i).frames(1);
         
         % End frame
         frEnd(j,1) = F(i).frames(end);        
         
         % Advance index
         j = j + 1;     
     end 
end

% Loop trhu feet to assign bar color
for i = 1:length(armNum)
    
    % Color for bar
    clrFt(i,1) = interp1(cval,cmap(:,1),armPos(i)./max(armPos));
    clrFt(i,2) = interp1(cval,cmap(:,2),armPos(i)./max(armPos));
    clrFt(i,3) = interp1(cval,cmap(:,3),armPos(i)./max(armPos));
end

clear cmap cval

% Width of body along x-axis
normLen = range(Body.xArmT(1,:));

% Index of body values to include
iBodStart = find(Body.frames==min(frStart),1,'first');
iBodEnd   = find(Body.frames==max(frEnd),1,'first');

% Full index of Body values 
iBod = iBodStart:iBodEnd;

% Start time
tStart = Body.t(iBodStart);

% Mean x and y arm coordinates
meanArmx = 0.8*hBar .* Body.xArmTL(iBodStart,:)'./normLen;
meanArmy = 0.8*hBar .* Body.yArmTL(iBodStart,:)'./normLen;

% Angular position of arms
armAng = atan2(Body.yArmTL(iBodStart,:),Body.xArmTL(iBodStart,:));

% Indicies for arm numbers, sorted by angular position
[tmp,iArm] = sort(abs(armAng));

% Points for whole body
[starPts,armLen] = makeStar(meanArmx,meanArmy,0,0);

% Make legend for tube foot colors
if 0
    
    figure
    subplot(4,5,1)
    for i = 1:maxFt
        sSize = 0.25;
        
        yV = [i-sSize i+sSize i+sSize i-sSize i-sSize];
        xV = [-sSize -sSize sSize sSize -sSize];
        
        h = fill(xV, yV, cmap2(i,:));
        set(h,'EdgeColor','none');
        
        h = text(1.1*sSize,i,num2str(i));
        set(h,'FontSize',16,'HorizontalAlignment','left','VerticalAlignment','middle')
        
        hold on
    end
    set(gca,'XColor','none','YColor','none');
    %axis equal
    hold off
    title('Tube foot number')
    
    clear sSize h xV yV  
end

% Body diameter (pix)
bDiam = max([nanmax(nanmax(Body.xArmTL,[],2)-nanmin(Body.xArmTL,[],2)), ...
             nanmax(nanmax(Body.yArmTL,[],2)-nanmin(Body.yArmTL,[],2))]);

% Calculate speed
spd = hypot(diff(smooth(Body.xCntr)),diff(smooth(Body.yCntr)))./diff(Body.t);
spd = [spd(1); spd]./bDiam;

% Calculate heading
for i = 1:length(Body.xArmT)
    theta(i,1) = 180/pi*atan2(Body.yArmT(i,1)-Body.yT(i),Body.xArmT(i,1)-Body.xT(i));
end

% Plot speed and heading
figure(f);
subplot(4,5,2:5)
h(1) = plot(Body.t(iBod)-tStart,spd(iBod));ax1=gca;
% [ax1,h(1),h(2)] = plotyy(Body.t(iBod)-tStart,spd(iBod),...
% Body.t(iBod)-tStart,theta(iBod));
ylabel(ax1(1),'Speed (diameters/s)')
%ylabel(ax1(2),'Heading (deg)')

%xlabel('Time s)');
set(h(1),'Color','k','LineWidth',1.5)
set(ax1(1),'YColor','k')
% set(h(2),'Color',0.7.*[1 1 1],'LineWidth',1.5)
set(ax1,'TickDir','out')
%set(ax1(2),'YColor',0.5.*[1 1 1])
title(ttext)
grid on
% Get x-limits
xL = xlim;

% Create axes
h1 = subplot(4,5,6:5:16);
h2 = subplot(4,5,[7:10 12:15 17:20]);

% Order the arms according to their angle wrt the trajectory


% Loop thru arms and plot seastar bodies
for i = 1:numArms
    axes(h1)
    h = fill(starPts(:,1),starPts(:,2)+i-0.5,sClr,'EdgeColor','none');
    axis equal
    hold on
    
    h = line([0 meanArmx(iArm(i))],[i meanArmy(iArm(i))+i]-0.5,...
        'Color',0.4.*[1 1 1],'LineWidth',2);

    drawnow
end
ylim([0 5])


axes(h2)

bClr = 0.95.*[1 1 1];

arms = 2:2:6;

% Loop thru every-other arm to create light bars in background
for i = 1:length(arms)
    
    % Current arm number
    n = arms(i);
    
    % Coords for bar
    xVals = [xL(1) xL(2) xL(2) xL(1) xL(1)];
    yVals = [n n n+hBar n+hBar n];

    % Make bar
    h = fill(xVals,yVals,bClr,'EdgeColor','none');
    hold on
end


% Loop thru each footprint
for i = 1:length(armPos)

    % Current arm 
    n = find(iArm==armNum(i),1,'first');
    
    % Current foot number
    nFt = armPos(i);

    % Micro bar height
    mhBar = 0.05;
    
    % Start and enf times for the power stroke
    tStartFt = frStart(i)./frRate - tStart;
    tEndFt   = frEnd(i)./frRate - tStart;
    
    % Coords for bar
    xVals = [tStartFt tEndFt tEndFt tStartFt tStartFt];
    yVals = [n n n+mhBar n+mhBar n];
    
    % Offset y by foot position along arm
    yVals = yVals + (armPos(i)-min(armPos))./range(armPos);
    
    % Plot box
    axes(h2)
    h = fill(xVals,yVals,clrFt(i,:));
    set(h,'FaceAlpha',fAlpha,'EdgeColor','none')
    hold on
end

grid on
ylim([1 6])

% Position of top graph
pos3 = get(ax1(1),'Position');

% Position of gait graph
pos2 = get(h2,'Position');

% Position of body graphs
pos1 = get(h1,'Position');

% Adjust position of gait graph
% set(h2,'Position',[pos2(1)*0.8 pos2(2)*1.3 pos2(3) pos1(4)*.88],...
%     'YColor','none','TickDir','out');
 set(h2,'YColor','none','TickDir','out');

% Adjust position of top graph
set(h1,'XColor','none','YColor','none','Position',[pos1(1)*1.4 pos1(2:4)]);

set(ax1,'XColor','none')

% set(ax1,'Position',[pos2(1)*0.8 pos3(2) pos2(3) pos3(4)])
xlabel('Time (s)')


function gaitAngPlot(Body,F,ttext,numPairs)

% Number of arms
numArms = size(Body.xArmT,2);

% Height of all bars for one arm
hBar = 1;

% Color of the seastar body
sClr = 0.8.*[1 1 1];

% Make figure window
f = figure('Units','normalized');

% Transparency of bars
fAlpha = 0.8;

% Pull default colormap
cmap = colormap(f,'cool');
cval = linspace(0,1,size(cmap,1));

% Initialize index
j = 1;

% Frame rate
frRate = 1./mean(diff(Body.t));

% Starting time of the video
%tMin = F(1).frames(1)./frRate;

% Loop thru feet, saving properties
for i = 1:length(F)
    
     % If there are coodrinates (excluding feet that start at first frame)
     if max(~isnan(F(i).armNum)) %&& F(i).frames(1)>F(1).frames(1)

         % Mean position along arm
         armPos(j,1) = mean(F(i).xA);
        
         % Arm number
         armNum(j,1)     = F(i).armNum(1);
         
         % Index of current foot
         ftIdx(j,1) = i;
         
         % Start frame
         frStart(j,1) = F(i).frames(1);
         
         % End frame
         frEnd(j,1) = F(i).frames(end);   
         
         minAng(j,1) = min(F(i).ftAng);
         maxAng(j,1) = max(F(i).ftAng);
         
         % Advance index
         j = j + 1;     
     end 
end

% Calc normalized foot angle
for i = 1:length(armNum)
    angNorm{i} = (F(ftIdx(i)).ftAng-minAng(i))./(maxAng(i)-minAng(i));
end


% Width of body along x-axis
iS = find(~isnan(Body.xArmTL(:,1)),1);
normLen = 1*range(Body.xArmTL(iS,:));

% Index of body values to include
iBodStart = find(Body.frames==min(frStart),1,'first');
iBodEnd   = find(Body.frames==max(frEnd),1,'first');

% Full index of Body values 
iBod = iBodStart:iBodEnd;

% Start time
tStart = Body.t(iBodStart);

% Mean x and y arm coordinates
meanArmx = 0.8*hBar .* Body.xArmTL(iBodStart,:)'./normLen;
meanArmy = 0.8*hBar .* Body.yArmTL(iBodStart,:)'./normLen;

% Min number of feet to include in polar order parameter
minFeet = 3;

% Loop thru time values
for i = 1:length(iBod)
    
    % Current index in Body data
    cIdx = iBod(i);
    
    % Current frame
    cFrame = Body.frames(cIdx);
    
    % Current time
    cTime = cFrame./frRate - tStart;
    
    % Indicies of foot data include current frame
    iFeet = cFrame>=frStart & cFrame<=frEnd;
    
    % Indicies from F structure
    iF = ftIdx(iFeet);  
    
    % initialize order parameter to zero 
    z = 0;  
    
    if ~isempty(iF) && length(iF)>minFeet
        % Loop thru each matching foot
        for j = 1:length(iF)
            
            % Time values for current foot
            tVal = F(iF(j)).frames./frRate - tStart;
            
%             cFtAng(j,1) = 2*pi*interp1(tVal,cAngVal,cTime);
            
            % Current normalized power stroke
            cAngVal = 2.*pi.*(F(iF(j)).ftAng-min(F(iF(j)).ftAng))./...
                      range(F(iF(j)).ftAng);
            
            % Find current angle for this foot
            cFtAng(j,1) = interp1(tVal,cAngVal,cTime);
            
            % Add to order parameter
            z = z + exp(1i*cFtAng(j));
        end    
        
        % Normalize x by number of feet
        z = z/length(iF);
        
        % compute absolute value which is the polar order parameter
        OP(i,1) = abs(z); 
        
    else
        OP(i,1) = nan;
    end
    
    clear iF iFeet cIdx cFrame
end


% Angular position of arms
armAng = atan2(Body.yArmTL(iBodStart,:),Body.xArmTL(iBodStart,:));

% Body diameter (pix)
bDiam = max([nanmax(nanmax(Body.xArmTL,[],2)-nanmin(Body.xArmTL,[],2)), ...
             nanmax(nanmax(Body.yArmTL,[],2)-nanmin(Body.yArmTL,[],2))]);

% Indicies for arm numbers, sorted by angular position
[tmp,iArm] = sort(abs(armAng));

% Points for whole body
[starPts,armLen] = makeStar(meanArmx,meanArmy,0,0);

% Calculate speed
spd = hypot(diff(smooth(Body.xCntr)),diff(smooth(Body.yCntr)))./diff(Body.t);
spd = [spd(1); spd] ./ bDiam;

% Calculate heading
for i = 1:length(Body.xArmT)
    theta(i,1) = 180/pi*atan2(Body.yArmT(i,1)-Body.yT(i),Body.xArmT(i,1)-Body.xT(i));
end

% Plot speed and heading
figure(f);
subplot(4,5,2:5)
h(1) = plot(Body.t(iBod)-tStart,spd(iBod));ax1=gca;
% [ax1,h(1),h(2)] = plotyy(Body.t(iBod)-tStart,spd(iBod),...
% Body.t(iBod)-tStart,theta(iBod));
ylabel(ax1(1),'Speed (diameter/s)')
%ylabel(ax1(2),'Heading (deg)')

%xlabel('Time s)');
set(h(1),'Color','k','LineWidth',1.5)
set(ax1(1),'YColor','k')
% set(h(2),'Color',0.7.*[1 1 1],'LineWidth',1.5)
set(ax1,'TickDir','out')
%set(ax1(2),'YColor',0.5.*[1 1 1])
title(ttext)
% grid on
% Get x-limits
xL = xlim;

% Create axes
h1 = subplot(4,5,6:5:16);
h2 = subplot(4,5,[7:10 12:15 17:20]);

% Loop thru arms and plot seastar bodies
for i = 1:numArms
    axes(h1)
    h = fill(starPts(:,1),starPts(:,2)+i-0.5,sClr,'EdgeColor','none');
    axis equal
    hold on
    
    h = line([0 meanArmx(iArm(i))],[i meanArmy(iArm(i))+i]-0.5,...
        'Color',0.4.*[1 1 1],'LineWidth',2);

    drawnow
end
ylim([0 5])


axes(h2)

bClr = 0.95.*[1 1 1];

arms = 2:2:6;

% Loop thru every-other arm to create light bars in background
for i = 1:length(arms)
    
    % Current arm number
    n = arms(i);
    
    % Coords for bar
    xVals = [xL(1) xL(2) xL(2) xL(1) xL(1)];
    yVals = [n n n+hBar n+hBar n];

    % Make bar
    h = fill(xVals,yVals,bClr,'EdgeColor','none');
    hold on
end

if 1
% Loop thru each footprint
for i = 1:length(armPos)

    % Current arm 
    n = find(iArm==armNum(i),1,'first');
    
    % Current foot number
    nFt = armPos(i);

    % Micro bar height
    mhBar = 0.05;
    
    % Start and enf times for the power stroke
    tStartFt = frStart(i)./frRate - tStart;
    tEndFt   = frEnd(i)./frRate - tStart;
    tDurFt   = tEndFt - tStartFt;
    
    % Coords for bar
    yVals = [n n n+mhBar n+mhBar]';
    
    % Offset y by foot position along arm
    yVals = yVals + (armPos(i)-min(armPos))./range(armPos);
       
    % Number of divisions in each bar
    nDiv = 6;
    
    % Break up each interval into 4 parts, loop thru parts
    for j = 1:(nDiv)
       
        tVals(1,1) = tStartFt + tDurFt/nDiv*(j-1);
        tVals(2,1) = tStartFt + tDurFt/nDiv*(j);
        tVals(3,1) = tVals(2);
        tVals(4,1) = tVals(1);
        
        % Interpolate for angular values at ends of patch
        angVal1 = interp1(F(i).frames./frRate-tStart,angNorm{i},tVals(1));
        angVal2 = interp1(F(i).frames./frRate-tStart,angNorm{i},tVals(2));
        
        % Interpolate first color
        col(1,1) = interp1(cval,cmap(:,1),angVal1);
        col(1,2) = interp1(cval,cmap(:,2),angVal1);
        col(1,3) = interp1(cval,cmap(:,3),angVal1);
        
        % Interpolate second color
        col(2,1) = interp1(cval,cmap(:,1),angVal2);
        col(2,2) = interp1(cval,cmap(:,2),angVal2);
        col(2,3) = interp1(cval,cmap(:,3),angVal2);
        
        % Define other colors
        col(3,:) = col(2,:);
        col(4,:) = col(1,:);
        
        % Verticies
        v = [tVals yVals];
        
        % Draw patch
        patch('Faces',1:4,'Vertices',v,'FaceVertexCData',col,...
              'FaceColor','interp','EdgeColor','none','FaceAlpha',fAlpha);  
    end
end
end

grid on
ylim([1 6])

% Position of top graph
pos3 = get(ax1(1),'Position');

% Position of gait graph
pos2 = get(h2,'Position');

% Position of body graphs
pos1 = get(h1,'Position');

% Adjust position of gait graph
% set(h2,'Position',[pos2(1)*0.8 pos2(2)*1.3 pos2(3) pos1(4)*.88],...
%     'YColor','none','TickDir','out');
 set(h2,'YColor','none','TickDir','out');

% Adjust position of top graph
set(h1,'XColor','none','YColor','none','Position',[pos1(1)*1.4 pos1(2:4)]);

set(ax1,'XColor','none')

% set(ax1,'Position',[pos2(1)*0.8 pos3(2) pos2(3) pos3(4)])
xlabel('Time (s)')


figure;
subplot(4,5,2:5)
h2(1) = plot(Body.t(iBod)-tStart,OP,'k-');
set(gca,'TickDir','out')
ylabel('Polar order parameter')
xlabel('time (s)')
% grid on
ylim([0 1])
set(h2(1),'Color','k','LineWidth',1.5)

ttt=2






