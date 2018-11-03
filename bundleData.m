function S = bundleData(iC,Centroid,Rotation,H,A,v)
% Bundles together the various elements of the acquired data and reports on
% its characteristics


%% Parameters

% Number of points to define the region of interest
numroipts = 500;

% Starting duration (for fitting trajectory
S.startDur = 40;


%% Initial variables

% Region of interest for first frame
roi0 = giveROI('define','circular',numroipts,iC.r,iC.x,iC.y);

% Create coordinate transformation structure
S_t = defineSystem2d('roi',roi0,Centroid,Rotation);

if ~isfield(H,'ft') || (length(H.ft)==1 && sum(~isnan(H.ft(1).xBase))==0)
    error('There are no manual points yet selected')
end

% Find bounds of frame numbers
frMin = max([H.frames(1)  S_t.frames(1)]);
frMax = min([H.frames(end) S_t.frames(end)]);

% Definitive frame vector in common btwn Centroid and Rotation
frames = frMin:frMax;

% Report frame interval
disp(['Start frame: ' num2str(frames(1)) '    End frame: ' num2str(frames(end))])


%% Adjust H structure

unEqual = 0;

% Loop thru feet
for i = 1:length(H.ft)
    % If mismatch in length of foot data and frames
    if length(H.ft(1).xBase)~=length(H.frames)
        warning(['Foot ' num2str(i) ' does not have data of the same length as the number of frames'])
        unEqual = 1;        
%         H.ft(i).xBase = H.ft(i).xBase(H.frames);
%         H.ft(i).yBase = H.ft(i).yBase(H.frames);
%         H.ft(i).xTip  = H.ft(i).xTip(H.frames);
%         H.ft(i).yTip  = H.ft(i).yTip(H.frames);
    end
end

if unEqual==1
   error('See if these unequal numbers can be corrected') 
end


%% Analyze trajectory

% Add time vector
S_t.t = (S_t.frames-S_t.frames(1)) ./ v.FrameRate;

% Include all points in linear fit of trajectory
%iDur = 1:length(S_t.xCntr);

% Include early linear fit of trajectory
iDur = S_t.t < S.startDur;

% Linear fit to trajectory
cX = polyfit(S_t.t(iDur),S_t.xCntr(iDur),1);
cY = polyfit(S_t.t(iDur),S_t.yCntr(iDur),1);

% Starting and end points
S.pStart = [polyval(cX,min(S_t.t(iDur))) polyval(cY,min(S_t.t(iDur)))];
S.pEnd   = [polyval(cX,max(S_t.t)*1.1) polyval(cY,max(S_t.t)*1.1)];

% Cntr points in global FOR
CntrPnts = [S_t.xCntr S_t.yCntr];

% Cntr points in trajectory FOR
CntrPtsT = transCoord2d('xax G2L',S.pStart,S.pEnd,CntrPnts);

% Show trajectory fit
if 0
    figure
    plot(CntrPnts(:,1),CntrPnts(:,2),'o',...
         [S.pStart(1) S.pEnd(1)],[S.pStart(2) S.pEnd(2)],'-k')
    axis equal 
    hold on
    
end

clear cX cY CntrPnts


%% Translate body center points into local and trajectory FOR

% Loop trhu common frames
for j = 1:length(frames)
    
    % Index for current frame in S_t data
    iS_t = find(S_t.frames==frames(j),1,'first');
    
    if isempty(iS_t)
        error('No match to S strcuture!')
    end
    
    % Copy over from S_t to S
    S.frames(j)    = frames(j);
    S.xCntr(j)     = S_t.xCntr(iS_t);
    S.yCntr(j)     = S_t.yCntr(iS_t);
    S.xCntrT(j)    = CntrPtsT(iS_t,1);
    S.yCntrT(j)    = CntrPtsT(iS_t,2);
    S.ang(j)       = S_t.ang(iS_t);
    S.roi(j)       = S_t.roi(iS_t);
    
    clear iS_t
    
    % Calcuate tform from rotation angle
    cOrigin = [S.xCntr(j) S.yCntr(j)];
    xAxis   = [cOrigin(1)+cos(S.ang(j)/180*pi) cOrigin(2)+sin(S.ang(j)/180*pi)];
    tform = defineSystem2d('x-axis',cOrigin,xAxis);
    S.tform(j) = tform.tform;  
    
    clear cOrigin xAxix tform
end

clear CntrPntsT 

%% Add time vector

S.t = (S.frames-S.frames(1)) ./ v.FrameRate;


%% Translate arm coorinates into local and traj FOR

% Loop thru frames with arm coordinates
for j = 1:length(A.frames)
    
    % Index of current frame in H data
    iFrameH = find(H.frames==A.frames(j),1,'first');
    
    % Index of current frame in S data
    iFrameS = find(S.frames==A.frames(j),1,'first');
    
    if isempty(iFrameS)
        error('Frame numbers for arms does not fall in range of centroid data')
    end
    
    % Get origin
    cOrigin = [S.xCntr(iFrameS) S.yCntr(iFrameS)];
    
    % Loop trhu arms
    for k = 1:length(A.arm)
        
        % Global coordinates for the arms in current frame
        armG(k,:) = [A.arm(k).x(j) A.arm(k).y(j)];
        
    end
    
    % Get local coordinates for arms
    %armL(:,:,j) = G2L(cOrigin,-S.ang(iFrameS)/180*pi,armG);
    armL(:,:,j) = transCoord2d('ang G2L',cOrigin,-S.ang(iFrameS)/180*pi,armG);
 
    clear armG iFrameH iFrameS cOrigin
end

% Find mean arm position
armL = mean(armL,3);

% Step thru frames to transform arms in local FOR into global FOR
for i = 1:length(S.xCntr)
    
    % Current origin
    cOrigin = [S.xCntr(i) S.yCntr(i)];
    
    % Calcuate the global position of the base
    %armG = G2L(cOrigin,-S.ang(i)/180*pi,armL);
    armG = transCoord2d('ang L2G',cOrigin,-S.ang(i)/180*pi,armL);
    
    % Arm points in trajectory FOR
    armT = transCoord2d('xax G2L',S.pStart,S.pEnd,armG);
    
    % Find trajectory points in current local FOR
    pntL = transCoord2d('ang G2L',cOrigin,-S.ang(i)/180*pi,[S.pStart;S.pEnd]-S.pStart+cOrigin);
    
    % x-axis point, relative to the body center in local FOR
    pntTL(i,:) = pntL(2,:);
    
    % Store coordinates
    S.arm(i).x    = armG(:,1);
    S.arm(i).y    = armG(:,2);
    S.armT(i).x   = armT(:,1);
    S.armT(i).y   = armT(:,2);
    
    clear cOrigin armG armT pntL
end

% Store local mean coordinates
S.armL.x = armL(:,1);
S.armL.y = armL(:,2);

% Include early linear fit of trajectory
iDur = S.t < S.startDur;

% Mean angle of the trajectory in local FOR
meanAng = mean(atan2(pntTL(iDur,2),pntTL(iDur,1)));

% Rotate local points wrt trajectory
tmp = transCoord2d('ang G2L',[0 0],meanAng,armL);

% Store rotated arm points
S.armLT.x = tmp(:,1);
S.armLT.y = tmp(:,2);

clear armL pntTL tmp armL


%% Translate tip and base points into local and trajectory FOR

% Loop thru tube feet
for k = 1:length(H.ft)
    
    % Loop trhu common frames
    for j = 1:length(frames)

        % Index for current frame for manual data
        iH = find(H.frames==frames(j),1,'first');
        
        % Transfer data
        S.ft(k).xBase(j,1)   = H.ft(k).xBase(iH);
        S.ft(k).yBase(j,1)   = H.ft(k).yBase(iH);
     
        % If base point is a nan . . .
        if isnan(S.ft(k).xBase(j))

            S.ft(k).xBaseL(j,1)  = nan;
            S.ft(k).yBaseL(j,1)  = nan;
            S.ft(k).xBaseLT(j,1) = nan;
            S.ft(k).yBaseLT(j,1) = nan;
            
        % If there is a tip point
        else
            
            % Current origin
            cOrigin = [S.xCntr(j) S.yCntr(j)];
            
            % Current base point
            baseG = [S.ft(k).xBase(j) S.ft(k).yBase(j)];
            
            % Base points in trajectory FOR
            baseT = transCoord2d('xax G2L',S.pStart,S.pEnd,baseG);

            % Local base point
            baseL = transCoord2d('ang G2L',cOrigin,-S.ang(j)/180*pi,baseG);
            
            % Local tip point wrt traj
            baseLT = transCoord2d('ang G2L',[0 0],meanAng,baseL);
            
            % Store local coordinates
            S.ft(k).xBaseL(j,1)   = baseL(1);
            S.ft(k).yBaseL(j,1)   = baseL(2);
            S.ft(k).xBaseT(j,1)   = baseT(1);
            S.ft(k).yBaseT(j,1)   = baseT(2);
            S.ft(k).xBaseLT(j,1)  = baseLT(1);
            S.ft(k).yBaseLT(j,1)  = baseLT(2);
            
            clear baseL baseLT
        end
            
        % Transfer data
        S.ft(k).xTip(j,1)    = H.ft(k).xTip(iH);
        S.ft(k).yTip(j,1)    = H.ft(k).yTip(iH);
     
        % If tip point is a nan . . .
        if isnan(S.ft(k).xTip(j))
            
            S.ft(k).xTipL(j,1)   = nan;
            S.ft(k).yTipL(j,1)   = nan;
            S.ft(k).xTipLT(j,1)  = nan;
            S.ft(k).yTipLT(j,1)  = nan;
            
        % If there is a tip point
        else
            
            % Current origin
            cOrigin = [S.xCntr(j) S.yCntr(j)];
            
            % Current tip point
            tipG = [S.ft(k).xTip(j) S.ft(k).yTip(j)];           
            
            % Local tip point
            %tipL = G2L(cOrigin,-S.ang(j)/180*pi,tipG);
            tipL = transCoord2d('ang G2L',cOrigin,-S.ang(j)/180*pi,tipG);
            
            % Tip points in trajectory FOR
            tipT = transCoord2d('xax G2L',S.pStart,S.pEnd,tipG);          
            
            % Local tip point wrt traj
            tipLT = transCoord2d('ang G2L',[0 0],meanAng,tipL);
            
            
            % Store local coordinates
            S.ft(k).xTipL(j,1)    = tipL(1);
            S.ft(k).yTipL(j,1)    = tipL(2);
            S.ft(k).xTipT(j,1)    = tipT(1);
            S.ft(k).yTipT(j,1)    = tipT(2);
            S.ft(k).xTipLT(j,1)   = tipLT(1);
            S.ft(k).yTipLT(j,1)   = tipLT(2);
            
            clear tipL tipG tipLT
        end

        clear iH
    end
end


%% Fill out base coordinates from position measurement

% Loop thru tube feet (write coordinates for base)
for j = 1:length(H.ft)
    
    % Find indicies for contact and release
    S.ft(j).iStart = find(~isnan(S.ft(j).xTip),1,'first');
    S.ft(j).iEnd   = find(~isnan(S.ft(j).xTip),1,'last');
    
    
    % If no foot number . . .
    if ~isfield(H,'ft') || sum(~isnan(H.ft(j).footNum))==0
        warning(['No foot number given for foot ' num2str(j) ', assuming 2'])
        S.ft(j).footNum      = 2;
        S.ft(j).footLet      = 'A';
    else
        S.ft(j).footNum      = H.ft(j).footNum;
        S.ft(j).footLet      = H.ft(j).footLet;
    end
    
    % If there's a base point
    if sum(~isnan(S.ft(j).xBase))>0
        
        % Find local position of base
        idx     = find(~isnan(S.ft(j).xBase),1,'first');
        baseL   = [S.ft(j).xBaseL(idx) S.ft(j).yBaseL(idx)];
        baseLT  = [S.ft(j).xBaseLT(idx) S.ft(j).yBaseLT(idx)];
        
        % Step thru frames to transform into global FOR
        for k = 1:length(S.xCntr)
            
            % Current origin
            cOrigin = [S.xCntr(k) S.yCntr(k)];
            
            % If there is a tip point
            if ~isnan(S.ft(j).xTip(k))
                
                % Calcuate the global position of the base
                %baseG = G2L(cOrigin,-S.ang(k)/180*pi,baseL);
                baseG = transCoord2d('ang L2G',cOrigin,-S.ang(j)/180*pi,baseL);
                
                % Overwrite the global position
                S.ft(j).xBase(k,1)   = baseG(1);
                S.ft(j).yBase(k,1)   = baseG(2);
                S.ft(j).xBaseL(k,1)  = baseL(1);
                S.ft(j).yBaseL(k,1)  = baseL(2);
                S.ft(j).xBaseLT(k,1) = baseLT(1);
                S.ft(j).yBaseLT(k,1) = baseLT(2);
                
                clear baseG base L base LT
            end
        end
        
        % If no base point . . .
    elseif sum(~isnan(S.ft(j).xTip))>0
        warning(['    Foot ' num2str(j) ' in does not have a base, ' ...
                 'but does have tip coordinates'])
    end
end


%% Organize data wrt leading arm

% Find angle of each arm wrt trajectory
armAng = atan2(S.armLT.y,S.armLT.x);

% Index of arms, ordered wrt heading
[tmp,iArm] = sort(abs(armAng));

% Sort arm points wrt lead arm
S.armLT.x    = S.armLT.x(iArm);
S.armLT.y    = S.armLT.y(iArm);
S.armL.x     = S.armL.x(iArm);
S.armL.y     = S.armL.y(iArm);
S.armAngLT   = armAng(iArm);

% Loop thru frames
for i = 1:length(S.arm) 
    
    % Reorder all arm points wrt lead arm
    S.arm(i).x   = S.arm(i).x(iArm);
    S.arm(i).y   = S.arm(i).y(iArm);
    S.armT(i).x  = S.armT(i).x(iArm);
    S.armT(i).y  = S.armT(i).y(iArm);
    
end

clear tmp armAng iArm

% Loop thru feet, find arm number for each
for i = 1:length(S.ft)
        
     % Index for values
     idx = ~isnan(S.ft(i).xBase);
     
     % If there are coordinates . . .
     if sum(idx)>0

         % Angular positions of base
         baseAng = mean(atan2(S.ft(i).yBaseLT(idx),S.ft(i).xBaseLT(idx)));
         
         % Mean absolute difference between base and arm
         diffAng = abs(baseAng-S.armAngLT);
         
         % Index for arm belonging to current tube foot
         iCurrArm = find(diffAng==min(diffAng),1,'first');
         
         % Get points for current foot in trajectory FOR
         S.ft(i).armNum = iCurrArm;
         
         clear baseT baseAng
     end 
end


%% Visual check on calculations

% Visualize data
if 0
    
    aVal = 0.5;
    
    % Current video object
    v = defineVidObject(currVidPath);
    
    % Loop thru frames
    for j = 300:length(S.xCntr)
        
        % Get frame
        im = getFrame(currVidPath,v,S.frames(j));
        
        % Get roi image
        [im_roi,bw_mask] = giveROI('stabilized',im,S.roi(j),...
            0,S.tform(j));
        % Show frame
        subplot(1,2,1)
        imshow(im);
        hold on
        
        % Overlay data
        plot(S.xCntr(j),S.yCntr(j),'g+')
        plot(S.arm(j).x,S.arm(j).y,'go')
        plot(S.arm(j).x(1),S.arm(j).y(1),'g+')
        
        % Loop thru feet
        for k = 1:length(S.ft)
            % If there is a base point
            if ~isnan(S.ft(k).xBase(j))
                % Plot line between base and tip
                h = line([S.ft(k).xBase(j) S.ft(k).xTip(j)],...
                    [S.ft(k).yBase(j) S.ft(k).yTip(j)],'LineWidth',3,...
                    'Color',[1 0 0 aVal]);
                h = scatter(S.ft(k).xTip(j),S.ft(k).yTip(j),...
                    'MarkerEdgeColor','r','SizeData',100);
            end
        end
        
        hold off
           
        % Show roi
        subplot(1,2,2)
        imshow(im_roi);
        hold on
        
        % Arm points in roi
        armR(:,1) = S.armL.x + size(im_roi,1)/2;
        armR(:,2) = S.armL.y + size(im_roi,2)/2;
        
        % Overlay data
        plot(armR(:,1),armR(:,2),'go')
        plot(armR(1,1),armR(1,2),'g+')
        hold off

    end
end