function varargout = postProcess(pMode,varargin)
% Post-proccessing of tube foot kinematics


%% Process inputs

% Set mode-specific parameters
if strcmp(pMode,'find arms')
    
    Body = varargin{1};
    iC   = varargin{2};
    B_ft = varargin{3};
    
    
elseif strcmp(pMode,'add arms')
    
    Body = varargin{1};
    iC   = varargin{2};
    B2   = varargin{3};
    
elseif strcmp(pMode,'connect')
    
    Body        = varargin{1};
    iC          = varargin{2};
    B_ft        = varargin{3};
    dist_thresh = varargin{4};
   
elseif strcmp(pMode,'reorganize')
    
    Body        = varargin{1};
    iC          = varargin{2};
    B_ft        = varargin{3};
    dist_thresh = varargin{4};
  
elseif strcmp(pMode,'Traj body system')
    
    Body        = varargin{1};
    
    % Duration of walking considered the 'start' (s)
    startDur    = 30;
    
    
elseif strcmp(pMode,'Traj foot system')
    
    Body        = varargin{1};
    F           = varargin{2};

    startDur = Body.tStartT;
    
    % Distance from body to ground (wrt body diameter)
    bodHeight = 0.05;
    
    
elseif strcmp(pMode,'Remove bad feet')
    
    Body        = varargin{1};
    F           = varargin{2};

    startDur = Body.tStartT;
    
    % Minimum displacement in direction of traj
    minDisp = 30;
    
    % Minium duration of contact
    minFrames = 30;
 
    % Minimum distance from mouth, relative to body size
    minMouthDist = 0.05;
    
    
else
    error(['Do not recognize ' pMode])
end


%% Parameters

% Radial increment to partion membership to an arm
rPart = pi/5;

% Colormap values
cmap = cmap_vals;


%%  Find arm address and index for global points

if strcmp(pMode,'find arms')

    % Origin of roi
    xO = Body.Rotation.roi(1).rect(3)/2;
    yO = Body.Rotation.roi(1).rect(4)/2;
    
    % Loop thru arms
    for j = 1:5
        % Arm coordinates
        xA(j,1) = iC.xArms(j)-iC.x+xO;
        yA(j,1) = iC.yArms(j)-iC.y+yO;
        
        % Radial position of arm
        rA(j,1) = atan2(yA(j)-yO,xA(j)-xO);
    end
    
    % Body diameter
    bDiam = 2*max(hypot(iC.xArms-iC.x,iC.yArms-iC.y));
    
    % ASSIGN ARM NUMBERS ------------------------------------------------------
    
    % Loop thru all frames
    for i = 1:length(B_ft)
        
        B(i).fr_num = B_ft(i).fr_num;
        B(i).frIdx  = B_ft(i).frIdx;
        
        % Store properties of local blobs
        B(i).L = B_ft(i).propsL;
        
        % If there is data . . .
        if ~isempty(B_ft(i).frIdx)
            
            % Current transformation matrix
            tform = Body.Rotation.tform(i);
            
            % Current region of interest
            roi = Body.Rotation.roi(i);
            
            % Loop thru each blob
            for j = 1:length(B_ft(i).propsL)
                
                % Current local coordinates for foot
                xC = B_ft(i).propsL(j).Centroid(1);
                yC = B_ft(i).propsL(j).Centroid(2);
                rC = atan2(yC-yO,xC-xO);
                
                % Transform local to global points
                [xOr,yOr] = transformPointsInverse(tform,xC-roi.r,yC-roi.r);
                
                % Translate by center of roi
                xOr = xOr  + Body.xCntr(i);
                yOr = yOr  + Body.yCntr(i);
                
                % Store local values
                %B_ft(i).propsL(j).coordL = [xOr yOr];
                
                % Indicies for finding match to global
                iMatch  = 0;
                minDist = inf;
                
                % Loop thru global points
                for k = 1:length(B_ft(i).propsG)
                    
                    % Current distance
                    cDist = hypot(B_ft(i).propsG(k).Centroid(1)-xOr, ...
                                  B_ft(i).propsG(k).Centroid(2)-yOr);
                    
                    % Store index, if closer than previous
                    if cDist<minDist
                        % Update min distance
                        minDist = cDist;
                        
                        % Store index
                        iMatch = k;
                    end
                end                
                
                % Check by visualizing
                if 0
                    for k = 1:length(B_ft(i).propsG)
                        plot(B_ft(i).propsG(k).Centroid(1), ...
                             B_ft(i).propsG(k).Centroid(2),'ok');
                        hold on;
                    end
                    plot(B_ft(i).propsG(iMatch).Centroid(1), ...
                         B_ft(i).propsG(iMatch).Centroid(2),'r+');
                    plot(xOr,yOr,'ro');
                    axis equal; hold off
                end
                
                % Store global data that matches local data
                for fn = fieldnames(B_ft(i).propsG(iMatch))'
                    B(i).G(j,1).(fn{1}) = B_ft(i).propsG(iMatch).(fn{1});
                end
                
                % Matching coordinate
                xMatch = B_ft(i).propsG(iMatch).Centroid(1);
                yMatch = B_ft(i).propsG(iMatch).Centroid(2);
                
                cOffset(j,:) = [xOr-xMatch yOr-yMatch];
                
                % Store matching index from global coordinates
%                 B_ft(i).propsL(j).iG = iMatch;
                
                % Loop thru arms to find a match
                for k = 1:5
                    
                    % If within range of radial positions, assign arm number
                    if abs(angdiff2(rC,rA(k)))<rPart
                        
                        B(i).L(j).armNum = k;
                        B(i).G(j).armNum = k;
                        break
                    else
                        B(i).L(j).armNum = nan;
                        B(i).G(j).armNum = nan;
                    end
                end
            end
            
            B(i).cOffset = mean(cOffset);
            
            clear cOffset
        end
    end
    
    % Define output
    varargout{1} = B;
end


%% Add arm coordinates to Body

if strcmp(pMode,'add arms')

% Body.xArmL = iC.xArms;
% Body.yArmL = iC.yArms;

xArmL = Body.xArmL;
yArmL = Body.yArmL;
    
% Loop thru time
for i = 1:length(Body.frames)
    
    % Current roi and tform
    roi   = Body.Rotation.roi(i);
    tform = Body.Rotation.tform(i);
    
    % Current body center
    xCntr = Body.xCntr(i);
    yCntr = Body.yCntr(i);
    
    % Current offset
    if isempty(B2(i).cOffset)
        cOffset = [0 0];
    else
        cOffset = B2(i).cOffset;
    end
    
    % Find current global posiiton of arm tips
    [xArmsG,yArmsG] = local2global(xCntr,yCntr,cOffset,tform,roi,xArmL,yArmL);
    
    % Store result
    Body.xArmG(i,:) = xArmsG;
    Body.yArmG(i,:) = yArmsG;
    
    % Visual to test the transformation
    if 0
        % Loop trhu feet
        for j = 1:length(B2(i).G)
            % Current local feet 
            xFtL = B2(i).L(j).Centroid(1);
            yFtL = B2(i).L(j).Centroid(2);
            
            % Current global feet 
            xFtG = B2(i).G(j).Centroid(1);
            yFtG = B2(i).G(j).Centroid(2);
            
            % Attempt to match global
            [xFtG2,yFtG2] = local2global(xCntr,yCntr,cOffset,tform,roi,xFtL,yFtL);
            
            plot(xFtG,yFtG,'ko',xFtG2,yFtG2,'ro')
            hold on
            axis equal
            title(num2str(i))
        end
        plot(xArmsG,yArmsG,'k+')
        axis equal
        hold off
    end
end
    
% Define output
varargout{1} = Body;    
    
end


%% Connect points across frames

if strcmp(pMode,'connect')
    
    % Find starting frame
    for i = 1:length(B_ft)
        if ~isnan(B_ft(i).frIdx)
            iStart = i;
            break
        end
    end

    
    % Current frame's global properties
    currG = B_ft(iStart).G;
    
%     % Next frame's global properties
%     nextG = B_ft(i+1).G;
%   
%     % Find next points
%     [currG,nextG] = findNexts(currG,nextG,dist_thresh);
%     
%     % Advance current G
%     currG = nextG;
  
    % Loop thru time
    for i = iStart:(length(B_ft)-1)       
        
        % Next frame's global properties
        nextG = B_ft(i+1).G;
        
        % Find next points
        [currG,nextG] = findNexts(currG,nextG,dist_thresh);
        
        % Store results
        B_ft(i).G = currG;
        
        % Advance current G
        currG = nextG;
    end
 
    % Define output
    varargout{1} = B_ft;
    
    
    % Initialize used structure
    for i = 1:length(B_ft)
        for j = 1:length(B_ft(i).G)
            B_ft(i).G(j).used = 0;
        end
    end
    
    % Indicies for F
    iFoot  = 1;
    iFrame = 1;
    
    % Indicies for frame and foot in source
    i = iStart;
    j = 1;
    
    iColor = 1;
    
    while iStart<(length(B_ft)-1)
        
        % Store away coordinates and frame number
        F(iFoot).frames(iFrame,1) = B_ft(i).fr_num;
        F(iFoot).xG(iFrame,1)     = B_ft(i).G(j).Centroid(1);
        F(iFoot).yG(iFrame,1)     = B_ft(i).G(j).Centroid(2);
        F(iFoot).xL(iFrame,1)     = B_ft(i).L(j).Centroid(1);
        F(iFoot).yL(iFrame,1)     = B_ft(i).L(j).Centroid(2);
        F(iFoot).clr(iFrame,:)    = cmap(iColor,:);
        F(iFoot).armNum(iFrame,1) = B_ft(i).L(j).armNum;
        
        % Score as used
        B_ft(i).G(j).used = 1;
        
        % If index for next exists . . .
        if ~isnan(B_ft(i).G(j).iNext)
            % Advance index for frame on A
            iFrame = iFrame + 1;
            
            % Jump to index for next foot
            j = B_ft(i).G(j).iNext;
            
            % Advance frame on source
            i = i + 1;
            
        % Otherwise, move onto next foot
        else
            
            % Marker variable
            isfoot = 0;
            
            % Reset frame and foot indicies
            i = nan;
            j = nan;
            
            % Loop trhu frames
            for k = iStart:length(B_ft)
                % Loop thru feet
                for l = 1:length(B_ft(k).G)
                    
                    % If not used already, 
                    if ~B_ft(k).G(l).used
                        
                        % Update marker
                        isfoot = 1;
                        
                        % Update indicies from B_ft
                        i = k;
                        j = l;
                        
                        break
                    end
                end
                
                % If there is an unused foot, break
                if isfoot
                    break
                    
                % If no match starting at iStart
                else
                    
                    % advance to skip next time
                    iStart = k+1;                       
                end
            end
            
            % Advance indicies in A for foot
            iFoot  = iFoot + 1;
            iFrame = 1;
            
            % Advance color for plotting later
            if iColor==size(cmap,1)
                iColor = 1;
            else
                iColor = iColor + 1;
            end
            
            % Update
            %disp(['iFoot = ' num2str(iFoot)]) 
        end
        %disp(['iFoot = ' num2str(iFoot), ' iFrame = ' num2str(iFrame)]) 
    end
    
    % Define output
    varargout{1} = F;
end


%% Add coordinate for trajectory to Body

if strcmp(pMode,'Traj body system')
    
    % Include early linear fit of trajectory
    iDur = Body.t-min(Body.t) < startDur;
    
    % Linear fit to trajectory
    cX = polyfit(Body.t(iDur),Body.xCntr(iDur),1);
    cY = polyfit(Body.t(iDur),Body.yCntr(iDur),1);
    
    % Starting and end points
    pStart = [polyval(cX,min(Body.t(iDur))) polyval(cY,min(Body.t(iDur)))];
    pEnd   = [polyval(cX,max(Body.t)*1.1) polyval(cY,max(Body.t)*1.1)];
    
    % Cntr points in trajectory FOR
    CntrPtsT = transCoord2d('xax G2L',pStart,pEnd,[Body.xCntr Body.yCntr]);
    
    % Loop trhu arm points, define in T system
    for i = 1:size(Body.xArmL,2)
        
        % Arm points in trajectory FOR
        ArmPtsT = transCoord2d('xax G2L',pStart,pEnd,[Body.xArmG(:,i) Body.yArmG(:,i)]);
        
        % Store
        Body.xArmT(:,i)  = ArmPtsT(:,1);
        Body.yArmT(:,i)  = ArmPtsT(:,2);
        Body.xArmTL(:,i) = ArmPtsT(:,1) - CntrPtsT(:,1);
        Body.yArmTL(:,i) = ArmPtsT(:,2) - CntrPtsT(:,2);
        
        clear ArmPtsT
    end
    
    % Store body center data
    Body.tStartT   = startDur;
    Body.xT        = CntrPtsT(:,1);
    Body.yT        = CntrPtsT(:,2);
    Body.startT    = pStart;
    Body.endT      = pEnd;
    
    % Show trajectory fit
    if 0
        figure
        plot(CntrPnts(:,1),CntrPnts(:,2),'o',...
            [S.pStart(1) S.pEnd(1)],[S.pStart(2) S.pEnd(2)],'-k')
        axis equal
        hold on   
    end
    
    % Define output
    varargout{1} = Body;   
end


%% Add coordinate syestems for trajectory to feet

if strcmp(pMode,'Traj foot system')
    
    % Body size
    bSize = max([range(Body.xArmL) range(Body.yArmL)]);
    
    % Loop thru feet
    for i = 1:length(F)
    
        % Loop trhu frames
         for j = 1:length(F(i).frames)
            
             % Current frame
             cFrame = F(i).frames(j);
             
             % Current arm number
             armNum = F(i).armNum(j);
             
             % Index for frame in Body
             iFrame = find(Body.frames==cFrame,1,'first');
             
             % Center point
             cntrPnt = [Body.xCntr(iFrame) Body.yCntr(iFrame)];
             
             % Foot points in trajectory FOR
             ptsT = transCoord2d('xax G2L',Body.startT,...
                 Body.endT,[F(i).xG(j) F(i).yG(j)]);
             
             % Coords on body, along trajectory 
             ptsTL(1) = ptsT(1) - Body.xT(iFrame);
             ptsTL(2) = ptsT(2) - Body.yT(iFrame);
             
             % Store
             F(i).xT(j,1)  = ptsT(1);
             F(i).yT(j,1)  = ptsT(2);
             F(i).xTL(j,1) = ptsTL(1);
             F(i).yTL(j,1) = ptsTL(2);
             
             if ~isnan(armNum)
                 % Current arm point
                 armPnt = [Body.xArmG(iFrame,armNum) ...
                           Body.yArmG(iFrame,armNum)];
                 
                 % Foot points in arm FOR
                 ptsA = transCoord2d('xax G2L',cntrPnt,...
                     armPnt,[F(i).xG(j) F(i).yG(j)]);
                 
                 % Store
                 F(i).xA(j,1) = ptsA(1);
                 F(i).yA(j,1) = ptsA(2);
                 
             else
                 % Store
                 F(i).xA(j,1) = nan;
                 F(i).yA(j,1) = nan;
                 
             end 
         end
         
         % Base position taken as the mean
        xBase = mean(F(i).xTL);
        yBase = mean(F(i).yTL);
        
        % Distance from center
        ftDist = hypot(F(i).xTL-xBase,F(i).yTL-yBase).*sign(F(i).xTL-xBase);
        
        % Angular position of foot 
        F(i).ftAng(:,1) = atan2(bodHeight*bSize,ftDist);
        
        clear xBase yBase
    end
    
    % Define output
    varargout{1} = F;   
end


%% Remove erroneous feet

if strcmp(pMode,'Remove bad feet')
    
    % Body size
    bSize = max([range(Body.xArmL) range(Body.yArmL)]);
    
    % Start index
    k = 1;
    
    % Loop thru feet
    for i = 1:length(F)
     
         % Displacement in direction of trajectory 
         dispT = F(i).xTL(1) - F(i).xTL(end);
         
         % Position along arm
         armPos = mean(F(i).xA);
         
         % Position along the width of the arm
         %meanYA   = abs(mean(F(i).yA));

        % If displacement is above threshold
        if (dispT > minDisp) && (armPos>(minMouthDist*bSize)) && ...
                length(F(i).frames)>minFrames
            
             % Store contents of current F
                for fn = fieldnames(F(i))'
                    nF(k).(fn{1}) = F(i).(fn{1});
                end
            k = k + 1;
        end

         clear dispT armPos
    end
    
    
    % Visual check
    if 0
        
        for i = 1:length(nF)
            subplot(1,2,1)
            plot(nF(i).xA./bSize,nF(i).yA./bSize,'-',...
                 nF(i).xA(1)./bSize,nF(i).yA(1)./bSize,'o');
            axis equal
            hold on
            title('Arm system (A)')
            
            subplot(1,2,2)
            plot(nF(i).xTL,nF(i).yTL,'-',nF(i).xTL(1),nF(i).yTL(1),'o');
            axis equal
            hold on
            title('Local system, along T (TL)')
            
            % Displacement in direction of trajectory 
            dispT = nF(i).xTL(1) - nF(i).xTL(end);
            
            disp(['frames = ' num2str(length(nF(i).frames)) ...
                ', dispT =' num2str(dispT) ...
                ])
        end
    end
    
    % Define output
    varargout{1} = nF;   
end



%% Reorganize feet

% if strcmp(pMode,'reorganize')
%     
%     while 
%     
%     
% end



function [xG,yG] = local2global(xCntr,yCntr,cOffset,tform,roi,x,y)

% Transform local to global points
[xG,yG] = transformPointsInverse(tform,x-roi.r,y-roi.r);

% Translate by center of roi
xG = xG  + xCntr - cOffset(1);
yG = yG  + yCntr - cOffset(2);


function [G,Gn] = findNexts(G,Gn,dist_thresh)
% Identifies corresponding point in next frame for each current point

% Enter default values for 'used'
used = zeros(length(Gn),1);

% Loop thru blobs
for i = 1:length(G)    
    
    % Set default for next index
    G(i).iNext = nan;
    
    % Set iPrev to nan, if not defined
    if ~isfield(G(i),'iPrev') || isempty(G(i).iPrev)
        G(i).iPrev = nan;
    end
    
    % Current corresponding global coordinate
    posG = G(i).Centroid;

    % Initialize min distance
    minDist = inf;
    
    % Loop thru all next points to find closest distance
    for j = 1:length(Gn)
        
        % Centroid of points in next frame
        posGnext = Gn(j).Centroid;
        
        % Current distance
        cDist = hypot(posG(1)-posGnext(1),posG(2)-posGnext(2));
        
        % If current is lower than prior minimum, store index
        if G(i).armNum==Gn(j).armNum && cDist<minDist && ~used(j)
            % Upate candidate index 
            iNext = j;
            
            % Update min distance
            minDist = cDist;
        end
    end
    

    % If minimum is below threshold . . .
    if minDist<dist_thresh
        
        % Store index for match to next frame
        G(i).iNext = iNext;
        
        % The 'previous' index on the next is the current index
        Gn(iNext).iPrev = i;
        
        % Mark that coordinates on next have been used
        used(iNext) = 1;
    end
end


function  delta = angdiff2(alpha, beta)
   validateattributes(alpha, {'numeric'}, {'vector', 'nonsparse'}, 1);
   if nargin > 1 
      validateattributes(beta, {'numeric'}, {'vector', 'nonsparse', 'size', size(alpha)}, 2);
      alpha = beta - alpha;
   else
      alpha = diff(alpha);
   end
   delta = mod(alpha + pi, 2*pi) - pi;  %constrain to [-pi, pi[. will prefer -pi over +pi
   delta(delta == -pi & alpha >= 0) = pi;         %so force -pi to +pi (only when the original angle was positive)


function cmap = cmap_vals
% Colormap of line colors
cmap(1,:) = [0         0.4470    0.7410];
%cmap(2,:) = [0.8500    0.3250    0.0980];
cmap(7,:) = [252 191 126]./256;
cmap(3,:) = [0.9290    0.6940    0.1250];
%cmap(4,:) = [0.4940    0.1840    0.5560];
cmap(4,:) = [243 182 249]./256;
cmap(5,:) = [0.4660    0.6740    0.1880];
cmap(6,:) = [0.3010    0.7450    0.9330];
%cmap(7,:) = [0.6350    0.0780    0.1840];
cmap(7,:) = [242 155 155]./256;


%     if 0
%         figure
%         for i = 1:length(currG)
%             h = plot(currG(i).Centroid(1),currG(i).Centroid(2),'o');
%             set(h,'MarkerEdgeColor',0.8.*[1 1 1])
%             
%             hold on
%         end
%         
%         for i = 1:length(nextG)
%             h = plot(nextG(i).Centroid(1),nextG(i).Centroid(2),'o');
%             set(h,'MarkerEdgeColor',0.5.*[1 1 1])
%             hold on
%         end
%         
%         iClr = 1;
%         
%         for i = 1:length(currG)
%             
%             idx = currG(i).iNext;
%             
%             h = plot(currG(i).Centroid(1),currG(i).Centroid(2),'o');
%             set(h,'MarkerEdgeColor',cmap(iClr,:))
%             
%             h = plot(nextG(idx).Centroid(1),nextG(idx).Centroid(2),'o');
%             set(h,'MarkerEdgeColor',cmap(iClr,:))
%             
%             if iClr==size(cmap,1)
%                 iClr = 1;
%             else
%                 iClr = iClr + 1;
%             end
%             hold on
%         end
%         
%         %cmap(currG(i).armNum,:)
%         axis equal
%     end
    
    %         A(i).xF = nans(length(B_ft),1);



    