function varargout = postProcess(pMode,varargin)
% Post-proccessing of tube foot kinematics


%% Process inputs

% Set mode-specific parameters
if strcmp(pMode,'find arms')
    
    Body = varargin{1};
    iC   = varargin{2};
    B_ft = varargin{3};
    
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
            
            % Loop thru each blob
            for j = 1:length(B_ft(i).propsL)
                
                % Current local coordinates for foot
                xC = B_ft(i).propsL(j).Centroid(1);
                yC = B_ft(i).propsL(j).Centroid(2);
                rC = atan2(yC-yO,xC-xO);
                
                % Current transformation matrix
                tform = Body.Rotation.tform(i);
                
                % Current region of interest
                roi = Body.Rotation.roi(i);
                
                % Trasform local to global points
                [xOr,yOr] = transformPointsInverse(tform,xC- roi.r,yC-roi.r);
                
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
                
                % Store matching index from global coordinates
%                 B_ft(i).propsL(j).iG = iMatch;
                
                % Loop thru arms to find a match
                for k = 1:5
                    
                    % If within range of radial positions, assign arm number
                    if abs(angdiff(rC,rA(k)))<rPart
                        
                        B(i).L(j).armNum = k;
                        B(i).G(j).armNum = k;
                        break
                    else
                        B(i).L(j).armNum = nan;
                        B(i).G(j).armNum = nan;
                    end
                end
            end
        end
    end
    
    % Define output
    varargout{1} = B;
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


%% Reorganize feet

% if strcmp(pMode,'reorganize')
%     
%     while 
%     
%     
% end




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



    