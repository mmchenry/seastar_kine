function varargout = postProcess(pMode,varargin)
% Post-proccessing of tube foot kinematics


%% Process inputs

% Set mode-specific parameters

if strcmp(pMode,'package data')
    Body   = varargin{1};
    iC     = varargin{2};
    bPath  = varargin{3};
    fName  = 'foot_blobs';
    
elseif strcmp(pMode,'find offset')
    
    Body   = varargin{1};
    iC     = varargin{2};
    bPath  = varargin{3};
    B2     = varargin{4};

elseif strcmp(pMode,'post rotation')
    
    Body      = varargin{1};
    v         = varargin{2};
    iC        = varargin{3};
    maxSize   = varargin{4};
    
elseif strcmp(pMode,'add arms')
    
    Body = varargin{1};
    iC   = varargin{2};
    B2   = varargin{3};
    
elseif strcmp(pMode,'connect')
    
    Body        = varargin{1};
    iC          = varargin{2};
    B2          = varargin{3};
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
    minFrames = 25;
 
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

% Maximum size of an image dimension (for downsampling)
% maxSize = 350;


%% Post-processing after rotataion tracking 

% function Body = rotationPostProcess(Body,v,iC)

if strcmp(pMode,'post rotation')
    
    % Define time
    Body.t = Body.frames'./v.FrameRate;
    
    % Define local coordinates for arms
%     Body.xArmL = iC.xArms' - Body.x(1) + Body.Rotation.roi(1).r;
%     Body.yArmL = iC.yArms' - Body.y(1) + Body.Rotation.roi(1).r;
    
    % Store raw Rotation data
    Body.Rotation_raw = Body.Rotation;
    
    % Loop thru frames
    for i = 1:length(Body.Rotation.tform)
        
        % Scaling factor for downsampling
        Body.Rotation.roi(i).imFactor = maxSize./Body.Rotation.roi(i).rect(3);
        
%         % Compensate tform for downsampling
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
        Body.xCntr(i,1)  = Body.x(i)-roi.r + xCntr;
        Body.yCntr(i,1)  = Body.y(i)-roi.r + yCntr;
        Body.xOther(i,1) = Body.x(i)-roi.r + xOther;
        Body.yOther(i,1) = Body.y(i)-roi.r + yOther;
        
        ttt=3;
        
        %     Body.xArmL(i,:)  = Body.x(i)-roi.r + xArm';
        %     Body.yArmL(i,:)  = Body.y(i)-roi.r + yArm';
        
        % Adjust roi coordinates
%         Body.Rotation.roi(i).rect(1) = Body.x(i)-roi.r+round(xOr);
%         Body.Rotation.roi(i).rect(2) = Body.y(i)-roi.r+round(yOr);
%         Body.Rotation.roi(i).xCntr   = Body.xCntr(i,1);
%         Body.Rotation.roi(i).yCntr   = Body.yCntr(i,1);
%         Body.Rotation.roi(i).xPerimG = Body.Rotation.roi(i).xPerimG + xOr;
%         Body.Rotation.roi(i).yPerimG = Body.Rotation.roi(i).yPerimG + yOr;
        
        clear xOr yOr xCntr yCntr
    end
    
    
    % Define output
    varargout{1} = Body;
end


%% Package data into B2

if strcmp(pMode,'package data')
    
    % LOAD DATA FROM B_ft files  ------------------------------------------
    disp(' ')
    disp('postProcess(package data)')
    disp('     Loading B_ft files into B2 structure . . .')
    
    % Get listing of blob data files
    a = blobList(bPath,'foot_blobs');
    
    % Initialize storage index
    n = 1;
    
    % Loop thru all frames
    for i = 1:length(a)   
        
        % Load B_ft for current frame
        load([bPath filesep a(i).name])

        % If there is data in B_ft . . .
        if ~isempty(B_ft.frIdx) && isfield(B_ft,'propsG') &&  ...
                isfield(B_ft,'propsL') && ~isempty(B_ft.propsG) && ...
                ~isempty(B_ft.propsL)
      
            % Remove unneeded local fields
            if isfield(B_ft.propsL,'PixelIdxList')
                B_ft.propsL = rmfield(B_ft.propsL,'PixelIdxList');
                B_ft.propsL = rmfield(B_ft.propsL,'PixelList');
                B_ft.propsL = rmfield(B_ft.propsL,'MajorAxisLength');
                B_ft.propsL = rmfield(B_ft.propsL,'MinorAxisLength');
            end
            
            % Remove unneeded global fields
            if isfield(B_ft.propsG,'PixelIdxList')
                B_ft.propsG = rmfield(B_ft.propsG,'PixelIdxList');
                B_ft.propsG = rmfield(B_ft.propsG,'PixelList');
                B_ft.propsG = rmfield(B_ft.propsG,'MajorAxisLength');
                B_ft.propsG = rmfield(B_ft.propsG,'MinorAxisLength');
            end
            
            % Current frame number
            B2(n).fr_num = a(i).frNum;
            
            % Store properties of local blobs
            B2(n).L = B_ft.propsL;
            
            % Store properties of global blobs
            B2(n).G = B_ft.propsG;
            
            % Remove B_ft data
            clear B_ft
            
            % Advance index
            n = n + 1;
        end
    end
    
    disp('                                                  . . . done.')
 
    clear n i a
    
    % Define output
    varargout{1} = B2;
end


%%  Find arm address and index for global points of feet


if strcmp(pMode,'find offset')
    
    % PAIR LOCAL AND GLOBAL COORDS ----------------------------------------
    disp(' ')
    disp('postProcess(find offset)')
%     disp('     Loading B_ft files into B2 structure . . .')
    
    % Loop thru frames
    for i = 1:length(B2)
        
        % Index for current frame in 'Body'
        iFrame = find(Body.frames == B2(i).fr_num,1);
        
        % Current transformation matrix
        tform = Body.Rotation.tform(iFrame);
        
        % Current region of interest
        roi = Body.Rotation.roi(iFrame);
        
        % Global center
        xCntrG = Body.x(iFrame);
        yCntrG = Body.y(iFrame);
        
        % Extract local and global coords
        L = B2(i).L;
        G = B2(i).G;
        
        % Number of feet to track
        nPts = min([length(G) length(L)]);
        
        % Loop thru local coordinates & transform into global FOR into L2
        for j = 1:length(L)
            
            % Local coordinates for current foot
            xC(j,1) = L(j).Centroid(1);
            yC(j,1) = L(j).Centroid(2);
            
            % Rotated to orientation of current frame
            [xOrT(j,1),yOrT(j,1)] = transformPointsInverse(tform,xC(j),yC(j));
            
            % Translate by center of roi
            L2(j).Centroid(1) = xOrT(j) - roi.r  + xCntrG;
            L2(j).Centroid(2) = yOrT(j) - roi.r  + yCntrG;
        end
        
        % Define containers, based on relative length
        if length(G)>=length(L2)
            Cs = L2;
            Cl = G;
            Gbig = 1;    
        else
            Cs = G;
            Cl = L2;
            Gbig = 0;
        end
        
        % Vectors specifying whether used
%         sUsed = zeros(length(Cs),1);
        

        % Matrix of coordinates for short data
        for j = 1:length(Cs)
            Xs(1,j) = Cs(j).Centroid(1);
            Ys(1,j) = Cs(j).Centroid(2);
        end
        Xs = repmat(Xs,length(Cl),1);
        Ys = repmat(Ys,length(Cl),1);
        
        % Matrix of coordinates for long data
        for j = 1:length(Cl)
            Xl(j,1) = Cl(j).Centroid(1);
            Yl(j,1) = Cl(j).Centroid(2);
        end
        Xl = repmat(Xl,1,length(Cs));
        Yl = repmat(Yl,1,length(Cs));
        
        if Gbig
            % Matrix of coordinates for long data
            for j = 1:length(L)
                XL(1,j) = L(j).Centroid(1);
                YL(1,j) = L(j).Centroid(2);
            end
            XL = repmat(XL,length(Cl),1);
            YL = repmat(YL,length(Cl),1);
        else
            % Matrix of coordinates for long data
            for j = 1:length(L)
                XL(j,1) = L(j).Centroid(1);
                YL(j,1) = L(j).Centroid(2);
            end
            XL = repmat(XL,1,length(Cs));
            YL = repmat(YL,1,length(Cs));
        end
        
        % Keep track of what is used
        mUsed = zeros(length(Cl),length(Cs))==1;
        
        % Distance matrix
        dists = hypot(Xs-Xl,Ys-Yl);
        
        % Shorage index
        n = 1;
        
        % Step thru each column (i.e. short data) in distance matrix
        for j = 1:size(dists,2)
            
            % Set used coorddinates to infinity
            dists(mUsed) = inf;
            
            % distances in current column
            cDist = dists(:,j);
            
            % Indicies of distances along column, in order
            [~,iDist] = sort(cDist);
            
            % Step thru each distance value in order
            for k = 1:length(iDist)
                
                % If current distance equals the min for the row . . .
                if ~isinf(cDist(iDist(k))) && cDist(iDist(k))==min(dists(iDist(k),:))
                    
                    % Store values for short and long data
                    Xval_s(n,1) = Xs(iDist(k),j);
                    Yval_s(n,1) = Ys(iDist(k),j);
                    Xval_l(n,1) = Xl(iDist(k),j);
                    Yval_l(n,1) = Yl(iDist(k),j);
                    Xval_L(n,1) = XL(iDist(k),j);
                    Yval_L(n,1) = YL(iDist(k),j);
                    
                    i_l(n,1)    = iDist(k);
                    
                    % Mark the two coordinates as used
                    mUsed(iDist(k),:)  = ones(1,size(dists,2));
                    mUsed(:,j)         = ones(size(dists,1),1);
                
                    % Advance storage index
                    n = n + 1;
                    
                    break
                end
            end
        end
        
        
                % Difference, relative to
%                 ddiff = abs(dists(j,:)-min(dists(j,:)));
                
%                 % Index of column with min value
%                 idx = find(ddiff==min(ddiff),1);
%                 
%                 if dists(j,idx)==min(dists(:,idx))
%                     
%                     
%                     % Mark the two coordinates as used
%                     mUsed(j,:)   = ones(1,size(dists,2));
%                     mUsed(:,idx) = ones(size(dists,1),1);
%                     
%                     n = n + 1;
%                 else
%                     
%                 end
%             end
%         end
        
        
%         % Loop thru all coords of shorter data
%         for j = 1:length(Cs)
%              
%             % Loop thru all coords for longer data
%             for k = 1:length(Cl)
%                 % Find distance bwtn coords
%                 cDist(k,1) = hypot(Cs(j).Centroid(1)-Cl(k).Centroid(1),...
%                                    Cs(j).Centroid(2)-Cl(k).Centroid(2));
%             end      
%             
%             % Set used values to infinate
%             cDist(lUsed==1) = inf;
%             
%             % Index for min distance
%             iMatch(j,1) = find(cDist==min(cDist),1);
%             
%             % Log used long point
%             lUsed(iMatch) = 1;
%         end
        
        % Transfer centroid values
%         for j = 1:length(Cs)
%             
%             % If G is the smaller . . .
%             if ~Gbig
%                 G3(j).Centroid   = [Xval_s(j) Yval_s(j)];
%                 L3(j).Centroid   = L(i_l(j)).Centroid;
%                 Lloc(j).Centroid = L2(i_l(j)).Centroid;
%                 
%             % If L is the smaller
%             else
%                 Lloc(j).Centroid = L2(j).Centroid;
%                 G3(j).Centroid   = G(i_l(j)).Centroid;
%                 L3(j).Centroid   = [Xval_s(j) Yval_s(j)];
%             end         
%         end
        
       for j = 1:length(Xval_L)
            
           L3(j).Centroid   = [Xval_L(j) Yval_L(j)];
           
            % If G is the smaller . . .
            if ~Gbig
                G3(j).Centroid   = [Xval_s(j) Yval_s(j)];
                Lloc(j).Centroid = [Xval_l(j) Yval_l(j)];
                
            % If L is the smaller
            else
                G3(j).Centroid   = [Xval_l(j) Yval_l(j)];
                Lloc(j).Centroid = [Xval_s(j) Yval_s(j)];
            end         
        end


        %clear Cl Cs L G L2
        
        % Loop thru coords
        for j = 1:length(L3)
            % Offset
            cOffset(j,:) = [G3(j).Centroid(1) - Lloc(j).Centroid(1)  ...
                            G3(j).Centroid(2) - Lloc(j).Centroid(2)];
            % Distance between points
            dDist(j,1)   = hypot(G3(j).Centroid(1) - Lloc(j).Centroid(1),  ...
                                 G3(j).Centroid(2) - Lloc(j).Centroid(2));
        end
        
        % Index for excluding offset outliers
        idx = dDist < quantile(dDist,0.9);
        
        % Store mean offset
        B2(i).cOffset = nanmean(cOffset(idx,:),1);
        
        if sum(sum(dDist > quantile(dDist,0.95)))>0
            idx = find(dDist > quantile(dDist,0.95));
            
            for j = 1:length(idx)
                G3(idx(j)).Centroid    = [nan nan];
                L3(idx(j)).Centroid    = [nan nan];
                Lloc(idx(j)).Centroid  = [nan nan];
            end
        end
                
        
        % Put coordinate back into B2
        B2(i).G = G3;
        B2(i).L = L3;

        % Check by visualizing
        if 0 %i==1 || i==800
            
            figure
            clrs = colormap('lines');
            
            % Loop thru feet
            for j = 1:length(B2(i).G)
                
                 % Local coordinates for current foot
                xG = B2(i).G(j).Centroid(1);
                yG = B2(i).G(j).Centroid(2);
                
                % Local coordinates for current foot
                xC = B2(i).L(j).Centroid(1);
                yC = B2(i).L(j).Centroid(2);
                
                % Rotated to orientation of current frame
                [xOrT,yOrT] = transformPointsInverse(tform,xC,yC);
                
                % Translate by center of roi
                xOrT = xOrT - roi.r  + xCntrG + B2(i).cOffset(1);
                yOrT = yOrT - roi.r  + yCntrG + B2(i).cOffset(2);
                
                % Cuyrrent color
                cclr = clrs(ceil(rand(1)*size(clrs,1)),:);
                
                % Plot
                plot(xG,yG,'o','Color',cclr);
                hold on
                plot(xOrT,yOrT,'+','Color',cclr);
                plot([xG xOrT],[yG yOrT],'-','Color',cclr);
                axis equal
                title(num2str(i))
            end
            ttt=3;
        end
        
        clear LG nPts Gbig xC yC xOrT yOrT tform roi iFrame j idx cOffset
        clear cDist dDist Cl Cs L G L2 ttt xCntrG yCntrG G3 L3 Lloc
        clear XL Xl Xs Xval_L Yval_L Xval_s Yval_s Xval_l Yval_l YL Yl Ys
        clear dists i_l k mUsed cclr clrs 
        tt=3;
    end
    
    
    
%     % FIND OFFSET FOR COORD. TRANS.  --------------------------------------
%     
%     disp(' ')
%     disp('     Finding offset . . .')
%     
%     % Loop thru all frames
%     for i = 1:length(B2)
%         
%         % Index for current frame in 'Body'
%         iFrame = find(Body.frames == B2(i).fr_num,1);
%         
%         % Current transformation matrix
%         tform = Body.Rotation.tform(iFrame);
%         
%         % Current region of interest
%         roi = Body.Rotation.roi(iFrame);
%         
%         % Global center
%         xCntrG = Body.x(iFrame);
%         yCntrG = Body.y(iFrame);
%         
%         % Local center
%         xCntrL = roi.r;
%         yCntrL = roi.r;
%         
%         % Loop thru each local blob
%         for j = 1:length(B2(i).L)
%             
%             % Coordinates for current foot
%             xC(j,1) = B2(i).L(j).Centroid(1);
%             yC(j,1) = B2(i).L(j).Centroid(2);
%             
%             % Radial position for foot
%             rC(j,1) = atan2(yC(j)-iC.r,xC(j)-iC.r);
%             
%             % Rotate local to global points
%             [xOrT,yOrT] = transformPointsInverse(tform,xC(j),yC(j));
%             
%             % Translate by center of roi
%             xOrT = xOrT - roi.r  + xCntrG;
%             yOrT = yOrT - roi.r  + yCntrG;
%             
%             % If there are global points
%             if ~isempty(B2(i).G)
%                 
%                 % Indicies for finding match to global
%                 iMatch  = 0;
%                 minDist = inf;
%                 
%                 % Loop thru global points
%                 for k = 1:length(B2(i).G)
%                     
%                     % Current distance
%                     cDist = hypot(B2(i).G(k).Centroid(1)-xOrT, ...
%                                   B2(i).G(k).Centroid(2)-yOrT);
%                     
%                     % Store index, if closer than previous
%                     if cDist<minDist
%                         % Update min distance
%                         minDist = cDist;
%                         
%                         % Store index
%                         iMatch = k;
%                     end
%                 end
%                 
%                 % Matching coordinate
%                 xOr(j,1)    = xOrT;
%                 yOr(j,1)    = yOrT;
%                 xMatch(j,1) = B2(i).G(iMatch).Centroid(1);
%                 yMatch(j,1) = B2(i).G(iMatch).Centroid(2);
%                 
%             else
%                 xOr(j,1)    = nan;
%                 yOr(j,1)    = nan;
%                 xMatch(j,1) = nan;
%                 yMatch(j,1) = nan;
%             end
%             
%             clear xOrT yOrT
%         end
%         
%         if exist('xC','var') && exist('xMatch','var')
%             
%             % All offsets
%             cOffset = hypot(xOr-xMatch,yOr-yMatch);
%             
%             % Index for excluding outliers
%             idx = cOffset < quantile(cOffset,0.9);
%             
%             % Store mean offset
%             B2(i).cOffset = nanmean(cOffset(idx),1);
% 
%             % Check by visualizing
%             if 0
%                 subplot(2,1,1)
%                 plot(xC,yC,'ko',iC.r,iC.r,'ro')
%                 title(['Local FOR (frame ' num2str(B2(i).fr_num) ')'])
%                 axis equal
%                 
%                 subplot(2,1,2)
%                 plot(xOr,yOr,'ko',xMatch+cOffset(1),yMatch+cOffset(2),'r+')
%                 hold on
%                 
%                 for m = 1:length(xOr)
%                     plot([xOr(m) xMatch(m)],[yOr(m) yMatch(m)],'r-')
%                 end
%                 hold off
%                 axis equal
%                 ttt = 3;
%             end
%             
%         else
%             B2(i).cOffset = [0 0];
%         end  
%     end
%     
%     disp('                                                  . . . done.')
%     
    % Define output
    varargout{1} = B2;
end


%% Add arm coordinates to Body, find arm assignments

if strcmp(pMode,'add arms')
    
disp(' ')
disp('postProcess(add arms) . . .')    
 
% Set up placeholders
Body.xArmG = nan(length(Body.frames),5);
Body.yArmG = nan(length(Body.frames),5);

% Position of arms in local FOR
xArmsL = iC.xArms - iC.x + iC.r;
yArmsL = iC.yArms - iC.y + iC.r;
    
% Angular positions of arms
angArms0 = atan2(iC.yArms - iC.y,iC.xArms - iC.x);

% Loop thru time
for i = 1:length(B2)
    
    % Index for current frame
    iFrame = find(B2(i).fr_num == Body.frames,1);
    
    % Current roi and tform
    roi   = Body.Rotation.roi(iFrame);
    tform = Body.Rotation.tform(iFrame);
    
    % Current body center
    xCntr = Body.x(iFrame);
    yCntr = Body.y(iFrame);
    
    % Current offset
    if isempty(B2(i).cOffset) || isnan(B2(i).cOffset(1))
        cOffset = [0 0];
    else
        cOffset = B2(i).cOffset;
    end
    
    % Rotate arm points from starting orientation to current frame
    [xArmsG,yArmsG] = transformPointsInverse(tform,xArmsL,yArmsL);
    
    % Translate by center of roi
    xArmsG = xArmsG - roi.r  + xCntr;
    yArmsG = yArmsG - roi.r  + yCntr;
    
    % Store result
    Body.xArmG(iFrame,:) = xArmsG;
    Body.yArmG(iFrame,:) = yArmsG;
    
    % Visual to test the transformation
    if 0 %B2(i).fr_num==9 || B2(i).fr_num==800
        
        figure;
        xFtL = []; yFtL = [];
        xFtG = []; yFtG = [];
        
        % Loop trhu local blobs
        for j = 1:length(B2(i).L)
            xL(j,1)    = B2(i).L(j).Centroid(1);
            yL(j,1)    = B2(i).L(j).Centroid(2);
        end
        
        % Loop trhu global blobs
        for j = 1:length(B2(i).G)
            xG(j,1)    = B2(i).G(j).Centroid(1);
            yG(j,1)    = B2(i).G(j).Centroid(2);
        end
        
        % Attempt to match global
%         [xFtG2,yFtG2] = local2global(xCntr,yCntr,cOffset,tform,roi,xFtL,yFtL);
        
        subplot(2,1,1)
        plot(xL,yL,'ko')
        hold on
        plot(xArmsL,yArmsL,'k+',iC.r,iC.r,'r+')
        axis equal
        title(['Frame ' num2str(B2(i).fr_num)])       
        axis equal
        hold off

        subplot(2,1,2)
        plot(xG,yG,'ko')
        hold on
        axis equal
        title(num2str(i))
        
        plot(xArmsG,yArmsG,'k+',xCntr,yCntr,'r+')
        axis equal
        hold off
    end
    
    % Angular position of arms in current frame         
     angArms = atan2(yArmsG-yCntr,xArmsG-xCntr);
     
     
    % Loop thru feet in local FOR
    if ~isempty(B2(i).L)
        for j = 1:length(B2(i).L)
            % Angular position of foot
            angFt = atan2(B2(i).L(j).Centroid(2) - iC.r,...
                          B2(i).L(j).Centroid(1) - iC.r);
                      
            % Diff in ang wrt arms
            for k = 1:length(angArms) 
                dAng(k,1) =  abs(angdiff(angFt,angArms0(k)));
            end
   
            % Store arm number
            B2(i).L(j).armNum = find(dAng==min(dAng),1);
            
            % Store angular position
            B2(i).L(j).ang = angFt;
        end
    end
    

    % Loop thru feet in global FOR
    if ~isempty(B2(i).G)
        for j = 1:length(B2(i).G)
            % Angular position of foot
            angFt = atan2(B2(i).G(j).Centroid(2) - yCntr,...
                          B2(i).G(j).Centroid(1) - xCntr);
            
            % Diff in ang wrt arms          
            for k = 1:length(angArms)
                dAng(k,1) =  abs(angdiff(angFt,angArms(k)));
            end
            
            % Store arm number
            B2(i).G(j).armNum = find(dAng==min(dAng),1);
            
            % Store angular position
            B2(i).G(j).ang = angFt;
        end
    end
    
    % Visual test of arm assignment
    if 0 %B2(i).fr_num==9 || B2(i).fr_num==800
        figure
        
        clrs = colormap('lines');
        
        % Loop trhu local blobs
        for j = 1:length(B2(i).L)
            xL(j,1)    = B2(i).L(j).Centroid(1);
            yL(j,1)    = B2(i).L(j).Centroid(2);
            angL(j,1)  = B2(i).L(j).ang;
            armL(j,1)  = B2(i).L(j).armNum; 
        end
        
        % Loop trhu global blobs
        for j = 1:length(B2(i).G)
            xG(j,1)    = B2(i).G(j).Centroid(1);
            yG(j,1)    = B2(i).G(j).Centroid(2);
            angG(j,1)  = B2(i).G(j).ang;
            armG(j,1)  = B2(i).G(j).armNum; 
        end
         
        % Loop thru and plot arms
        for j = 1:5
            
            subplot(2,1,1)
            iL = armL==j;
            plot(xL(iL),yL(iL),'o','Color',clrs(j,:))
            hold on
            plot([iC.r xArmsL(j)],[iC.r yArmsL(j)],'-','Color',clrs(j,:))
            axis equal
            title(['Frame ' num2str(B2(i).fr_num)])
            
            subplot(2,1,2)
            iG = armG==j;
            plot(xG(iG),yG(iG),'o','Color',clrs(j,:))
            hold on
            plot([xCntr xArmsG(j)],[yCntr yArmsG(j)],'-','Color',clrs(j,:))
        end
        
        ttt=2;
    end
end

% Define output
varargout{1} = Body;   
varargout{2} = B2; 
   
disp('                                                  . . . done.')
end


%% Connect points across frames

if strcmp(pMode,'connect')
    
    % Loop thru time
    for i = 1:(length(B2)-1)       
        
        % Current frame's global properties
        currG = B2(i).G;
    
        % Next frame's global properties
        nextG = B2(i+1).G;
        
        if isfield(currG,'Centroid') && isfield(nextG,'Centroid')
            % Find next points
            currG = findNexts(currG,nextG,dist_thresh);
            
            % Store results
            B2(i).G = currG;
        end  
    end

    % Initialize used structure
    for i = 1:length(B2)
        for j = 1:length(B2(i).G)
            B2(i).G(j).used = 0;
        end
    end
    
    % Indicies for F
    iFoot  = 1;
    iFrame = 1;
    
    % Indicies for frame and foot in source
    i = 1;
    
    % Index for foot in G field of B2 structure
    j = 1;
    
    % Index for starting frame at current foot
    iStart = 1;
    
    iColor = 1;
    
    while iStart<(length(B2)-1)
        
        % If there are centroid points for current selecttion, store them
        if ~isempty(B2(i).G) && isfield(B2(i).G(j),'Centroid') && ...
                ~isnan(B2(i).G(j).Centroid(1))
            
            % Store away coordinates and frame number
            F(iFoot).frames(iFrame,1) = B2(i).fr_num;
            F(iFoot).xG(iFrame,1)     = B2(i).G(j).Centroid(1);
            F(iFoot).yG(iFrame,1)     = B2(i).G(j).Centroid(2);
            F(iFoot).xL(iFrame,1)     = B2(i).L(j).Centroid(1);
            F(iFoot).yL(iFrame,1)     = B2(i).L(j).Centroid(2);
            F(iFoot).clr(iFrame,:)    = cmap(iColor,:);
            F(iFoot).armNum(iFrame,1) = B2(i).G(j).armNum;
        else
            tt=2;
        end
        
        % Score as used
        B2(i).G(j).used = 1;
        
        % If index for next exists . . .
        if ~isempty(B2(i).G) && isfield(B2(i).G(j),'iNext') && ...
                ~isnan(B2(i).G(j).iNext)
            
            % Advance index for frame on A
            iFrame = iFrame + 1;
            
            % Jump to index for next foot
            j = B2(i).G(j).iNext;
            
            % Advance frame on source
            i = i + 1;
            
        % Otherwise, move onto next foot
        else
            
            % Marker variable
            isfoot = 0;
            
            jLast = j;
            
            % Reset frame and foot indicies
            i = nan;
            j = nan;
            
            % Loop trhu frames
            for k = iStart:length(B2)
                
                % Loop thru feet
                for l = 1:length(B2(k).G)
                    
                    % If not used already,
                    if ~B2(k).G(l).used
                        
                        % Update marker
                        isfoot = 1;
                        
                        % Update indicies from B2
                        i = k;
                        j = l;
                        
                        break
                    end
                end
               
                % If there is an unused foot, break
                if isfoot || (i>length(B2)-1)
                    
%                     % If new foot
%                     if j~=jLast
%                         disp(['   Foot = ' num2str(j)])
%                     end
                    
%                     % If next iStart . . .
%                     if k~=iStart
%                         % Update status
%                         disp(['iStart = ' num2str(iStart) ', ' num2str((length(B2)-1)) ...
%                             ' frames'])
%                     end
                    
                    % advance to skip next time
                    iStart = k;
                    
                    % Break loop trhu frames
                    break
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
        end
        if isnan(i) && isnan(j)
            error('nans here for frame and foot');
        end
        %disp(['iFoot = ' num2str(iFoot), ' iFrame = ' num2str(iFrame)])
%         [i j iFoot iFrame]
    end
    
    
    if exist('F','var')
        % Define output
        varargout{1} = F;
    else
        varargout{1} = [];
    end
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
    
    Body.xArmT = nan(length(Body.frames),5);
    Body.yArmT = nan(length(Body.frames),5);
    Body.xArmTL = nan(length(Body.frames),5);
    Body.yArmTL = nan(length(Body.frames),5);
    
    % Loop trhu arm points, define in T system
    for i = 1:size(Body.xArmL,2)
        
        % Arm points in trajectory FOR
        ArmPtsT = transCoord2d('xax G2L',pStart,pEnd,...
                               [Body.xArmG(:,i) Body.yArmG(:,i)]);
        
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
        plot(Body.xT,Body.yT,'ko')
        axis equal
        hold on   
        plot(Body.xT(iDur),Body.yT(iDur),'r+')
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
         
         if isfield(F(i),'xTL') & ~isempty(F(i).xTL)
             % Base position taken as the mean
             xBase = mean(F(i).xTL);
             yBase = mean(F(i).yTL);
             
             % Distance from center
             ftDist = hypot(F(i).xTL-xBase,F(i).yTL-yBase).*sign(F(i).xTL-xBase);
             
             % Angular position of foot
             F(i).ftAng(:,1) = atan2(bodHeight*bSize,ftDist);
         else
             F(i).ftAng(:,1) = nan;
         end
        
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
        
        if ~isempty(F(i).xTL)
            % Displacement in direction of trajectory
            dispT = F(i).xTL(1) - F(i).xTL(end);
            
            % Position along arm
            armPos = mean(F(i).xA);
            
            % Position along the width of the arm
            %meanYA   = abs(mean(F(i).yA));
            
            % If displacement is above threshold
%             if (dispT > minDisp) && (armPos>(minMouthDist*bSize)) && ...
%                     length(F(i).frames)>minFrames
                
            if length(F(i).frames)>minFrames
                
                % Store contents of current F
                for fn = fieldnames(F(i))'
                    nF(k).(fn{1}) = F(i).(fn{1});
                end
                k = k + 1;
            end
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
            ttt=3;
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


% function [L,G] = matchCoord(L,G,tform,roi,xCntrG,yCntrG)
% 
% % Number of feet to track
% nPts = min([length(G) length(L)]);
% 
% % Identify which is bigger
% if length(G)>=length(L)
%     Gbig = 1;
% else
%     Gbig = 0;
% end
% 
% % Loop thru local coordinates to find local version
% for i = 1:length(L)
%     
%     % Local coordinates for current foot
%     xC(i,1) = L(i).Centroid(1);
%     yC(i,1) = L(i).Centroid(2);
% 
%     % Rotated to orientation of current frame
%     [xOrT(i,1),yOrT(i,1)] = transformPointsInverse(tform,xC(i),yC(i));
% 
%     % Translate by center of roi
%     xOrT(i) = xOrT(i) - roi.r  + xCntrG;
%     yOrT(i) = yOrT(i) - roi.r  + yCntrG;
% 
%     
% end


function [xG,yG] = local2global(xCntr,yCntr,cOffset,tform,roi,x,y)

% Transform local to global points
[xG,yG] = transformPointsInverse(tform,x-roi.r,y-roi.r);
% [xG,yG] = transformPointsInverse(tform,x,y);

% Translate by center of roi
xG = xG  + xCntr - cOffset(1);
yG = yG  + yCntr - cOffset(2);


% function [G,Gn] = findNexts(G,Gn,dist_thresh)
% % Identifies corresponding point in next frame for each current point
% 
% % Enter default values for 'used'
% used = zeros(length(Gn),1);
% 
% % Loop thru blobs
% for i = 1:length(G)    
%     
%     % Set default for next index
%     G(i).iNext = nan;
%     
%     % Set iPrev to nan, if not defined
%     if ~isfield(G(i),'iPrev') || isempty(G(i).iPrev)
%         G(i).iPrev = nan;
%     end
%     
%     % Current corresponding global coordinate
%     posG = G(i).Centroid;
% 
%     % Initialize min distance
%     minDist = inf;
%     
%     % Loop thru all next points to find closest distance
%     for j = 1:length(Gn)
%         
%         if isfield(G(i),'Centroid') && isfield(Gn(j),'Centroid') && ...
%                 ~isnan(Gn(j).Centroid(1)) && ~isnan(G(i).Centroid(1))
%             % Centroid of points in next frame
%             posGnext = Gn(j).Centroid;
%             
%             % Current distance
%             cDist = hypot(posG(1)-posGnext(1),posG(2)-posGnext(2));
%             
%             % If current is lower than prior minimum, store index
%             if G(i).armNum==Gn(j).armNum && cDist<minDist && ~used(j)
%                 % Upate candidate index
%                 iNext = j;
%                 
%                 % Update min distance
%                 minDist = cDist;
%             end
%         else
%             posGnext = nan;
%             cDist    = inf;
%             minDist   = nan;
%         end
% 
%         
%     end
%     
% 
%     % If minimum is below threshold . . .
%     if minDist<dist_thresh
%         
%         % Store index for match to next frame
%         G(i).iNext = iNext;
%         
%         % The 'previous' index on the next is the current index
%         Gn(iNext).iPrev = i;
%         
%         % Mark that coordinates on next have been used
%         used(iNext) = 1;
%     end
% end


function G = findNexts(G,Gn,dist_thresh)


% Matrix of coordinates for current data
for i = 1:length(G)
    Xc(1,i) = G(i).Centroid(1);
    Yc(1,i) = G(i).Centroid(2);
end
Xc = repmat(Xc,length(Gn),1);
Yc = repmat(Yc,length(Gn),1);

% Matrix of coordinates for long data
for i = 1:length(Gn)
    Xn(i,1) = Gn(i).Centroid(1);
    Yn(i,1) = Gn(i).Centroid(2);
end
Xn = repmat(Xn,1,length(G));
Yn = repmat(Yn,1,length(G));

% Keep track of what is used
mUsed = zeros(length(Gn),length(G))==1;

% Distance matrix
dists = hypot(Xc-Xn,Yc-Yn);

% Step thru each column (i.e. Current data) in distance matrix
for i = 1:size(dists,2)
    
    % Set used coorddinates to infinity
    dists(mUsed) = inf;
    
    % distances in current column
    cDist = dists(:,i);
    
    % Indicies of distances along column, in order
    [~,iDist] = sort(cDist);
    
    ismatch = 0;
    
    % Step thru each distance value in order
    for k = 1:length(iDist)
        
        % If current distance equals the min for the row . . .
        if ~isinf(cDist(iDist(k))) && cDist(iDist(k))==min(dists(iDist(k),:))
            
            % Store values for short and long data
%             Xval_c(n,1) = Xc(iDist(k),i);
%             Yval_c(n,1) = Yc(iDist(k),i);
%             Xval_n(n,1) = Xn(iDist(k),i);
%             Yval_n(n,1) = Yn(iDist(k),i);
            
            % Current and next indicies
%             i_c(i,1)    = i;
            i_n(i,1)    = iDist(k);
            
            % Mark the two coordinates as used
            mUsed(iDist(k),:)  = ones(1,size(dists,2));
            mUsed(:,i)         = ones(size(dists,1),1);
            
            % Break, if there's a match
            ismatch = 1;
            break
        end  
    end
    
    % Fill with nans, if no match
    if ~ismatch
%         i_c(i,1)    = nan;
        i_n(i,1)    = nan;
    end
end

for i = 1:length(G)
   G(i).iNext = i_n(i);
end

% for i = 1:length(Gn)
%    Gn(i).iPrev= i_c;
% end



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


function [a,frNums] = blobList(fPath,fPrefix)

frNums = [];

% File listing
a = dir([fPath filesep fPrefix '*']);

% Loop trhu files
for i = 1:length(a)
    
    % Index of separator
    iSep = find(a(i).name=='_',1,'last');
    
    % Get frame number
    a(i).frNum = str2num(a(i).name((iSep+1):end-4));
    
    % Listing of frame numbers
    frNums = [frNums; a(i).frNum];
end


    