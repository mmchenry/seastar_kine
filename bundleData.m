function S = bundleData(iC,Centroid,Rotation,H,A,v)
% Bundles together the various elements of the acquired data and reports on
% its characteristics


%% Parameters

% Number of points to define the region of interest
numroipts = 500;


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
frMax = min([H.frames(end) length(H.ft(1).xBase) S_t.frames(end)]);

% Definitive frame vector in common btwn Centroid and Rotation
frames = frMin:frMax;

% Report frame interval
disp(['Start frame: ' num2str(frames(1)) '    End frame: ' num2str(frames(end))])



%% Adjust H structure

% Loop thru feet
for i = 1:length(H.ft)
    % If mismatch in length of foot data and frames
    if length(H.ft)~=length(H.frames)
        H.ft(i).xBase = H.ft(i).xBase(H.frames);
        H.ft(i).yBase = H.ft(i).yBase(H.frames);
        H.ft(i).xTip  = H.ft(i).xTip(H.frames);
        H.ft(i).yTip  = H.ft(i).yTip(H.frames);
    end
end


%% Translate tip and base points into local FOR

% Loop trhu common frames
for j = 1:length(frames)
    
    % Index for current frame in S_t data
    iS_t = find(S_t.frames==frames(j),1,'first');
    
    % Copy over from S_t to S
    S.frames(j)    = frames(j);
    S.xCntr(j)     = S_t.xCntr(iS_t);
    S.yCntr(j)     = S_t.yCntr(iS_t);
    S.ang(j)       = S_t.ang(iS_t);
    S.roi(j)       = S_t.roi(iS_t);
    
    clear iS_t
    
    % Calcuate tform from rotation angle
    cOrigin = [S.xCntr(j) S.yCntr(j)];
    xAxis   = [cOrigin(1)+cos(S.ang(j)/180*pi) cOrigin(2)+sin(S.ang(j)/180*pi)];
    tform = defineSystem2d('x-axis',cOrigin,xAxis);
    S.tform(j) = tform.tform;
    
    % Loop thru tube feet
    for k = 1:length(H.ft)
        
        % Index for current frame for manual data
        iH = find(H.frames==frames(j),1,'first');
        
        % Transfer data
        S.ft(k).xBase(j,1)   = H.ft(k).xBase(iH);
        S.ft(k).yBase(j,1)   = H.ft(k).yBase(iH);
        S.ft(k).xTip(j,1)    = H.ft(k).xTip(iH);
        S.ft(k).yTip(j,1)    = H.ft(k).yTip(iH);
        S.ft(k).footNum    = H.ft(k).footNum;
        S.ft(k).footLet    = H.ft(k).footLet;
        
        % If tip point is a nan . . .
        if isnan(S.ft(k).xTip(j))
            
            S.ft(k).xTipL(j,1)   = nan;
            S.ft(k).yTipL(j,1)   = nan;
            S.ft(k).xBaseL(j,1)  = nan;
            S.ft(k).yBaseL(j,1)  = nan;
            
            % If there is a tip point
        else
            % Current tip point
            tipG = [S.ft(k).xTip(j) S.ft(k).yTip(j)];
            
            % Current base point
            baseG = [S.ft(k).xBase(j) S.ft(k).yBase(j)];
            
            % Local tip point
            %tipL = G2L(cOrigin,-S.ang(j)/180*pi,tipG);
            tipL = transCoord2d('ang G2L',cOrigin,-S.ang(j)/180*pi,tipG);
            
            % Local base point
            %baseL = G2L(cOrigin,-S.ang(j)/180*pi,baseG);
            baseL = transCoord2d('ang G2L',cOrigin,-S.ang(j)/180*pi,baseG);
            
            % Store local coordinates
            S.ft(k).xTipL(j,1)   = tipL(1);
            S.ft(k).yTipL(j,1)   = tipL(2);
            S.ft(k).xBaseL(j,1)  = baseL(1);
            S.ft(k).yBaseL(j,1)  = baseL(2);
            
            clear tipL tipG baseL
        end
        
        % Tip point in region of interest
        %tipROI(1,:) = tipL(1,:) + size(im_roi,1)/2;
        %tipROI(2,:) = tipL(2,:) + size(im_roi,2)/2;
        
        clear iH
    end
end


%% Translate arm coorinates into local FOR

% Loop thru frames with arm coordinates
for j = 1:length(A.frames)
    
    % Index of current frame in H data
    iFrameH = find(H.frames==A.frames(j),1,'first');
    
    % Index of current frame in S data
    iFrameS = find(S.frames==A.frames(j),1,'first');
    
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
    
    clear armG
end

% Find mean arm position
armL = mean(armL,3);

% Step thru frames to transform arms into global FOR
for j = 1:length(S.xCntr)
    
    % Current origin
    cOrigin = [S.xCntr(j) S.yCntr(j)];
    
    % Calcuate the global position of the base
    %armG = G2L(cOrigin,-S.ang(j)/180*pi,armL);
    armG = transCoord2d('ang L2G',cOrigin,-S.ang(j)/180*pi,armL);
    
    % Store coordinates
    S.arm(j).x = armG(:,1);
    S.arm(j).y = armG(:,2);
    
    clear cOrigin armG
end

% Store local coordinates
S.armL.x = armL(:,1);
S.armL.y = armL(:,2);

clear armL


%% Fill out base coordinates from position measurement

% Loop thru tube feet (write coordinates for base)
for j = 1:length(H.ft)
    
    % Find indicies for contact and release
    S.ft(j).iStart = find(~isnan(S.ft(j).xTip),1,'first');
    S.ft(j).iEnd   = find(~isnan(S.ft(j).xTip),1,'last');
    
    % If there's a base point
    if sum(~isnan(S.ft(j).xBase))>0
        
        % Find local position of base
        idx = find(~isnan(S.ft(j).xBase),1,'first');
        baseL = [S.ft(j).xBaseL(idx) S.ft(j).yBaseL(idx)];
        %                 baseR(1,1) = baseL(1,:) + size(im_roi,1)/2;
        %                 baseR(2,1) = baseL(2,:) + size(im_roi,2)/2;
        
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
                
                clear baseG
            end
        end
        
        % If no base point . . .
    else
        warning(['    Foot ' num2str(j) ' in does not have a base'])
    end
end


%% Add time vector

S.t = S.frames ./ v.FrameRate;


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