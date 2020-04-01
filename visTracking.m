function visTracking(vid_path,v,data_path,imInvert,opType,varargin)


% Downsample
dSample = 0;



%% Handle inputs

if strcmp(opType,'Basic')
    
    Body     = varargin{1}; 
    imVis    = varargin{2};
    F        = varargin{3};
    iFrames  = varargin{4};
    iC       = varargin{5};
    nFrames  = varargin{6};

    % Frame numebrs to analyze
    frames   = Body.frames(iFrames)';
    
    % Approximate interval between frames
    apporxInt = round(length(frames)/nFrames);
    
    % Frame numbers to visualize
    visFrames = round(linspace(min(frames)+apporxInt,max(frames)-apporxInt,nFrames));
    
    % Number of frames for streak image
    nStreak = 50;
    
    % Get coord transformation
    S = Body.Rotation;
end


%% Extract coordinates fir each frame

for i = 1:length(frames)
    
    % Number of feet in frame
    n = 0;
    
    % Index for body data
    iMatchBody = Body.frames==frames(i);
    
    % Arm coordinates
    xArm(i,:)   = Body.xArmG(iMatchBody,:);
    yArm(i,:)   = Body.yArmG(iMatchBody,:);
    
    % Body center
    xCntr(i,:)   = Body.xCntr(iMatchBody);
    yCntr(i,:)   = Body.yCntr(iMatchBody);
    
    % Loop thru feet
    for j = 1:length(F)
        
        % Index of current frame to F data
        iMatch = F(j).frames==frames(i);
        
        % Check for repeats
        if sum(iMatch)>1
            error('Multiple instances of same frame for a foot');
        end
        
        % If there is a match, store coordinates
        if max(iMatch)
            % Advance index
            n = n + 1;
            
            % Store coordinates and colors
            x{i}(n,1)      = F(j).xG(iMatch);
            y{i}(n,1)      = F(j).yG(iMatch);
            clr{i}(n,:)    = F(j).clr(1,:);
        end
    end
    
    % Add nans, if no matching feet
    if n==0
        x{i}(1)      = nan;
        y{i}(1)      = nan;
        clr{i}(1,:)  = F(j).clr(1,:);
    end
end


%% Loop thru frames

% im(:,:,1) = double(getFrame(vid_path,v,visFrames(1),imInvert,'gray'));

% Loop thru frames
for i = 1:length(visFrames)
    
    % Current frames
    cFrames = [visFrames(i)-round(nStreak/2):5:(visFrames(i)+round(nStreak/2)-1)];
    
    im = motionImage(vid_path,v,'mean','none',...
                  imInvert,Body,dSample,cFrames,iC);
              
              
    % Index for current frame in the data
    iFrame = find(visFrames(i)==Body.frames,1','first');     
    
%     imROI = giveROI('unstabilized',im,S.roi(iFrame),dSample);
    
%     im(:,:,2) = double(getFrame(vid_path,v,visFrames(i),imInvert,'gray'));
    
    
    
%     imM = uint8(mean(im,3));
    
    imshow(im)
    
%     im = getFrame(vid_path,v,cFrames(i),imInvert,'gray');
    
    
    % Read images to make streak image
%     im = imStreak(vid_path,v,visFrames(i),imInvert,nStreak);
    
    hold on
    
    for j = 1:length(x{iFrame})
        %                 h = scatter(x{i}(j),y{i}(j),...
        %                     'MarkerEdgeColor',clr{i}(j,:),'SizeData',200,...
        %                     'MarkerEdgeAlpha',0.8);
        h = scatter(x{iFrame}(j),y{iFrame}(j),...
            'MarkerEdgeColor','r','SizeData',150,...
            'MarkerEdgeAlpha',0.5,'LineWidth',1);
    end
    
%     im(:,:,1) = im(:,:,2);

ttt=4;
    
end




function im = imStreak(vid_path,v,frNum,imInvert,nStreak)

% Current frames
cFrames = [frNum-round(nStreak/2) :5:  (frNum+round(nStreak/2)-1)];

% for i = 1:length(cFrames)
%     % Current whole frame
%     imC(:,:,i) = getFrame(vid_path,v,cFrames(i),imInvert,'gray');
% end


ttt=3;