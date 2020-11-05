function generateMaskMovie(vid_path,data_path,v,imInvert,varargin)
% Exports movies of individual arms of sea stars that can be used for
% tracking individual tube feet

%% Parse inputs

Body        = varargin{1};
iC          = varargin{2};
echoFrames  = varargin{3};

% Extract fields from Body
S        = Body.Rotation;
frames   = Body.frames;
xArmG    = Body.xArmG;
yArmG    = Body.yArmG;

% Index of nans
iNan = isnan(xArmG(:,1));

% If there are more than 2 blocks of nans, throw error
if sum(abs(diff(iNan)))>2
    error('More than 2 blocks of nans in the data');
end

% Remove frames with nan data
frames = frames(~iNan);

% Range of arm coordinates
xMin = min(xArmG(:));
xMax = max(xArmG(:));
yMin = min(yArmG(:));
yMax = max(yArmG(:));

% Border for crop
b_len = nanmax(range(xArmG,2))*.25;

% ROI for cropping
crop_rect = [xMin-b_len yMin-b_len xMax-xMin+2*b_len yMax-yMin+2*b_len]; 


clear Body iNan xMin xMax yMin yMax


%% Parameter values

% Number of agular points to use
nAng = 100;


%% Set up for importing and exporting video

% Path for mean images
mPath = [data_path filesep 'mean_images'];

% Directory for arm videos
iLast     = find(vid_path==filesep,1,'last');
armDir    = [vid_path(1:(iLast-1)) filesep v.Name(1:(end-4)) '_masked'];

% Make folder, if not there
if ~isfolder(armDir)
    mkdir(armDir);
end

% Loop thru arms
for i = 1:5
    
    % File name of output videos
    outMovie{i}.fName = [v.Name(1:(end-4)) '_masked' num2str(i)];
    
    % Path for output video
    outMovie{i}.fPath = [armDir filesep outMovie{i}.fName];
    
    % Set up output video file
    %outMovie{i}.v = VideoWriter([outMovie{i}.fPath '.avi'],'Motion JPEG AVI');
    outMovie{i}.v = VideoWriter([outMovie{i}.fPath '.mp4'],'MPEG-4');
    outMovie{i}.v.Quality = 90;
    
    % Initialize video for writing
    open(outMovie{i}.v)
end

clear iLast armDir


%% Add arms in global FOR

% Add arm points
%Body = addArms(Body,iC);


%% Perform export for each frame

% Get mean image for current frame
imM = getMeanImage(mPath);

% Loop thru frames
for i = 1:length(frames)
    
    % Current frame
    cFrame = frames(i);
    
    % Read current image
    im = getFrame(vid_path,v,cFrame,imInvert,'gray',[],iC.r);
    
    % Get current roi
    roi = S.roi(i);   
    
    % Center coordinates
    xCntr = S.roi(cFrame).xCntr;
    yCntr = S.roi(cFrame).yCntr;
    
    % Get arm coordinates
    xArm = xArmG(cFrame,:);
    yArm = yArmG(cFrame,:);
    
    % Masked version of stable image
    %         imStableM = imStable;
    
    % Radius of ROI
    r = 1.5.*hypot(yArm-yCntr,xArm-xCntr);
    
    % Distance for backing off of center
    rBack = 0.2.*hypot(yArm-yCntr,xArm-xCntr);
    
    % Angle of arm
    angArm = atan2(yArm-yCntr,xArm-xCntr);

    % Loop thru arms
    for j = 1:5
        
        % Current image
        imC = im;
        
        % Center point for mask
        xC = xCntr - rBack(j).*cos(angArm(j));
        yC = yCntr - rBack(j).*sin(angArm(j));

        % Coordinates for ROI
        xROI = [xC; xC+r(j).*cos(linspace(angArm(j)-pi/5,angArm(j)+pi/5,nAng)'); xC]; 
        yROI = [yC; yC+r(j).*sin(linspace(angArm(j)-pi/5,angArm(j)+pi/5,nAng)'); yC]; 
        
        % binary ROI
        bwROI = roipoly(im,xROI,yROI);
        
        % Apply mask
        imC(~bwROI) = 255;
        
        % Crop image
        imC = imcomplement(adapthisteq(imcrop(imC,crop_rect)));
        
        % Write video frame
        drawnow
        writeVideo(outMovie{j}.v,imC);
        
        % Visual check
        if 0
            imshow(imStableM,'InitialMag','fit');
            hold on
            plot(xArmsL,yArmsL,'ro',xCntrL,yCntrL,'r+')
            plot(xROI,yROI,'g-')
            hold off
            title(['Frame ' num2str(i) ', arm ' num2str(j)])
            
            pause(0.001);
        end
    end
    
    % Store mask data
    maskData.frame(i,1) = cFrame;
    

    if echoFrames
        disp(['  Masked movie export: Frame ' num2str(i) ' of ' num2str(length(frames))])
    end
end


maskData.rect = crop_rect;



%% Finish up

% Loop thru movies
for i = 1:length(outMovie)
    % Initialize video for writing
    close(outMovie{i}.v)
end

% Save mask data
save([data_path filesep 'maskData.mat'],'maskData')



function imM = getMeanImage(mPath)

% Listing of mean images
a = dir([mPath filesep 'mean*']);

for i = 1:length(a)
    
    % Index of separators
    iSep = find(a(i).name=='_');
    
    % Get start frame
    imM.frStart(i) = str2num(a(i).name((iSep(2)+1):(iSep(3)-1)));
    
    % Get end frame
    imM.frEnd(i) = str2num(a(i).name((iSep(3)+1):(end-4)));
    
    % Load imean image data
    load([mPath filesep a(i).name])
    
    % Store images
    imM.mean{i} = roiM.im;
    imM.std{i}  = roiM.imStd;
end




function Body = addArms(Body,iC)

% Set up placeholders
Body.xArmG = nan(length(Body.frames),5);
Body.yArmG = nan(length(Body.frames),5);

% Position of arms in local FOR
xArmsL = iC.xArms - iC.x + iC.r;
yArmsL = iC.yArms - iC.y + iC.r;
    
% Angular positions of arms
angArms0 = atan2(iC.yArms - iC.y,iC.xArms - iC.x);

% Loop thru time
for i = 1:length(Body.frames)
    
    % Current roi and tform
    roi   = Body.Rotation.roi(i);
    tform = Body.Rotation.tform(i);
    
    % Current body center
    xCntr = Body.x(i);
    yCntr = Body.y(i);
    
%     % Current offset
%     if isempty(B2(i).cOffset) || isnan(B2(i).cOffset(1))
%         cOffset = [0 0];
%     else
%         cOffset = B2(i).cOffset;
%     end
    
    % Rotate arm points from starting orientation to current frame
    [xArmsG,yArmsG] = transformPointsInverse(tform,xArmsL,yArmsL);
    
    % Translate by center of roi
    xArmsG = xArmsG - roi.r  + xCntr;
    yArmsG = yArmsG - roi.r  + yCntr;
    
    % Store result
    Body.xArmG(i,:) = xArmsG;
    Body.yArmG(i,:) = yArmsG;
    
    % Visual to test the transformation
    if 0
        
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
    end
end



