function generateArmMovies(vid_path,data_path,v,imInvert,varargin)
% Exports movies of individual arms of sea stars that can be used for
% tracking individual tube feet

%% Parse inputs

Body        = varargin{1};
iC          = varargin{2};
echoFrames  = varargin{3};


% Extract fields from Body
S        = Body.Rotation;
frames   = Body.frames;

clear Body


%% Parameter values

% Number of agular points to use
nAng = 100;


%% Set up for importing and exporting video

% Path for mean images
mPath = [data_path filesep 'mean_images'];

% Directory for arm videos
iLast     = find(vid_path==filesep,1,'last');
armDir    = [vid_path(1:(iLast-1)) filesep v.Name(1:(end-4)) '_arms'];

% Make folder, if not there
if ~isfolder(armDir)
    mkdir(armDir);
end

% Loop thru arms
%TODO: Set this up for 5 arms
for i = 1:5
    
    % File name of output videos
    outMovie{i}.fName = [v.Name(1:(end-4)) '_arm' num2str(i)];
    
    % Path for output video
    outMovie{i}.fPath = [armDir filesep outMovie{i}.fName];
    
    % Set up output video file
    outMovie{i}.v = VideoWriter([outMovie{i}.fPath '.avi'],'Motion JPEG AVI');
    outMovie{i}.v.Quality = 90;
    
    % Initialize video for writing
    open(outMovie{i}.v)
end

clear iLast armDir


%% Add arms in global FOR


% Add arm points
%Body = addArms(Body,iC);





%% Perform export for each frame

% Loop thru frames
for i = 1:length(frames)
    
    % Current frame
    cFrame = frames(i);
    
    % Get mean image for current frame
    [imRoiMean,imRoiStd] = getMeanImage(cFrame,mPath);
    
    % Read current image
    im = getFrame(vid_path,v,cFrame,imInvert,'gray',[],iC.r);
    
    % Get current roi
    roi = S.roi(i);
    
    % Stabilized image
    imStable =  giveROI('stabilized',im,roi,0,S.tform(i));
    
    % Brighten mean image
    imRoiMean = imadjust(imRoiMean,...
                         [min(imRoiMean(:))/255; max(imRoiMean(:))/255],[50/255; 1]);
    
    % Subtract mean image, adjust contrast
    imStable = adapthisteq(imsubtract(imRoiMean,imStable));
    
    % Position of arms in local FOR
    xArmsL = iC.xArms - iC.x + iC.r;
    yArmsL = iC.yArms - iC.y + iC.r;
    
    % Center of body in local FOR
    xCntrL = roi.rect(4)/2;
    yCntrL = roi.rect(3)/2;
    
    % Loop thru arms
    for j = 1:5
        
        % Masked version of stable image
        imStableM = imStable;
        
        % Radius of ROI
        r = 1.4*hypot(yArmsL(j)-yCntrL,xArmsL(j)-xCntrL);
        
        % Distance for backing off of center
        rBack = 0.1*hypot(yArmsL(j)-yCntrL,xArmsL(j)-xCntrL);
        
        % Angle of arm in local FOR
        angArm = atan2(yArmsL(j)-yCntrL,xArmsL(j)-xCntrL);
        
        % Center point for mask
        xC = xCntrL - rBack*cos(angArm);
        yC = yCntrL - rBack*sin(angArm);
        
        % Coordinates for ROI
        xROI = [xC; xC+r.*cos(linspace(angArm-pi/5,angArm+pi/5,nAng)'); xC]; 
        yROI = [yC; yC+r.*sin(linspace(angArm-pi/5,angArm+pi/5,nAng)'); yC]; 
        
        % binary ROI
        bwROI = roipoly(imStableM,xROI,yROI);
        
        % Apply mask
        imStableM(~bwROI) = 0;
        
        % Write video frame
        drawnow
        writeVideo(outMovie{j}.v,imStableM);
        
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

    if echoFrames
        disp(['  Arm movie export: Frame ' num2str(i) ' of ' num2str(length(frames))])
    end
end




%% Finish up

% Loop thru movies
for i = 1:length(outMovie)
    % Initialize video for writing
    close(outMovie{i}.v)
end



function [imRoiMean,imRoiStd] = getMeanImage(cFrame,mPath)

% Listing of mean images
a = dir([mPath filesep 'mean*']);

for i = 1:length(a)
    
    % Index of separators
    iSep = find(a(i).name=='_');
    
    % Get start frame
    frStart = str2num(a(i).name((iSep(2)+1):(iSep(3)-1)));
    
    % Get end frame
    frEnd = str2num(a(i).name((iSep(3)+1):(end-4)));
    
    if cFrame>=frStart && cFrame<=frEnd
        % Load imean image data
        load([mPath filesep a(i).name])
        
        % Define 
        imRoiMean = roiM.im;
        imRoiStd  = roiM.imStd;

        break
    end
end

% Check for definition
if ~exist('imRoiMean','var')
    error(['No match for cFrame = ' num2str(cFrame)])
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
     
     
%     % Loop thru feet in local FOR
%     if ~isempty(B2(i).L)
%         for j = 1:length(B2(i).L)
%             % Angular position of foot
%             angFt = atan2(B2(i).L(j).Centroid(2) - iC.r,...
%                           B2(i).L(j).Centroid(1) - iC.r);
%                       
%             % Diff in ang wrt arms
%             for k = 1:length(angArms) 
%                 dAng(k,1) =  abs(angdiff(angFt,angArms0(k)));
%             end
%    
%             % Store arm number
%             B2(i).L(j).armNum = find(dAng==min(dAng),1);
%             
%             % Store angular position
%             B2(i).L(j).ang = angFt;
%         end
%     end
%     
% 
%     % Loop thru feet in global FOR
%     if ~isempty(B2(i).G)
%         for j = 1:length(B2(i).G)
%             % Angular position of foot
%             angFt = atan2(B2(i).G(j).Centroid(2) - yCntr,...
%                           B2(i).G(j).Centroid(1) - xCntr);
%             
%             % Diff in ang wrt arms          
%             for k = 1:length(angArms)
%                 dAng(k,1) =  abs(angdiff(angFt,angArms(k)));
%             end
%             
%             % Store arm number
%             B2(i).G(j).armNum = find(dAng==min(dAng),1);
%             
%             % Store angular position
%             B2(i).G(j).ang = angFt;
%         end
%     end
    
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



