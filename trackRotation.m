function trackRotation(vid_path,v,currDataPath,method,varargin)
% Tracks the motion of an object in a video sequence
%   vid_path - path to video file or image sequence
%   v - structure of info about video (generated by defineVidObject)
%   imInvert - logical that indicates whether to invert the images
%   method - type of tracking approach ('threshold translation' or 'body
%   rotation')
%
% Rotation = tracker(vid_path,v,imInvert,'body rotation',roi0,Centroid,frames,imMean) 
% uses advanced algoriothm to return just centroid cooridinates
%   Rotation - Structure with rotation data
%
% Rotation = tracker(vid_path,v,imInvert,'advanced rotation',roi0,Centroid,frames,imMean,xMask,yMask) 
% uses advanced algoriothm to return just centroid cooridinates
%   Rotation - Structure with rotation data

%
%  imInvert - choose 0 for light on dark field, 1 for dark on light field
%  visTracking - Logical to visualize the tracking
%  frames - listing of frame numbers to analyze
%  Centroid - structure with fields x_pix, y_pix, & r_pix for roi coord & radius
%   
%
% Developed by McHenryLab at UC Irvine


%% General parameters

% Number of analyzed frames at which data is save
saveInterval = 10;

% Maximum size of an image dimension (for downsampling)
maxSize = 250;

currDataFile = 'Body.mat';

% Predator & prey colors
clrs{2} = [41 171 226]./255;
clrs{1} = [241 90 36]./255;

% Level to brighten the movie
bLevel = 0.5;

if strcmp(method,'advanced') 
    
    % Initialize image registration parameters
    [optimizer, metric]  = imregconfig('monomodal');
     optimizer.MaximumStepLength = 5e-4;
    optimizer.MaximumIterations = 1500;
    optimizer.RelaxationFactor  = 0.2;

end


%% Load/create data files

if isfile([currDataPath filesep 'Body.mat'])
    load([currDataPath filesep 'Body.mat'])
else
    % Load Centroid
    load([currDataPath filesep 'Centroid.mat'])
    
    % Transfer data to Body structure
    Body.frames = Centroid.frames;
    Body.x      = Centroid.x;
    Body.y      = Centroid.y;
    Body.y_flip = Centroid.y_flip;
    
    clear Centroid
end

% Load iC
load([currDataPath filesep 'Initial conditions'])


%% Parse inputs

% Color mode 
clrMode = 'gray';
%clrMode = 'rgb';

% Transfer Tank Mask
r = iC.r;

% Number of points to define roi
numroipts = 400;

% Downsample
dSample = varargin{1};


if length(varargin)>3
    visSteps = varargin{2};
    iminvert = varargin{3};
else
    visSteps = 0;
    iminvert = 0;
end

clr = clrs{1};


%% Determine frames to analyze


% Analyze all possible frames
ana.auto = ones(length(Body.frames),1);

% Set current frame to first for reference image
cFrame = Body.frames(1);

iFrame = find(Body.frames==cFrame,1,'first');
x0     = Body.x(iFrame);
y0     = Body.y(iFrame);
r      = iC.r;

% Enlarge r a bit
%r = 1.25*r;


% Downsmampling used only for image registartion
%dSample = 0;

clear idx


%% Define reference image
    
% First image
% im0 = getFrame(vid_path,v,1,iminvert,clrMode);
im0 = getFrame(vid_path,v,cFrame,iminvert,clrMode);

% Current roi
roi = giveROI('define','rectangular',numroipts,r,x0,y0);

% Log first roi
roi0 = roi;

% Focus on roi
[im_roi0,bw_mask] = giveROI('unstabilized',im0,roi,dSample);

% High-res version of reference image
[im_roi0HR,bw_maskHR] = giveROI('unstabilized',im0,roi,0);

% Adjust image
%im_roi0 = im_modify(im_roi0);
%im_roi0HR = im_modify(im_roi0HR);

% Counter for saving data
nSave = 1;

% If re-running analysis . . . 
if ~isfield(Body,'Rotation')
    
    % Set index at the start
    iFrame = 1;   
    
% If resuming analysis . . .
else
    % Jump to next unanalyzed frame
    iFrame = length(Body.Rotation.tform)+1;
end

if iFrame~=1
    disp(' '); 
    disp('========== Resuming analysis of sequence ===============');
    disp(' ')
end

%iFrame = find(Body.frames,1,'first');


%% Tracking object

% Loop thru frames
while true
        
    % Stop, if beyond duration
    if iFrame > length(Body.frames)
        break
    end
    
    % Current image
    cFrame = Body.frames(iFrame);
    
    im = getFrame(vid_path,v,cFrame,iminvert,clrMode);
    
    % Center point
    cX = Body.x(iFrame);
    cY = Body.y(iFrame);
    
    % Current roi
    roi = giveROI('define','rectangular',numroipts,r,cX,cY);
    
    % Give image, unstabilized (i.e. unrotated)
    [im_roi,bw_mask_roi] = giveROI('unstabilized',im,roi,dSample);
    
    % Modify image
   % im_roi = im_modify(im_roi);
    
    % If masking . . .
    if strcmp(method,'advanced rotation with mask')
        im_roi = addmask(im_roi);
        
    elseif strcmp(method,'advanced rotation with simple mask')
        im_roi = addsimplemask(im_roi,tVal);
    end
    
    % If first frame...
    if iFrame == 1
        
        % tform is the identity matrix
        tform = affine2d(eye(3));
        
        % Mark keyframe
        Body.Rotation.ref_frame(iFrame,1) = 1;
        
        % Angular rotation up to this point
       % Body.Rotation.rot_ang(iFrame,1) = 0;
        
        % Store tform
        Body.Rotation.tform(iFrame) = tform;
        Body.Rotation.roi(iFrame) = roi;
        
        % Current roi image
        im_roi0curr = im_roi0;
        
        %iFrame = 0;
        
    % If after first frame . . .
    else
        
        % last transformation and angle
        tform_last = Body.Rotation.tform(iFrame-1);
        
        % Round angle
        anglD = round(atan2(tform_last.T(1,2),tform_last.T(1,1)) * 180/pi,6);
    
        % Insert new values into the T-matrix
        tform_last.T(1:2,1:2) = [cosd(anglD) sind(anglD) ;-sind(anglD) cosd(anglD)];
        
        % Round the positional data too
        tform_last.T(3,1:2) = round(tform_last.T(3,1:2),6);
        
        % Transformation object for displacement since last frame
        tform = imregtform(im_roi,im_roi0,'rigid',optimizer,metric,...
            'InitialTransformation',tform_last);
            
        clear im_roiHR bw_mask_roiHR tform_last ang_last
   
        % Store transformation matrix
        Body.Rotation.tform(iFrame) = tform;
        Body.Rotation.roi(iFrame)   = roi;
        
        %clear displ disp_net tformInv Drot_ang tform netRot_ang T Drot_rad
    end
    
    % Update status
    disp(['trackRotation (' method ') : done frame ' num2str(Body.frames(iFrame)) ...
          ' of ' num2str(max(Body.frames))])
    
    %Visualize rotation, for debugging
    if visSteps

        % Use roi from current frame
        %roiCurr = giveROI('define','rectangular',numroipts,r,Body.x(iFrame),Body.y(iFrame));
        
        %imStable = giveROI('stabilized',im,roi,dSample,tform);
        imStable =  giveROI('stabilized',im,Body.Rotation.roi(iFrame),...
                            dSample,Body.Rotation.tform(iFrame));
        
        % Modify image
        %imStable = im_modify(imStable);
        
        warning off
        
        subplot(2,2,1)
        imshow(im_roi,'InitialMag','fit')
        title(['Current: frame ' num2str(cFrame)])
        hold on
        %plot(roi.xPerimL,roi.yPerimL,'k-')
        line(roi.xPerimL,roi.yPerimL,'Color',[1 0 0 0.2],'LineWidth',3);
        hold off
        %brighten(-0.7)
        
        subplot(2,2,2)
        imshow(imStable,'InitialMag','fit')
        title('Stabilized roi')
        hold on
        %plot(roi.xPerimL,roi.yPerimL,'k-')
        line(roi.xPerimL,roi.yPerimL,'Color',[1 0 0 0.2],'LineWidth',3);
        hold off
        %brighten(-0.7)
        
        subplot(2,2,3)
        imshow(im_roi0,'InitialMag','fit')
        title('Reference roi')
        hold on
        line(roi.xPerimL,roi.yPerimL,'Color',[1 0 0 0.2],'LineWidth',3);
        hold off

        if iFrame>1
            subplot(2,2,4)
            imshowpair(im_roi0,imStable)
            title('Image comparison')
        end
        
        warning on
    end
    
    
    % Pause briefly to render
    pause(0.001)
    
    
    % Visualize centroid, for debugging
    if 0
        
        % Title text
        t_txt = ['Frame ' num2str(cFrame) '/' num2str(frames(end))];
        
        imshow(im,'InitialMag','fit')
        if 1
            brighten(-0.8)
        end
        % brighten(-0.7)
        hold on
        plot(cX,cY,'g+')
        title(t_txt)
        hold off
        
        % Pause briefly to render
        pause(0.001)
    end
    
    %%  % Clear for next iteration
    clear im_roi tform_roi imStable xC yC h imMask im_roi t_txt
    
    if nSave > saveInterval
        
        save([currDataPath filesep currDataFile],'Body')
        nSave = 1;
    else
        nSave = nSave + 1;
        
    end
    
    % Advance index
    iFrame = iFrame + 1;
    
    clear bwOut bw_mask_roi currX currY currAng  roi
end

% Save data
save([currDataPath filesep currDataFile],'Body')


%% Define outputs

% Threshold method
if strcmp(method,'threshold translation') || ...
        strcmp(method,'threshold roi') || ...
        strcmp(method,'thresh trans advanced') || ...
        strcmp(method,'thresh trans advanced with mask')
    
        varargout{1} = Centroid;
            
% Body rotation method    
elseif strcmp(method,'body rotation') || ...
       strcmp(method,'advanced rotation') || ...
       strcmp(method,'advanced rotation with mask') || ...
       strcmp(method,'advanced rotation with simple mask')
    
%     varargout{1} = Rotation;
    
elseif strcmp(method,'visualize') 
    
    if makeVid
        varargout{1} = M;
    else
        varargout{1} = fig;
    end
end

function im = im_modify(im)

% Create blank image
im0 = im*0;

% Position of image center
xC = size(im,1)/2;
yC = size(im,2)/2;

% Coordinates of each pixel
[X,Y] = meshgrid(1:size(im,1),1:size(im,2));

% Distance from center, scaled to 255
pVal = hypot(X-size(im,1)/2,Y-size(im,2)/2)^2;
pVal = uint8(round(pVal./max(pVal)*255));

% Adjust contrast of input image
imRef = imadjust(im);

% Add two images 
im = imadd(pVal,imRef);


function im = addmask(im)
% Mask with distance map (trims fins)

% Enhance contrast
im1 = imadjust(im);

% Binary image
bw = im2bw(im1,graythresh(im1));

% Distance map
bwD = bwdist(bw);

% Max distance
maxDist = max(bwD(:));

% Threshold distance
threshDist = maxDist/3;

% Refine binary as above-threshold images
bw = bwD>threshDist;

% Dilate, white out outside
se = strel('disk',3,4);
bw = imdilate(bw,se);
im(~bw) = 255;

function im = addsimplemask(im,tVal)
% Adds a regular mask to image
%im1 = imadjust(im);
im1 = im;
bw = im2bw(im1,tVal);
se = strel('disk',3,4);
bw = imdilate(~bw,se);
im(~bw) = 255;
