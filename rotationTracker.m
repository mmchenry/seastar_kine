function varargout = rotationTracker(vid_path,v,imInvert,method,varargin)
% Tracks roation in the motion of an object in a video sequence
% Rotation = tracker(vid_path,v,imInvert,'advanced rotation',roi0,Centroid,frames) 
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


%TODO: Update rotation code

% Default mask coordinates
xMask = [];
yMask = [];


%% Parse inputs

if strcmp(method,'advanced rotation') 

    r          = varargin{1};
    Centroid   = varargin{2};
    framesIn   = varargin{3};
    savePath   = varargin{4};
    
    % Over how many frames data is saved
    saveInterval = 20;
    
    imMean = [];
    numroipts = 400;
    dSample = 1;
    
     % Threshold rotation at which a reference frame is created (deg)
    rot_thresh = 15;
    
    % Starting rotation angle
    rot_ang_net = 0;
    
    frameInterval = 20;
    
    
elseif strcmp(method,'advanced rotation with mask') || ...
       strcmp(method,'advanced rotation with simple mask')

    r          = varargin{1};
    Centroid   = varargin{2};
    framesIn      = varargin{3};
    
    % Mean image
    if nargin<8
        imMean = [];
    else
        imMean = varargin{4};
    end
    
    % Mask coordinates
    if nargin>8
        xMask = varargin{5};
        yMask = varargin{6};
    end
    
    % Number of points to define roi
    numroipts = 400;
    
    if nargin>10
        dSample = varargin{7};
    else
        % Downsample images for image registration
        dSample = 0;
    end
    
    if strcmp(method,'advanced rotation with simple mask')
        tVal = varargin{8};
    end
    
    % Threshold rotation at which a reference frame is created (deg)
    rot_thresh = 15;
    
    % Starting rotation angle
    rot_ang_net = 0;

% If no match on method
else
    error('requested method not recognized');
    
end

%% Determine frames to start analysis

if ~isempty(dir(savePath))
    
    % Load existing Rotation data
    load(savePath,'Rotation');
    
    % Make frames vector start from end of last data
    iFrame = length(Rotation.rot_ang)+1;
    
    % Frame vector
    frames = framesIn(iFrame:end);

else
    iFrame = 1;
    frames = framesIn;
    
end

clear framesIn 


%% Parameter defaults

% If threshold method . . .
if strcmp(method,'advanced rotation') || ...
       strcmp(method,'advanced rotation with mask') || ...
       strcmp(method,'advanced rotation with simple mask')
    
    % Initialize image registration parameters
    [optimizer, metric]  = imregconfig('monomodal');
    optimizer.MaximumStepLength = 5e-4;
    optimizer.MaximumIterations = 1500;
    optimizer.RelaxationFactor  = 0.2;
    
    % First image
    im0 = getFrame(vid_path,v,frames(1),imInvert,'gray',imMean);
    
    % Make binary mask of tank region
    bwMask = roipoly(im0,xMask,yMask);
    
    % Eliminate outside of roi
    im(~bwMask) = 255;
    
    % Current roi
    roi = giveROI('define','circular',numroipts,r,Centroid.x(1),Centroid.y(1));
    
    % Log first roi
    roi0 = roi;
    
    % Focus on roi
    %[im0,roi_mask,roi_rect] = isolate_roi(im0,Centroid.x_pix(1),Centroid.y_pix(1),Centroid.r_pix,theta);
    [im_roi0,bw_mask] = giveROI('unstabilized',im0,roi,dSample);
   
    % Add mask
    if strcmp(method,'advanced rotation with mask')
        im_roi0 = addmask(im_roi0);
        
    elseif strcmp(method,'advanced rotation with simple mask')
        im_roi0 = addsimplemask(im_roi0,tVal);
    end
    
    % Set starting roi structure
   % roi = roi0;
end


%% Tracking object

% Initialize counter
nFrames = 1;

% Loop thru frames
for i = 1:length(frames)
    
    % Current frame
    cFrame = frames(i);
    
    % Current image
    im = getFrame(vid_path,v,cFrame,imInvert,'gray',imMean);
    
    % Eliminate outside of mask
    if ~isempty(xMask)
        
        % Make binary mask of tank region
        bwTank = roipoly(im,xMask,yMask);
        
        % Apply mask
        im(~bwTank) = 255;
        
        clear bwTank
    end
    
    if strcmp(method,'body rotation')
        
%         % Find rotation matrix
%         tform = findRot(im,im0,Centroid.x_pix(i),Centroid.y_pix(i),...
%                         Centroid.r_pix,theta,optimizer,metric);
        numroipts = length(roi0.theta);
        
        roi = giveROI('define','circular',numroipts,roi0.r,...
                      Centroid.x(iFrame),Centroid.y(iFrame));
                  
         % Give image, unstabilized (i.e. unrotated)
         [im_roi,bw_mask] = giveROI('unstabilized',im,roi,dSample);
    
         % Focus on roi
%         [im_roi,imMask,roi_rect] = isolate_roi(im,Centroid.x_pix(i),...
%             Centroid.y_pix(i),Centroid.r_pix,theta);

         % Transformation object to stablize head wrt im0
         tform = imregtform(im_roi,im_roi0,'rigid',optimizer,metric);

        % Store results
        Rotation(iFrame).tform_roi = tform;
       % Rotation(iFrame).roi_rect = roi_rect;
         
        % Clear for next iteration
        clear im_roi tform_roi imStable xC yC h imMask im_roi       
        
    % body rotation method: 
    elseif strcmp(method,'advanced rotation') || ...
           strcmp(method,'advanced rotation with mask') || ...
           strcmp(method,'advanced rotation with simple mask')
        
        % Current center points
        cX = Centroid.x(iFrame);
        cY = Centroid.y(iFrame);
        
        if ~isempty(xMask)
            % Make binary mask of tank region
            bwMask = roipoly(im,xMask,yMask);
            
            % Eliminate outside of mask
            im(~bwMask) = 255;
        end

        % Current roi
        roi = giveROI('define','circular',numroipts,r,cX,cY);

        % Give image, unstabilized (i.e. unrotated)
        [im_roi,bw_mask_roi] = giveROI('unstabilized',im,roi,dSample);
           
        % If masking . . .
        if strcmp(method,'advanced rotation with mask')          
            im_roi = addmask(im_roi);
            
        elseif strcmp(method,'advanced rotation with simple mask')
            im_roi = addsimplemask(im_roi,tVal);
        end
        
         % If first frame . . .
         if iFrame==1
             
             % tform is the identity matrix
             tform = affine2d(eye(3));
             
             % Mark keyframe
             Rotation.ref_frame(iFrame,1) = 1;            

             % Angular rotation up to this point
             Rotation.rot_ang(iFrame,1) = 0;
             
             % Current roi image
             im_roi0curr = im_roi0;
             
         % If after first frame . . .
         else
  
             % Rotate reference image to last rotated angle
             im_roi0curr = giveROI('stabilized',im0,roi0,dSample,-Rotation.rot_ang(iFrame-1));
             
             % If masking . . .
             if strcmp(method,'advanced rotation with mask')
                 im_roi0curr = addmask(im_roi0curr);
             end
             
             % Transformation object to stablize head wrt im0
             tform = imregtform(im_roi,im_roi0curr,'rigid',optimizer,metric);
             
             % Get change in angular rotation
             tformInv       = invert(tform); 
             Drot_ang       = atan2(tformInv.T(2,1),tformInv.T(1,1))*180/pi;
             
             % Zero out Drot_ang, if too big
             if abs(Drot_ang)>120                
                Drot_ang = 0;                
             end
             
             % Angular rotation up to this point
             Rotation.rot_ang(iFrame,1)  = Rotation.rot_ang(iFrame-1,1) + Drot_ang
             
             % No keyframe
             Rotation.ref_frame(iFrame,1) = 0;           
         end
    end
    
    % Update status
    disp(['tracker (' method ') : done ' num2str(iFrame) ' of ' num2str(length(frames))])   
    
    % Visualize rotation, for debugging 
    if 1
        
        %imStable = giveROI('stabilized',im,roi,dSample,tform);
        imStable = giveROI('stabilized',im,roi,dSample,Rotation.rot_ang(iFrame,1));
        theta = linspace(0,2*pi,400);
        
        includeRot = 1;
        
        % Title text
        t_txt = ['Frame ' num2str(cFrame) '/' num2str(frames(end))];
            
        % If rotation data included . . .
        if includeRot
            
            subplot(2,2,1)
            imshow(im,'InitialMag','fit')
           % brighten(-0.7)
            hold on
            %plot(Centroid.x(iFrame),Centroid.y(iFrame),'g+',xMask,yMask,'k-')
            line(roi.xPerimG,roi.yPerimG,'Color',[1 0 0 0.2],'LineWidth',1);
            hold off
            title(t_txt)
            
            subplot(2,2,2)
            imshow(imStable,'InitialMag','fit')
            hold on
            %plot(roi.xPerimL,roi.yPerimL,'k-')
            line(roi.xPerimL,roi.yPerimL,'Color',[1 0 0 0.2],'LineWidth',3);
            hold off        
            %brighten(-0.7)
            
            subplot(2,2,3)
            imshowpair(im_roi,im_roi0,'falsecolor')
            title('Unrotated reference')
            
            subplot(2,2,4)
            %imshowpair(giveROI('stabilized',im,roi,dSample,tform),im_roi0)
            imshowpair(im_roi,im_roi0curr,'falsecolor')
            title('Rotated reference')
            
%             visTrack(im,Centroid.x(iFrame),Centroid.y(iFrame),r,theta,...
%                 Rotation(iFrame).tform_roi,t_txt);
            
            % Pause briefly to render
            pause(0.001)
        end
    end   

     % Clear for next iteration
     clear im_roi tform_roi imStable xC yC h imMask im_roi t_txt

     if nFrames >= frameInterval
        save(savePath,'Rotation')
        
        nFrames = 1;
     else
         nFrames = nFrames + 1;
     end
     
     iFrame = iFrame + 1;
end


%% Define outputs
       
% Body rotation method    
if strcmp(method,'body rotation') || ...
       strcmp(method,'advanced rotation') || ...
       strcmp(method,'advanced rotation with mask') || ...
       strcmp(method,'advanced rotation with simple mask')
    
    varargout{1} = Rotation;
    
elseif strcmp(method,'visualize') 
    
    if makeVid
        varargout{1} = M;
    else
        varargout{1} = fig;
    end
end


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


