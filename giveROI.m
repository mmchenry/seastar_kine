
function varargout = giveROI(action,varargin)
% Returns features for a region-of-interest (ROI)
% giveROI(action)
%  action - particular operation relating to a ROI ('define')
%
% roi = giveROI('define',shape)
% Defines a region of interest
%   roi - structure of roi properties 
%   shape - string indicating the shape of the roi ('circular')
%
% roi = giveROI('define','circular',num_pts,r)
%   num_pts - number of points for bounding the roi
%   r - radius of circular roi
%   Centroid - structure definining the center coordinates for the roi
%
% roi = giveROI('define','circular',num_pts,r,xC,yC)
%   xC - roi x-center
%   yC - roi y-center
%
% im_roi = giveROI('unstabilized',im,roi,dSample)
% returns roi image, without rotation correction
%   im - full frame image
%   roi - structure with stats on roi
%   dSample - logical that indicates whether to downsample the roi image
%
% im_roi = giveROI('stabilized',im,roi,dSample,tform,imMean)
% returns roi image, WITH rotation correction
%   tform - rotation matrix for current frame 
%   imMean - (optional) mean image, used for image subtraction
%
%[im_roi,bw_mask,roi_rect,bw_roi_mask,imStable] = giveROI(...
%   im_roi - image within the roi (downsampled, if dSampling==1)
%   bw_mask - binary image for the roi with dimensions of im
%   roi_rect - [x y width height] of the roi bounding box
%   bw_roi_mask - binary of mask in roi
%   imStable - image stablized against rotation
%
% Developed by McHenryLab at UC Irvine


%% Translate inputs 

if strcmp(action,'define')
    
    roi.shape       = varargin{1};
    num_pts         = varargin{2};
    roi.r           = varargin{3};
    
    % Define center points
    if nargin>4
        roi.xCntr = varargin{4};
        roi.yCntr = varargin{5};
    else
        roi.xCntr = 0;
        roi.yCntr = 0;
    end
    
    
elseif strcmp(action,'unstabilized')
    
    im      = varargin{1};
    roi     = varargin{2};
    dSample = varargin{3};
    
    if nargin > 4
        imMean  = varargin{4};
    else
        imMean  = [];
    end
        
elseif strcmp(action,'stabilized')
    
    im      = varargin{1};
    roi     = varargin{2};
    dSample = varargin{3};
    tform   = varargin{4};
    
    if nargin > 5
        imMean  = varargin{5};
    else
        imMean  = [];
    end
    
else
    error('Do not recognize action')
end


%% Parameters

if strcmp(action,'define')
    
    if strcmp(roi.shape,'circular')
        % Angular coordinates for circular roi
        roi.theta = linspace(0,2*pi,num_pts);
        
        % Circular coordinates for new roi in image FOR
        roi.xPerimG  = round(floor(roi.r-1).*cos(roi.theta) + roi.xCntr);
        roi.yPerimG  = round(floor(roi.r-1).*sin(roi.theta) + roi.yCntr);
        
        % Circular coordinates for new roi in roi FOR
        roi.xPerimL    = round(floor(roi.r-1).*cos(roi.theta) + roi.r);
        roi.yPerimL    = round(floor(roi.r-1).*sin(roi.theta) + roi.r);
        
        % Bounding rectange for the roi
        roi.rect = [round(roi.xCntr-roi.r) round(roi.yCntr-roi.r) ceil(roi.r*2) ceil(roi.r*2)];
        
    elseif strcmp(roi.shape,'rectangular')
        % No need for angular coordinates
        roi.theta = nan;
        
        % Bounding rectange for the roi
        roi.rect = [round(roi.xCntr-roi.r) round(roi.yCntr-roi.r) ...
                    ceil(roi.r*2) ceil(roi.r*2)];
        
        % roi in image FOR
        roi.xPerimG  = [roi.rect(1) roi.rect(1)+roi.rect(3) ...
                        roi.rect(1)+roi.rect(3) roi.rect(1) roi.rect(1)];
        roi.yPerimG  = [roi.rect(2) roi.rect(2) roi.rect(2)+roi.rect(4) ...
                        roi.rect(2)+roi.rect(4) roi.rect(2)];
                    
        % roi in local FOR
        roi.xPerimL    = [1 roi.rect(3) roi.rect(3) 1 1];
        roi.yPerimL    = [1 1 roi.rect(4) roi.rect(4) 1];
        
    else
        error(['Do not recognize ' roi.shape]);
    end
    
else
    
    % Maximum size of an image dimension (for downsampling)
    maxSize = 350;
    
end


%% Create images

if ~strcmp(action,'define')
    
    % Pad image to the left, if needed
    if roi.rect(1)<1
        n  = abs(roi.rect(1))+1;
        im = [ones(size(im,1),n) im];
        roi.rect(1) = 1;
    end
     
    % Pad image to the top, if needed
    if roi.rect(2)<1
        n           = abs(roi.rect(2))+1;
        
        if size(im,3)>1
            imT = im;
            clear im
            im(:,:,1)   = [ones(n,size(imT,2)); imT(:,:,1)];
            im(:,:,2)   = [ones(n,size(imT,2)); imT(:,:,2)];
            im(:,:,3)   = [ones(n,size(imT,2)); imT(:,:,3)];
            clear imT
        else
            im   = [ones(n,size(im,2)); im];
        end
        
        roi.rect(2) = 1;
    end
        
    % Pad image to the right, if needed
    if (roi.rect(1)+roi.rect(3)+1)>=size(im,2)
        
        n = ceil(roi.rect(1)+roi.rect(3)+2 - size(im,2));
        
        if size(im,3)>1
            imT = im;
            clear im
            im(:,:,1) = [imT(:,:,1) ones(size(imT,1),n)];
            im(:,:,2) = [imT(:,:,2) ones(size(imT,1),n)];
            im(:,:,3) = [imT(:,:,3) ones(size(imT,1),n)];
            clear imT
        else
            im = [im ones(size(im,1),n)];
        end
    end
    
    % Pad image to the bottom, if needed
    if (roi.rect(2)+roi.rect(4)+1)>=size(im,1)
        
        n = ceil(roi.rect(2)+roi.rect(4)+2 - size(im,1));
        
        if size(im,3)>1
            imT = im;
            clear im
            im(:,:,1) = [imT(:,:,1); ones(n, size(imT,2))]; 
            im(:,:,2) = [imT(:,:,2); ones(n, size(imT,2))]; 
            im(:,:,3) = [imT(:,:,3); ones(n, size(imT,2))]; 
            clear imT
        else
            im = [im; 255.*ones(n, size(im,2))];     
        end
    end
    
   clear n
        
    % Binary mask in image im    
    bw_mask = roipoly(size(im,1),size(im,2),...
        [roi.rect(1) roi.rect(1)+roi.rect(3)+1 ...
         roi.rect(1)+roi.rect(3)+1 roi.rect(1) ...
         roi.rect(1)],...
        [roi.rect(2) roi.rect(2) ...
         roi.rect(2)+roi.rect(4)+1 roi.rect(2)+roi.rect(4)+1 ...
         roi.rect(2)]);
    
    % Crop image
    im_roi  = imcrop(im,roi.rect);
    
    % Check that roi is square
    if size(im_roi,1)~=size(im_roi,2)
        error('Image not square and needs to be');
    end
    
    if ~isempty(imMean)
        im_roi = imadjust(imcomplement(imsubtract(imMean,im_roi)));
    end
    
    % White out pixels outside of circular roi
    %bw_roi_mask = roipoly(size(im_roi,1),size(im_roi,2),yCirc_roi,xCirc_roi);
    bw_roi_mask = roipoly(size(im_roi,1),size(im_roi,2),roi.yPerimL,roi.xPerimL);
    
    % White out area around roi
    if size(im_roi,3)==3
        tmp1 = im_roi(:,:,1);
        tmp1(~bw_roi_mask) = 255;
        
        tmp2 = im_roi(:,:,2);
        tmp2(~bw_roi_mask) = 255;
        
        tmp3 = im_roi(:,:,3);
        tmp3(~bw_roi_mask) = 255;
        
        im_roi(:,:,1) = tmp1;
        im_roi(:,:,2) = tmp2;
        im_roi(:,:,3) = tmp3;
        
        clear tmp1 tmp2 tmp3
    else
        im_roi(~bw_roi_mask) = 255;
    end
    
    % Reduce size (helps speed up image registration)
    if dSample && (length(im)>maxSize)
        
        % Factor by which to resize
        imFactor = maxSize/length(im);
        
        % Downsample image
        im_roi = imresize(im_roi,imFactor);
        
        bw_roi_mask = imresize(bw_roi_mask,imFactor);
        % White out pixels outside of circular roi
        %bw_roi = roipoly(size(im_roi,1),size(im_roi,2),round(yC*imFactor),round(xC*imFactor));
    end

    
    if strcmp(action,'stabilized')
        
        % Coordinate system for im_roi
    R = imref2d(size(im_roi));
    
            imStable = imwarp(im_roi,tform,'OutputView',R,...
                      'FillValues',255,'SmoothEdges',true);

        % Deliver imStable
        im_roi = imStable;
        
%         if length(tform)>1 && isfield(tform,'T') && isempty(tform.T)
%             error('You need to provide tform if you want imStable');
%         end
%         
%         % If tform exists, calculate angle . . .
%         if isprop(tform,'T')        
%             % Get angular rotation from tform
%             tformInv = invert(tform);
%             rot_ang  = atan2(tformInv.T(2,1),tformInv.T(1,1))*180/pi;
%             
%         % Rotation angle passed as tform . . .
%         elseif length(tform)==1
%             rot_ang = tform;
% 
%         end
%         
%         imStable = imrotate(im_roi,-rot_ang,'bilinear','crop');
%         
% %         % Coordinate system for im_roi
% %         R = imref2d(size(im_roi));
% %         
% %         % Adjust WorldLimits to restrict transformation to just rotation
% %         % around center
% %         R.XWorldLimits = R.XWorldLimits-mean(R.XWorldLimits);
% %         R.YWorldLimits = R.YWorldLimits-mean(R.YWorldLimits);
% % %         R.XIntrinsicLimits = R.XWorldLimits;
% % %         R.YIntrinsicLimits = R.YWorldLimits;
% %         
% %         
% %         % Stablize image
% %         imStable = imwarp(im_roi,R,tform,'OutputView',R,...
% %             'FillValues',255,'SmoothEdges',true);
% 
%         if size(imStable,3)==3
%             
%             tmp1 = imStable(:,:,1);
%             tmp1(~bw_roi_mask) = 255;
%             
%             tmp2 = imStable(:,:,2);
%             tmp2(~bw_roi_mask) = 255;
%             
%             tmp3 = imStable(:,:,3);
%             tmp3(~bw_roi_mask) = 255;
%             
%             imStable(:,:,1) = tmp1;
%             imStable(:,:,2) = tmp2;
%             imStable(:,:,3) = tmp3;
%             
%             clear tmp1 tmp2 tmp3
%         else
% 
%         % White out beyond roi
%         imStable(~bw_roi_mask) = 255;
        
%       end
        
        % Deliver imStable
        im_roi = imStable;
    end
end


%% Define output

% If defining roi . . .
if strcmp(action,'define')
    
    varargout{1} = roi;
    
% Otherwise . . .
else
    
    varargout{1} = im_roi;
    
    if nargout>1
        varargout{2} = bw_mask;
        
        if nargout>2
            varargout{3} = bw_roi_mask;
            
        end
    end
end


