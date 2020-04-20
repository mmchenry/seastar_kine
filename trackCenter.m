function varargout = tracker(vid_path,v,imInvert,method,varargin)
% Tracks the motion of an object in a video sequence
%
% Developed by McHenryLab at UC Irvine


%TODO: Update rotation code

% Default mask coordinates
xMask = [];
yMask = [];


%% Parse inputs

if strcmp(method,'threshold translation')
    
    roi0   = varargin{1};
    iC     = varargin{2};
    
    % Frame numbers to analyze
    if nargin<7 || isempty(varargin{3})  
        frames = [v.UserData.FirstFrame:v.UserData.LastFrame]';
    else
        frames = varargin{3};
    end
    
    % Mean image
    if nargin<8
        imMean = [];
    else
        imMean = varargin{4};
    end
    
    % Mask coordinates image
    if nargin>8
        xMask = varargin{5};
        yMask = varargin{6};
    end
    
    erodedialate = 1;
    
    
elseif strcmp(method,'thresh trans advanced') || ...
        strcmp(method,'thresh trans advanced with mask')
    
    iC      = varargin{1};
    
    % Frame numbers to analyze
    if nargin<6 || isempty(varargin{2})  
        frames = [v.UserData.FirstFrame:v.UserData.LastFrame]';
    else
        frames = varargin{2};
    end
    
    % Mean image
    if nargin<7
        imMean = [];
    else
        imMean = varargin{3};
    end
    
    % Mask coordinates image
    xMask = varargin{4};
    yMask = varargin{5};
    
    if strcmp(method,'thresh trans advanced with mask')
       Centroid_pd  = varargin{6};
       iCPd         = varargin{7};
    end
    
    numroipts = 400;
    
elseif strcmp(method,'threshold roi')

    roi0   = varargin{1};
    tVal   = varargin{2};
    
    % Frame numbers to analyze
    if nargin<7 || isempty(varargin{3})  
        frames = [v.UserData.FirstFrame:v.UserData.LastFrame]';
    else
        frames = varargin{3};
    end
    
    % Mean image
    if nargin<8
        imMean = [];
    else
        imMean = varargin{4};
    end

% If no match on method
else
    error('requested method not recognized');
    
end


%% Parameter defaults

% If threshold method . . .
if strcmp(method,'threshold translation') || ...
        strcmp(method,'threshold roi') 
        

    % Define frames
    frames = frames;   
    x      = nan(length(frames),1);
    y      = nan(length(frames),1);
    y_flip = nan(length(frames),1);
    
    % Initial central coordinates
    cX = roi0.xCntr;
    cY = roi0.yCntr;
    
    % Area bounds
    minArea = iC.area*0.8;
    maxArea = iC.area.*1.2;
    
elseif strcmp(method,'thresh trans advanced') || ...
       strcmp(method,'thresh trans advanced with mask')
    
    % Define frames
    frames = frames;   
    x      = nan(length(frames),1);
    y      = nan(length(frames),1);
    y_flip = nan(length(frames),1);
    
    % Initial central coordinates
    cX = iC.x;
    cY = iC.y;
    
end



%% Tracking object

% Set up pool for parallel processing
% parpool(8)

nFrames = 1;

% Loop thru frames
parfor i = 1:length(frames)
    
    % Current frame
    cFrame = frames(i);
    
    % Current image
    im = getFrame(vid_path,v,cFrame,imInvert,'gray',imMean,iC.r);
    
    % Eliminate outside of mask
    if ~isempty(xMask)
        
        % Make binary mask of tank region
        bwTank = roipoly(im,xMask,yMask);
        
        % Apply mask
        im(~bwTank) = 255;
        
       % clear bwTank
    end
    
    % Threshold method: find centroid coordinates
    if strcmp(method,'threshold translation')
        
        % Find blob at cX,cY
        [props,bwOut] = findBlobs(im,iC.tVal,'area single',minArea,maxArea,iC.area);

         % Store results
         x(i,1) = props.Centroid(1);
         y(i,1) = props.Centroid(2);
         y_flip(i,1) = size(im,1)-props.Centroid(2);

    end
    
    % Update status
    disp(['trackCenter (' method ') : done ' num2str(i) ' of ' num2str(length(frames))])   

        
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
end

% Store values
Centroid.x        = x;
Centroid.y        = y;
Centroid.y_flip   = y_flip;
Centroid.frames   = frames;


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


