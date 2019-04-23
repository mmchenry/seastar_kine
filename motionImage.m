function varargout = motionImage(vid_path,v,imType,varargin)
% Creates a single image from a video sequence
% im = motionImage(vid_path,v,imType,imProcessing)
%  im - output image
%  vid_path - path to video file or image sequence
%  v - structure of info about video (generated by defineVidObject)
%  imType - designates the type of image ('mean','mean and roi')
%  imProcess - manipulation of video images ('none','enhance contrast')
%
%  im = motionImage(vid_path,v,'mean',imProcess,imInvert,fr_num)
%  Returns mean image of provided frame numbers.
%  fr_num - column vector of frame numbers to be included in mean image,
%  or single integer for the max number of frames at an equal interval
%
%  im = motionImage(vid_path,v,'mean',imProcess,imInvert,fr_num)
%  Returns mean image 90% brightest intensity over recording. Works well
%  for dark object moving thru light background.
%  fr_num - column vector of frame numbers to be included in mean image,
%  or single integer for the max number of frames at an equal interval
%
%  im = motionImage(vid_path,v,'mean roi',imProcess,fr_num,S,dSample)
%  Returns mean image for roi at even intervals or frnums. 
%   S - structure for roi coorindate trasnformation (generated by
%       defineSystem2d('roi'))
%   dSample - logical that indicates whether to downsample the roi image
%
% im = motionImage(vid_path,v,'bw static',fr_num,B)
% Returns an image of the static elements in a binary image sequence
%   B - structure of blob data

% Developed by McHenryLab at UC Irvine


%% Process inputs

% imType
if strcmp(imType,'mean') || ...
   strcmp(imType,'mean bright') || ...
   strcmp(imType,'mean color')
    
    % Image processing
    imProcess = varargin{1};
    
    % Image processing
    imInvert = varargin{2};
    
    % Extract frame numbers
    if nargin>5
        if length(varargin{3})>1
            fr_num = varargin{3};
        else
            % Max number of frames to analyze
            maxFrames = varargin{3};
            fr_num = [];
        end
    else
        maxFrames = 100;
        fr_num = [];
    end
    

elseif strcmp(imType,'mean roi') 
    
    % Image processing
    imProcess = varargin{1};
    imInvert  = varargin{2};
    Body      = varargin{3};
    dSample   = varargin{4};
    
    % Extract frame numbers to analyze
    if length(varargin)>4
        
        fr_num = varargin{5};
   else
        maxFrames = 100;
        fr_num = [];
    end
    
    S      = Body.Rotation;
    frames = Body.frames;
    
    clear Body
   
elseif strcmp(imType,'im static') 
    
    fr_num   = varargin{1};
    S        = varargin{2};
    imInvert = varargin{3};
    
    imProcess = 'none';
    
elseif strcmp(imType,'bw static') 
    
    fr_num = varargin{1};
    B      = varargin{2};
    frames = varargin{3};
    imInvert = varargin{4};
    imProcess = [];
    %imInvert = 0;
    
    if length(fr_num)<2
        error('fr_num must be a vector of at least 2 numbers')
    end
    
elseif strcmp(imType,'mask static') 
    
    Body      = varargin{1};
    B         = varargin{2};
    imInvert  = varargin{3};
    iC        = varargin{4};
    imProcess = [];
    %imInvert = 0;
    
    % Extract data from Body
    fr_num  = Body.frames;
    frames  = Body.frames;
    
    % Get center points
    xCntr = Body.xCntr;
    yCntr = Body.yCntr;
    
    % Threshold for body
    tVal = iC.tVal;
    
    

else
    error('Do not recognize imType');
end

% imProcessing
if strcmp(imProcess,'none') || isempty(imProcess)
    %Do nothing
    
elseif strcmp(imProcess,'enhance contrast')
    %Do nothing
else
    error('Do not recognize imProcessing');
end


%% Preliminaries

% % Check for image sequence
% if ~isdir(vid_path)
%     error('This function requires that the video is saved as a series of images')
% end
   
% If frame numbers not provided . . .
if isempty(fr_num)
    
    % Number of frames
    numFrames = range(frames);
    
    % Define list of frame numbers, depending on max number of frames
    % requested
    if numFrames > maxFrames
        dframe    = floor(numFrames/maxFrames);
        fr_num    = frames(1):dframe:frames(end);
        
        clear dframe frame1
    else
        %fr_num = 1:numFrames;
        fr_num = v.UserData.FirstFrame:v.UserData.LastFrame;
    end
end
    
% if strcmp(imType,'bw static')
%     % Create sum image based on first frame
%     imCurr = getFrame(vid_path,v,fr_num(1),imInvert,'gray');    
% else     
%     % Create sum image based on first frame
%     imCurr = getFrame(vid_path,v,fr_num(1),imInvert,'gray');  
% end
% 
% clear imCurr 


%% Step thru frames

% Loop through frames
for i = 1:length(fr_num)

    % Current frame number
    cFrame = fr_num(i);
    
    % Index for current frame in the data
    iFrame = find(cFrame==frames,1','first');

    % Load current image
    if strcmp(imType,'bw static') || strcmp(imType,'mask static')
        % Create sum image based on first frame
        imCurr = getFrame(vid_path,v,cFrame,imInvert,'gray');    
            
    elseif strcmp(imType,'im static')
        % Create sum image based on first frame
        imCurr = getFrame(vid_path,v,cFrame,imInvert,'gray');    
        
    elseif strcmp(imType,'mean') || strcmp(imType,'mean bright')
           
        % Get current frame
        imCurr       = getFrame(vid_path,v,cFrame,imInvert,'gray');   
        
    elseif strcmp(imType,'mean color')
        
        % Get current frame
        imCurr       = getFrame(vid_path,v,cFrame,imInvert,'rgb');  
    
    elseif strcmp(imType,'mean roi') 
        % Get current frame
        im       = getFrame(vid_path,v,cFrame,imInvert,'gray');     
        
        % Get current roi
        imCurr = giveROI('stabilized',im,S.roi(iFrame),dSample,S.tform(iFrame));
        
    end
    
    % Convert to grayscale
    %imCurr = rgb2gray(imCurr);
    
    if strcmp(imType,'bw static') || strcmp(imType,'mask static') 
        % Start with blank
        currIm  = logical(zeros(size(imCurr)));

    end

    % Enhance contrast, if requested
    if strcmp(imProcess,'enhance contrast') 
        imCurr = imadjust(imCurr);
    end

    % If mean image(s) . . .
    if strcmp(imType,'mean') || ...
       strcmp(imType,'mean roi') || ...
       strcmp(imType,'mean bright')

        % Add current to total
        %imSum  = imSum + double(imCurr);    
        imStack(:,:,i) = double(imCurr);   
        
        if 0
            imshow(imCurr)
            title(['Frame ' num2str(cFrame)])
            pause(1)
        end
           
        
    elseif strcmp(imType,'mean color')
        
        imStack.R(:,:,i) = imCurr(:,:,1);
        imStack.G(:,:,i) = imCurr(:,:,2);
        imStack.B(:,:,i) = imCurr(:,:,3);       
        
    elseif strcmp(imType,'bw static') 
        
        % Loop thru blobs
        for k = 1:length(B(i).propsG)
            
            % Score pixels with blobs
            currIm(B(i).propsG(k).PixelIdxList) = 1;            
        end
       
        % Store resulting image
        bwStack(:,:,i) = currIm;
              
    elseif strcmp(imType,'mask static') 
        
        % Get binary of body
        bwMask = ~imbinarize(imCurr,tVal);
        
        se = strel('disk',12);
        bwMask = imdilate(bwMask,se);
        bwMask = imerode(bwMask,se);
        
        % Fill holes
        bwMask = imfill(bwMask,'holes');
        
        % Select body 
        bwMask = bwselect(bwMask,xCntr(iFrame),yCntr(iFrame));
        
        % Loop thru blobs
        for k = 1:length(B(iFrame).propsG)        
            % Score pixels with blobs
            currIm(B(iFrame).propsG(k).PixelIdxList) = 1;            
        end
        
        % White out all non-foot pixels
        im_tmp = imCurr;
        im_tmp(~currIm) = 255;
        
        % White out non-body pixels
        im_tmp(~bwMask) = 255;

        % Start with blank
%         tmpIm  = logical(zeros(size(imCurr)));
        
        
        % Display current frame
        if 0
            aLevel = 0.2;
            warning off
            h = imshowpair(currIm,imCurr);
            warning on
            
%             h = imshow(imCurr,'InitialMag','fit');
%             hold on
%             green = cat(3, zeros(size(imCurr)), ones(size(imCurr)), ...
%                         zeros(size(imCurr)));
%             h = imshow(green,'InitialMag','fit');
%             set(h, 'AlphaData',im_tmp.*aLevel)
            hold off
            %            title(['Frame ' num2str(cFrame)]);
            
        end
        
        % Store resulting image
        bwStack(:,:,i) = im_tmp;
        
    elseif strcmp(imType,'im static')  
        
        % Store resulting image
        bwStack(:,:,i) = imCurr;
        
    end
    
    % Clear fopr next loop
    clear imCurr currIm
    
    % Update status
    disp(['      motionImage (' imType ') : ' num2str(i) ' of ' num2str(length(fr_num))])
    
end

% Calculate mean image(s)
if strcmp(imType,'mean') || ...
   strcmp(imType,'mean roi')
    
    imOut  = uint8(round(mean(imStack,3)));
    stdOut = uint8(round(std(imStack,1,3)));
    
elseif strcmp(imType,'mean color')
    
    imOut(:,:,1) = uint8(round(mean(imStack.R,3)));
    imOut(:,:,2) = uint8(round(mean(imStack.G,3)));
    imOut(:,:,3) = uint8(round(mean(imStack.B,3)));
    
elseif strcmp(imType,'mean bright')
    
    imOut = uint8(round(quantile(imStack,0.9,3)));

  
elseif strcmp(imType,'im static')  
    
    error('Need to finish this line of code');
    
% Create image from motion score    
elseif strcmp(imType,'bw static') 
    
    % image of the pixels that are most static
%     imOut = imcomplement(uint8(sum(double(bwStack),3)./length(fr_num).*255));
imOut = bwStack;
    
% Create image from motion score    
elseif strcmp(imType,'mask static') 
    
    imOut = bwStack;
    
    % Get average image
    %imAvg = uint8(sum(double(bwStack),3)./length(fr_num));
    
    % Boost contrast
    %imOut = imcomplement(imadjust(imAvg));    
    
end


%% Output results

varargout{1} = imOut;

if nargout>1
    varargout{2} = stdOut;
end



% function imStable = roiImage(im,Centroid,Rotation,iNum,dSample)
% % Returns stabilized roi image
% 
% % Extract current parameters
% x = Centroid.x_pix(iNum);
% y = Centroid.y_pix(iNum);
% r = Centroid.r_pix;
% theta = Centroid.theta;
% tform = Rotation(iNum).tform_roi;
% 
% % Extract roi
% [bw_mask,im_roi,roi_rect,bw_roi,imStable] = giveROI('circular',...
%     im,x,y,r,theta,dSample,tform);
% 


