function im = getFrame(vid_path,v,fr_num,imInvert,clrMode,imMean,padLen)  
% Reads image file from a video. 
%   v - video structure
%   fr_num - frame number desired
%   vid_path - path to video (use path in v structure as default)
%   imInvert - logical that indicates whether to invert a grayscale image
%   clrMode - 'rgb' or 'gray' format for image 
%   imMean - (optional) mean image of video sequence

%
% Developed by McHenryLab at UC Irvine

%% Parse inputs

% Default for image invert
    if nargin < 4
        imInvert = 0;
    end
    
    % Default color mode
    if nargin < 5
        clrMode = 'rgb';
    end
    
    if nargin < 6
        imMean = [];
    end
    
    if nargin < 7
        padLen = 0;
    end

    % Check requested frame number
    if fr_num>v.UserData.LastFrame
        error(['Video sequence does not have a frame ' num2str(fr_num)]);
    end
    
    % Frame index
    iFrame = fr_num-v.UserData.FirstFrame + 1;
    %iFrame = fr_num;
    
    if isfield(v.UserData,'FileInfo')
        % Get filename and extension
        fName = v.UserData.FileInfo(iFrame).name;
        ext   = fName(find(fName=='.')+1:end);
        fName = fName(1:(find(fName=='.')-1));
        
        % Get frame number
        frNum = str2num(fName((find(fName=='_')+1):end));
        
        % Check for match
        if frNum~=fr_num
            error('file numbering does not match the video info');
        end
    end


%% Read frame

% If it is an image sequence . . .
if isfield(v.UserData,'FileInfo')

    % Read image
    im = imread([vid_path filesep v.UserData.FileInfo(iFrame).name]);    
        

% Read frame from a movie file . . .
else
  
%     warning off
%     % Adjust items to new 'v'
%     v = VideoReader(vid_path);
%     warning on
    
    % Set current time to last frame
    if fr_num/v.FrameRate==v.Duration
        v.CurrentTime = (fr_num-1)./v.FrameRate;
        
    % Set current time to other frame
    else
        v.CurrentTime = fr_num./v.FrameRate;
    end
    
     % Read next available frame
    im = readFrame(v);
    
end

%% Modify image

% Convert to grayscale
if strcmp(clrMode,'gray') && size(im,3)==3
    im = rgb2gray(im);
end

% Invert, if requested
if imInvert==1
    im = imcomplement(im);
end

% Subtract mean image, if present
if ~isempty(imMean)
    %im = imsubtract(im,imMean);
    
    % Get the compliment
    imMean  = imcomplement(imMean);
    im      = imcomplement(im);
    
    im = uint8(double(im)-double(imMean));
    im = imcomplement(im);
    
    %im = uint8(imadjust(imcomplement(imsubtract(imMean,im))));
    
end

% Add padding
if padLen~=0
    
    % Increase pad slightly
        padLen = padLen+2;
        
%     if strcmp(clrMode,'color')
        
    if mean(im(:))>(256/2)
        qVal = 0.75;
    else
        qVal = 0.25;
    end

    for i = 1:size(im,3)
        % Color of padding
        padClr = double(round(quantile(...
                        reshape(im(:,:,i),size(im,1)*size(im,2),1),qVal)));
        
        % Add pad to sides, top and bottom
        im2  = [padClr.*ones(size(im,1),padLen) im(:,:,i) ...
                      padClr.*ones(size(im,1),padLen)];
        im2 = [padClr.*ones(padLen,size(im2,2)); im2];
        im2 = [im2; padClr.*ones(padLen,size(im2,2))];
        
        % Store in im3
        im3(:,:,i) = im2;
        
        clear im2
    end
%     else
%         % Color of padding
%         padClr = double(round(quantile(im(:),0.75)));
%         
%         % Increase pad slightly
%         padLen = padLen+2;
%         
%         % Add pad to sides, top and bottom
%         im2 = [padClr.*ones(size(im,1),padLen) im padClr.*ones(size(im,1),padLen)];
%         im2 = [padClr.*ones(padLen,size(im2,2)); im2];
%         im2 = [im2; padClr.*ones(padLen,size(im2,2))];
%     end
    % Transfer image
     im = im3;
end


