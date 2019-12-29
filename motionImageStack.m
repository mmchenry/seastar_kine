function varargout = motionImageStack(vid_path,v,imType,varargin)
% Creates a series of images from a video sequence

% Developed by McHenryLab at UC Irvine


%% Process inputs

% imType
if strcmp(imType,'mask static') 
    
    Body        = varargin{1};
    blobPath    = varargin{2};
    imInvert    = varargin{3};
    iC          = varargin{4};
    motionPath  = varargin{5};
    winLen      = varargin{6};
    
    imProcess = [];
    %imInvert = 0;
    
    % Extract data from Body
%     fr_num  = Body.frames;
    frames  = Body.frames;
    
    % Get center points
    xCntr = Body.xCntr;
    yCntr = Body.yCntr;
    
    % Threshold for body
    tVal = iC.tVal;
    
    if ~isfolder(motionPath)
        mkdir(motionPath)
    end

    % File name
    fName = 'mask_static';
    
    % Half interval to survey for analysis
    halfIntvl = floor(winLen/2);
    
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

% Frame numbers to analyze
fr_num = (frames(1)+halfIntvl):(frames(end)-halfIntvl);


%% Step thru frames


imStack = [];
frCurr =[];

% Loop through frames
for i = 1:length(fr_num)
    
    % Current frame number
    cFrame = fr_num(i);
    
    % Current file name (for reading)
    cNumStr    = ['00000' num2str(cFrame)];
    blobName   = ['blobs_' cNumStr(end-5:end)];
    
    % Current file name (for writing)
    cName   = [fName '_' cNumStr(end-5:end)];
    
    % Index for current frame in the data
    iFrame = find(cFrame==frames,1','first');
    
    % Load B structure
    B = loadB([blobPath filesep blobName]);
    
    % Update stack
    [imStack,frCurr] = addToStack(vid_path,v,frCurr,imStack,cFrame,...
        halfIntvl,B,tVal,xCntr,yCntr,imInvert,iC);
    
    % Get average image
    imBlur.im = uint8(sum(double(imStack),3)./length(frCurr));
    
    % Store frame numbers
    imBlur.frames = frCurr;
    
    % Write to disk
%     saveImageData([motionPath filesep cName], imBlur)
    save([motionPath filesep cName],'-v7.3','imBlur');
    
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
    
    % Update status
    disp(['      motionImageSeries (' imType ') : ' num2str(i) ' of ' num2str(length(fr_num))])   
    
end



function [imStack,frCurr] = addToStack(vid_path,v,frCurr,imStack,cFrame,...
                                 halfIntvl,B,tVal,xCntr,yCntr,imInvert,iC)

% If no data yet . . .                             
if isempty(imStack)
    
    % Interval of frames to include
    startFrame = cFrame - halfIntvl;
    endFrame   = cFrame + halfIntvl;
    winFrames  = startFrame:endFrame;
    
    % Current frames include all winFrames
    frCurr = winFrames;
    
    % Load stack of images
    for i = 1:length(winFrames)
        
        % Create sum image based on first frame
        im = getFrame(vid_path,v,winFrames(i),imInvert,'gray',[],iC.r);
        
        % Modify image, store result
        imStack(:,:,i) = imModify(im,B,tVal,xCntr,yCntr);
    end
    
% If adding data . . .
else
    % New end frame
    endFrame   = cFrame + halfIntvl;
    
    % Create sum image based on first frame
    im = getFrame(vid_path,v,endFrame,imInvert,'gray',[],iC.r);
    
    % Modify image, store result
    imEnd = imModify(im,B,tVal,xCntr,yCntr);
    
    % Trim start frame
    imStack = imStack(:,:,2:end);
    
    % Add new frame to end
    imStack(:,:,end+1) = imEnd;
end

% function saveImageData(sPath, imAvg)
%     % Save data for current average image
%      save(sPath,'-v7.3','imAvg')   

function im = imModify(imCurr,B,tVal,xCntr,yCntr)

% Start with blank
currIm  = logical(zeros(size(imCurr)));

% Get binary of body
bwMask = ~imbinarize(imCurr,tVal);

% Dilation and erosion
se = strel('disk',12);
bwMask = imdilate(bwMask,se);
bwMask = imerode(bwMask,se);

% Fill holes
bwMask = imfill(bwMask,'holes');

% Select body
bwMask = bwselect(bwMask,xCntr,yCntr);

% Loop thru blobs
for k = 1:length(B.propsG)
    % Score pixels with blobs
    currIm(B.propsG(k).PixelIdxList) = 1;
end

% White out all non-foot pixels
im = imCurr;
im(~currIm) = 255;

% White out non-body pixels
im(~bwMask) = 255;


function saveData(sPath,imStack)

 % Write image data
 save(sPath,'-v7.3','imStack')


function  B = loadB(pathB)
% Load B (blob data)
load(pathB,'B'); 


