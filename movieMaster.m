function movieMaster(dataPath,vidPath,action)


%% Parameters

if strcmp(action,'deep movie')
    
    
else
    error(['Do not recognize ' action])
end

% Invert image
imInvert = 1;


%% Paths

% Root paths
paths = givePaths;

% Paths for current sequence
currDataPath = [paths.data filesep dataPath];
currVidPath  = [paths.vid filesep vidPath];


% Check video path 
if ~isfile(currVidPath)
    error(['Video file does not exist at ' currVidPath]);
end

% Check data path 
if ~isfolder(currDataPath)
    error(['Data folder does not exist at ' currDataPath]);
end

% Load video info (v)
v = defineVidObject(currVidPath);

 % Load body kinematics (Body)
load([currDataPath filesep 'Body.mat'])

S        = Body.Rotation;
% frames   = Body.frames;

aLevel = 0.2;


%% Movie for DeepLabCut


if strcmp(action,'deep movie')
    
% Downsample
dSample = 0;

% Path for mean images
mPath = [currDataPath filesep 'mean_images'];

% Load initial conditions (iC)
load([currDataPath filesep 'Initial conditions'])

% Path to save video
iLast = find(currVidPath==filesep,1,'last');
saveVidPath = [currVidPath(1:iLast) 'roi_' currVidPath((iLast+1):end-4)];

% Set up output video file
vOut = VideoWriter([saveVidPath '.mp4'],'MPEG-4');
vOut.Quality = 50;
open(vOut)

% Blob data files
%blobDir = dir([currDataPath filesep 'blobs' filesep 'blobs_*.mat']);
maskDir  = dir([currDataPath filesep 'mask_static' filesep 'mask_*.mat']);
%ftDir     = dir([currDataPath filesep 'foot_blobs' filesep 'foot_*.mat']);
ftDir     = dir([currDataPath filesep 'blobs' filesep 'blobs_*.mat']);

% Mean images for local FOR
meanDir = dir([currDataPath filesep 'mean_images' filesep 'mean_*.mat']);

% Get frame numbers from ftDir
for i = 1:length(ftDir)
    frames(i,1) = str2num(ftDir(i).name(end-9:end-4));   
end

% Sort frame numbers
[frames,idx] = sort(frames);

% Read frames corresponding to each mean image
for i = 1:length(meanDir)
    load([currDataPath filesep 'mean_images' filesep meanDir(i).name])
    meanDir(i).frames = roiM.frames;
end

disp(' ');disp(['Writing new roi video to ' saveVidPath  ' . . .'])

clear iLast saveVidPath i

% f = figure;


% Loop thru frames
for i = 1:length(frames)
    
    % Current frame
    cFrame = frames(i);
    
    % Load B_ft strcuture
    load([currDataPath filesep 'blobs' filesep  ftDir(idx(i)).name])
    
    roiM = [];
    
    % Find mean image (roiM)
    for j = 1:length(meanDir)
        if cFrame>=min(meanDir(j).frames) && cFrame<=max(meanDir(j).frames)
            load([currDataPath filesep 'mean_images' filesep meanDir(j).name])
        end
    end
    
    % Check that mean image found
    if isempty(roiM)
        error('No mean image for current frame')
    end
        
    % Get data from S structure
    [roi,tform] = returnS(S,i);
        
    % Current whole frame
    im = getFrame(currVidPath,v,cFrame,imInvert,'gray',[],iC.r);

    % Roi image, mean image subtracted
    [im_roi,bw_mask,bw_roi_mask] = giveROI('stabilized',im,...
            roi,dSample,tform);
        
    bw_roi = addmask(im_roi);
         
    % Subtract background
    im_roi = imsubtract(roiM.im,im_roi);
    
     % Remove background stuff
    im_roi(~bw_roi) = 0;
    
    % Check that the sea star is there
    if max(im_roi(:))==0
        error('No sea star is in the video frame');
    end
    
   % Subtract mean image, adjust contrast
   im_roi = imadjust(im_roi,[0 0.6*double(max(im_roi(:)))/255],[0 1]);    
   
  
    % Start with blank
%     currIm  = logical(zeros(size(im_roi)));

%     % Loop thru blobs
%     for k = 1:length(B.propsL)
%         % Score pixels with blobs
%         currIm(B.propsL(k).PixelIdxList) = 1;
%     end
    
%     h = imshow(im_roi,'InitialMag','fit');
%     hold on 
%     
%     % Make a truecolor all-green image, make non-blobs invisible
%     white = cat(3, ones(size(im_roi)), ones(size(im_roi)), ones(size(im_roi)));
%     h = imshow(white,'InitialMag','fit');
%     set(h, 'AlphaData', currIm.*aLevel)

    
    writeVideo(vOut,im_roi);
    
    disp(['    Frame ' num2str(i) ' of ' num2str(length(frames)) ' written'])
end


end






function bw = addmask(im)
% Mask with distance map (trims fins)

xC = size(im,2)/2;
yC = size(im,1)/2;

% Enhance contrast
im1 = imadjust(imcomplement(im));

% Binary image
bw = im2bw(im1,graythresh(im1));

% Fill holes
bw = imfill(bw,'holes');

% Dilate, white out outside
se = strel('disk',3,6);
bw = imerode(bw,se);
im(~bw) = 255;

% Select object in center
bw = bwselect(bw,xC,yC);
% 
% % Distance map
% bwD = bwdist(bw);
% 
% % Max distance
% maxDist = max(bwD(:));
% 
% % Threshold distance
% threshDist = maxDist/3;
% 
% % Refine binary as above-threshold images
% bw = bwD>threshDist;
% 
% % Dilate, white out outside
% se = strel('disk',3,4);
% bw = imdilate(bw,se);
% im(~bw) = 255;






function [roi,tform] = returnS(S,i)

roi    = S.roi(i);
tform  = S.tform(i);
