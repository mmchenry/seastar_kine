function convertVideos(action,varargin)
% Converts videos into mp4 format. 
% convertVideo('single',vidFile) performs conversion on single video
% convertVideo('single',vidFile,outPath) performs conversion on single
%             video, places result in outPath
% convertVideo('batch',rootPath) replictes a directory hierarchy,
%     starting at the inPath and converts all MOV files it can find.
% convertVideo('compare',vidFile) compares MOV file to its mp4 counterpart
%
% Only works on mac and linux and requires installation of ffmpeg


%% Parameters

% Default mode
if nargin<1
    action = 'single';
end

% Check action
if ~strcmp(action,'single') && ~strcmp(action,'batch') && ~strcmp(action,'compare') && ...
   ~strcmp(action,'single crop')
    error(['Do not recognize ' action]);
end

% File extension for video format to be converted
fileExt = 'MOV';
    

%% Single mode


if strcmp(action,'single') || ~strcmp(action,'single crop')
    
% Prompt for file path, if not given 
if nargin<2
   [inPath,inFile] =  uigetfile(['*.' fileExt],'Select a video file to be converted');
   
   if ~ischar(inPath) && ~ischar(inFile)
       return
   end
   
   % Define vidoe file
   vidFile = [inPath filesep inFile];
   
   % Name of output file
   outFile = [vidFile(1:end-4) '.mp4'];
   
   clear inPath inFile
   
% Check input, if given
else
    % Use inputs
    vidFile = varargin{1};
    
    % Check for presence
    if ~isfile(vidFile)
        error('Requested file is not present');
    end
    
    % Check file extension
    if ~strcmp(vidFile(end-2:end),fileExt)
        error(['vidFile is not of ' fileExt ' type']);
    end
    
    if length(varargin)>1

        % Destination folder
        outPath = varargin{2};
        
        % Extract file name
        vidFile2 = vidFile(find(vidFile==filesep,1,'last'):end-4);
        
        % Name of output file
        outFile = [outPath filesep vidFile2 '.mp4'];
    else
        % Name of output file
        outFile = [vidFile(1:end-4) '.mp4'];
    end   
end

% Constrcut string with ffmpeg command
% -an - removes audio channel
% -crf determines quality (0 is no compression, 23 is the default. The range is exponential, 
%             so increasing the CRF value +6 results in roughly half the bitrate / file size)

exStr = ['ffmpeg -i ' vidFile ' -an -crf 13 -vf format=gray ' outFile ';'];

%exStr = ['ffmpeg -i ' vid_file ' -an -crf 10 ' vid_file_out];    
  
% Execute ffmpeg at the command line
eval(['! ' exStr])


% Comapre output
% convertVideos('compare',vidFile)


end


%% Comparison of the quality of two video frames
% Useful for evaluating the effects of compression

if strcmp(action,'compare')

% Prompt for file path, if not given 
if nargin<2
   [inFile,inPath] =  uigetfile(['*.' fileExt],'Select MOV file to be compared');
   
   if ~ischar(inPath) && ~ischar(inFile)
       return
   end
   
   % Define vidoe file
   vidFile = [inPath filesep inFile];
   
   clear inPath inFile
   
% Check input, if given
else
    % Use inputs
    vidFile = varargin{1};
    
    % Check for presence
    if ~isfile(vidFile)
        error('Requested file is not present');
    end
    
    % Check file extension
    if ~strcmp(vidFile(end-2:end),fileExt)
        error(['vidFile is not of ' fileExt ' type']);
    end
end    

% Name of other video file
vidFile2 = [vidFile(1:end-4) '.mp4'];

if ~isfile(vidFile2)
    error(['File not found: ' vidFile2]);
end

% Video objects
v1 = VideoReader(vidFile);
v2 = VideoReader(vidFile2);
    
% Read first frames
im1 = rgb2gray(readFrame(v1));
im2 = rgb2gray(readFrame(v2));

imDiff = imsubtract(im1,im2);

figure

warning off
subplot(3,2,1)
imshow(im1)
title('MOV')
subplot(3,2,2)
imhist(im1)

subplot(3,2,3)
imshow(im2)
title('mp4')
subplot(3,2,4)
imhist(im2)

subplot(3,2,5)
imshow(imadjust(imDiff))
title('Difference (contrast enhanced)')
subplot(3,2,6)
imhist(imDiff)
warning on

end


%% Batch mode

if strcmp(action,'batch')

% Prompt for file path, if not given 
if nargin<2
    
   rootPath =  uigetdir(['*.' fileExt],'Select root directory');
   
   if ~ischar(rootPath) 
       return
   end
   
else    
     % Use inputs
    rootPath = varargin{1};

end

% Parent Path, containing both roots
% parentPath = rootPath(1:find(rootPath==filesep,1,'last')-1);

% Set current dirs
cPath  = rootPath;

% Get listing of all folders in path
dirs = dirList(rootPath);

% Replicate directory structure
rootPath_out = makeDirs(dirs,rootPath);

% Get listing of all MOV files
movFiles = movList(dirs,fileExt);

% Loop thru MOV files
for i = 1:length(movFiles)
    
    % Current MOV file
    cMOV = movFiles(i).fullName;
        
    % Current destination
    cDest = [rootPath_out filesep movFiles(i).relPath];
    
    % Current mp4 file
    cMp4 = [cDest filesep movFiles(i).name '.mp4'];
    
    % Run code, only if the file doesn't already exist
    if ~isfile(cMp4)
        
        disp(['Batch mode ' num2str(i) ' of ' num2str(length(movFiles)) ...
          ': running ffmpeg on ' cMOV])
        
        % Run current m-file in single mode
        convertVideos('single',cMOV,cDest); 
    end
end

end


function mFiles = movList(dirs,fileExt)
% Make listing of all mov files in path listing

mFiles = [];
iMov = 1;

for i = 1:length(dirs)
   
    % Get list of MOVs
    a = dir([dirs(i).path filesep '*.' fileExt]);
    
    % Loop thru MOV files, log ones that don't start with '.'
    for j = 1:length(a)  
        if ~strcmp(a(j).name(1),'.')
            
            mFiles(iMov).path      = a(j).folder;
            mFiles(iMov).relPath   = dirs(i).relPath;
            mFiles(iMov).name      = a(j).name(1:end-4);
            mFiles(iMov).ext       = fileExt;
            mFiles(iMov).fullName  = [a(j).folder filesep a(j).name];
            
            iMov = iMov + 1;
        end
    end
    
end



function rootPath_out = makeDirs(dirs,rootPath)
% Replicate directory structure

% Parent Path, containing both roots
parentPath = rootPath(1:find(rootPath==filesep,1,'last')-1);

% Root path stem dir
dirName = rootPath(find(rootPath==filesep,1,'last')+1:end);

% Output dir path
rootPath_out = [parentPath filesep dirName '_mp4'];

% Make output dir, if not there
if ~isfolder(rootPath_out)
    mkdir(rootPath_out);
end

% If there are children . . .
if length(dirs)>1
    
    % Loop thru children
    for i = 2:length(dirs)
        
        % Current directory
        cDir = [rootPath_out filesep dirs(i).relPath];
        
        % Make output dir, if not there
        if ~isfolder(cDir)
            mkdir(cDir);
        end
    end
end


function dirList = dirList(rootPath)
% Returns full listing of directories within the rootPath

% Parent path
parentPath = rootPath(1:find(rootPath==filesep,1,'last')-1);

% Start dList
dirList(1).path    = rootPath;
dirList(1).relPath = '';
dirList(1).parent  = parentPath;
dirList(1).haveChildren = nan;

% Index for logging folders
iLog = 2;

% Index for current directory
iCurr = 1;

% Index for whether to keep running while loop below
runCode = 1;

% Catalog all files to be converted in the hierarchy
while runCode
   
    % If unknown whether there are child . . .
    if isnan(dirList(iCurr).haveChildren)
        
        % Directory listing for current
        a = dir(dirList(iCurr).path);
        
        % Initialize index for child
        iCh = 1;
        
        % Step thru each file and folder
        for i = 1:length(a)
            
            % If folder . . .
            if a(i).isdir && ~strcmp(a(i).name(1),'.')
                
                tmp = [a(i).folder filesep a(i).name];
                
                % Log directory and file name
                chList(iCh).path      = [a(i).folder filesep a(i).name];
                chList(iCh).relPath   = tmp(length(rootPath)+2:end);
                chList(iCh).name      = a(i).name;
                
                % Advance index
                iCh = iCh + 1;
                
                clear tmp
            end
        end
        
        % Log, if no children
        if iCh==1
            
            dirList(iCurr).haveChildren = 0;
            
        % Otherwise . . .
        else
            
            dirList(iCurr).haveChildren = 1;
            
            % Loop thru children
            for i = 1:length(chList)
                
                % Declare the child for current path
                dirList(iCurr).child(i).name   = chList(i).name;
                dirList(iCurr).child(i).logged = 0;
                
                % Add child to dirList
                dirList(iLog).path = chList(i).path;
                dirList(iLog).relPath = chList(i).relPath;
                dirList(iLog).haveChildren = nan;
                
                % Advance index
                iLog = iLog + 1;
            end
        end
        
        clear chList
    end

    % Loop thru all directories collected so far
    for i = 1:length(dirList)
        
        % Stop at first instance of no children listed
        if isnan(dirList(i).haveChildren)
            
            % Set index for current path
            iCurr = i;
            
            
            break
            
        % Stop while-loop, if looked at all directores and all children
        % have been investigated
        elseif i == length(dirList) 
            
            runCode = 0;
        end
    end
    
    % If parent path not yet defined
    if ~isfield(dirList(iCurr),'parent') || isempty(dirList(iCurr).parent)
        
        dirList(iCurr).parentPath = dirList(iCurr).path(...
                             1:find(dirList(iCurr).path==filesep,1,'last')-1);
    end
end





