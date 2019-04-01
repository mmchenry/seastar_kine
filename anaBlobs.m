function varargout = anaBlobs(vid_path,v,opType,varargin)
% Analyze blobs that may be extracted from a video
% B = anaBlobs(vid_path,v,opType)
%  B - structure of blob data
%  vid_path - path to video file or image sequence
%  v - structure of info about video (generated by defineVidObject)
%  opType - designates the type of operation ('G&L props')
%
% B = anaBlobs(vid_path,v,'G&L props',S,blobParam,imInvert,imRoiMean)
% Finds the global and local parameters for blobs
%   im_seq - 3D sequence if images in the roi
%   fr_num - sequence of frame numbers to be included in mean image
%   Centroid - struction denoting the center of the roi (generated by tracker)
%   Rotation - rotation structure generated by tracker
%   blobParam - Structure of blob parameters
%   imInvert - logical that indicates whether to invert images
%
% % B = anaBlobs(vid_path,v,'filter motion',B,winLen,Centroid,...
%       Rotation,blobParam,tVal)
%   tVal - threshold value for blobs
%
% Developed by McHenryLab at UC Irvine


%% Parse inputs

if strcmp(opType,'G&L props')

    % Extract inputs
    %im_seq    = varargin{1};
    Body      = varargin{1};
    %Rotation  = varargin{4};
    blobParam = varargin{2};
    imInvert  = varargin{3};
    dSample   = varargin{4};
    roiM       = varargin{5};
    visSteps  = varargin{6};
 
elseif strcmp(opType,'filter motion')
    
    % Extract inputs
    B         = varargin{1};
    winLen    = varargin{2};
    Body        = varargin{3};
    blobParam = varargin{4};
%     tVal      = varargin{5};
    visSteps  = varargin{5};
    imInvert  = varargin{6};
    imStack   = varargin{7};
    
    dSample = 0;
    
    % Check frame numbers
    for i = 1:length(B)
       B_fr(i,1) = B(i).fr_num;      
    end
    if sum(B_fr'-Body.frames)~=0
        error('Frame numbers do not match between B and Body');
    end
    
    clear B_fr
    
else
    error('opType not recognized');
end

S        = Body.Rotation;
frames   = Body.frames;

if visSteps
    f = figure;
end


%% Loop thru frames ('G&L props')

if strcmp(opType,'G&L props')
    
    %B.frames        = frames;
    %B.blobParam     = blobParam;
    
    
    % Loop thru frames
    for i = 1:length(frames)
        
        % Current frame
        cFrame = frames(i);
        
        % Index for current frame in the data
        %iFrame = find(frames==cFrame);
        
        % Get current mean image
        for j = 1:length(roiM)
            if (cFrame >= roiM(j).frStart) && (cFrame <= roiM(j).frEnd)
                imRoiMean = roiM(j).im;
                imRoiStd  = roiM(j).imStd;
                break
            end
        end
        
        % Current whole frame
        im = getFrame(vid_path,v,cFrame,imInvert,'gray');

        
%         im1 = double(im_roi);
%         imMean = double(imRoiMean);
%         imStd = double(imRoiStd);
%         
%         im2 = im_roi<(imRoiMean+imRoiStd);
        
%         if imInvert
%             im = imcomplement(im);
%         end
        
        % Roi image, mean image subtracted
        [im_roi,bw_mask,bw_roi_mask] = giveROI('stabilized',im,...
            S.roi(i),dSample,S.tform(i),255);
 
        % Subtract mean image, adjust contrast
        im_roi2 = imadjust(imsubtract(imRoiMean,im_roi));
        
        % Find threshold value
        blobParam.tVal = graythresh(im_roi2);
        
        % Find blobs in roi
        [props,bw_roi,areas,xB,yB] = findBlobs(im_roi2,blobParam.tVal,...
            'area and circ',blobParam.areaMin,blobParam.areaMax,blobParam.AR_max);
        
        % Get roi data
        %[bw_mask,im_roi,roi_rect,bw_roi_mask] = giveROI('circular',im,x,y,r,theta,0);
        
        % Blobs in the G FOR
        bw_blobs_G = transCoord2d('bw L2G',S.tform(i),bw_roi,bw_mask,bw_roi_mask);
        
        % Survey blobs
        propsG = regionprops(bw_blobs_G,'Centroid','Area',...
            'MajorAxisLength','MinorAxisLength',...
            'PixelIdxList','PixelList');
        
        % Store blob data
%         B.bwG(:,:,i) = bw_blobs_G;
%         B.bwL(:,:,i) = bw_roi;
        B(i).fr_num = cFrame;
        B(i).propsG = propsG;
        B(i).propsL = props;
        
        if visSteps
            aLevel = 0.5;
            
            % Current whole frame
            imI = getFrame(vid_path,v,cFrame,~imInvert,'gray');

            subplot(2,1,1)
            h = imshow(imI,'InitialMag','fit');
            hold on
            
            % Start with blank
            currIm  = logical(zeros(size(im)));
        
            % Loop thru blobs
            for k = 1:length(propsG)
                % Score pixels with blobs
                currIm(propsG(k).PixelIdxList) = 1;
            end
            
            % Make a truecolor all-green image, make non-blobs invisible
            green = cat(3, zeros(size(im)), ones(size(im)), zeros(size(im)));
            h = imshow(green,'InitialMag','fit');
            %set(h, 'AlphaData', bw_blobs_G.*aLevel)
            set(h, 'AlphaData', currIm.*aLevel)
            title(['Frame ' num2str(cFrame)]);
            
            subplot(2,1,2)            
            h = imshow(imcomplement(im_roi),'InitialMag','fit');
            hold on
            
             % Make a truecolor all-green image, make non-blobs invisible
            green = cat(3, zeros(size(im_roi)), ones(size(im_roi)), zeros(size(im_roi)));
            h = imshow(green,'InitialMag','fit');
            set(h, 'AlphaData', bw_roi.*aLevel)

            pause(0.001)
            
        else
            disp(['anaBlobs (' opType ') : '  num2str(i) ' of ' num2str(length(frames))]);
        end
        
        clear props propsG bw_mask im_roi roi_rect bw_roi_mask bw_roi
        clear areas xB yB imRoiMean im_roi2 imRoiStd
    end   
end


%% Loop thru frames ('filter motion')

if strcmp(opType,'filter motion')
        
    Bin = B;
    clear B
    
    %tVal = blobParam.tVal;

    % Half interval to survey for analysis
    halfIntvl = floor(winLen/2);
    
    % Fill B with placeholders
    for i = 1:length(Bin)
        B(i).fr_num = Bin(i).fr_num;
        %B(i).frIdx  = Bin(i).frIdx;
        B(i).propsG = nan;
        B(i).propsL = nan;
    end
    
    % Frames to analyze
    anaFrames = frames((halfIntvl+1):(length(Bin)-halfIntvl-1));
        
%     % Produce image stack
%     imB = motionImage(vid_path,v,'mask static',frames,Bin,frames,imInvert);
    
    % Loop thru frames to analyze
    for i = 1:length(anaFrames)
        
        % Current frame
        cFrame = anaFrames(i);
        
        % Index for current frame in the data
        %iFrame = Bin(i).frIdx;
        iFrame = find(frames==cFrame,1,'first');

        % Window of frames to analyze
        startFrame    = max([1 cFrame-halfIntvl]);
        endFrame      = min([max(frames) cFrame+halfIntvl]);       
        winFrames     = startFrame:endFrame;
        
%         % Produce image that highlights static elements
%         imB = motionImage(vid_path,v,'mask static',winFrames,Bin,...
%                           frames,imInvert);
        
        % extract current frames 
        bwStack = imStack(:,:,winFrames);

        % Get average image
        imAvg = uint8(sum(double(bwStack),3)./length(winFrames));
        
        % Boost contrast
        %imB = imcomplement(imadjust(imAvg));
        imB = imcomplement((imAvg));

        % Get threshol value for feet
        tVal = 2*graythresh(imB);
        
        % Find blobs in image
        [propsG,bw_im] = findBlobs(imB,tVal,...
                         'area and circ',blobParam.areaMin,...
                         blobParam.areaMax,blobParam.AR_max);
        
        % Get roi data
        %[bw_mask,im_roi,roi_rect,bw_roi_mask] = giveROI('circular',imB,x,y,r,theta,0);
        im_roi = giveROI('stabilized',imB,S.roi(i),dSample, ...
                S.tform(i),0);
         
        % Find blobs in roi
        [propsL,bw_roi] = findBlobs(im_roi,tVal,...
                         'area and circ',blobParam.areaMin,...
                         blobParam.areaMax,blobParam.AR_max);
        
        % Store blob data
        B(iFrame).fr_num = cFrame;
        B(iFrame).frIdx  = iFrame;
        B(iFrame).propsG = propsG;
        B(iFrame).propsL = propsL;

        if visSteps
            
            % Alpha transparency
            aLevel = 0.5;
            
            subplot(2,1,1)
            % Current whole frame
            im = getFrame(vid_path,v,cFrame,imInvert,'gray');
            
            imshow(im,'InitialMag','fit');
            hold on
            
             % Make a truecolor all-green image, make non-blobs invisible
            green = cat(3, zeros(size(imB)), ones(size(imB)), ...
                           zeros(size(imB)));
            h = imshow(green,'InitialMag','fit');
            
            set(h, 'AlphaData', bw_im.*aLevel)
            title(['Frame ' num2str(cFrame)]); 
            for j = 1:length(B(iFrame).propsG)
               scatter(B(iFrame).propsG(j).Centroid(1),...
                       B(iFrame).propsG(j).Centroid(2),'SizeData',200,...
                       'MarkerEdgeColor','b')
            end
            hold off
            
            subplot(2,1,2)
            imshow(giveROI('stabilized',im,S.roi(i),dSample,S.tform(i)),...
                   'InitialMag','fit');
            hold on
             % Make a truecolor all-green image, make non-blobs invisible
            green = cat(3, zeros(size(bw_roi)), ones(size(bw_roi)), ...
                           zeros(size(bw_roi)));
            h = imshow(green,'InitialMag','fit');
            set(h, 'AlphaData', bw_roi.*aLevel)
            
            for j = 1:length(B(iFrame).propsL)
               scatter(B(iFrame).propsL(j).Centroid(1),...
                       B(iFrame).propsL(j).Centroid(2),'SizeData',400,...
                       'MarkerEdgeColor','b')
            end
            hold off
            title(['Frame ' num2str(cFrame)]); 
            
        end
        
        
        % Status report
%         disp(' ')
        disp(['anaBlobs (' opType ') : done ' num2str(i) ' of ' ...
              num2str((length(Bin)-halfIntvl))])
%         disp(' ');
        
        clear im_roi im bw_roi bw_im propsG propsL
    end 
        
end

%% Outputs

varargout{1} = B;



