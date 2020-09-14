function varargout = aniData(vid_path,v,data_path,movie_file,...
                             imInvert,opType,varargin)
% Animates blobs that may be extracted from a video
% M = anaBlobs(vid_path,v,imInvert,display,opType)
%  M - movie object
%  vid_path - path to video file or image sequence
%  v - structure of info about video (generated by defineVidObject)
%  opType - designates the type of operation ('blobs G&L')
%  imInvert - logical that designates whether to invert the image
%
% M = anaBlobs(vid_path,v,'blobs G&L',B,imVis)
% Overlays blob data onto video frames in local and global FOR
%   B - Blob data structure (created by anaBlobs)
%   imVis - logical indicating whether to display frames as they are put
%   together.
%
% M = anaBlobs(vid_path,v,'Centroid & Rotation',S,imVis)
%    S - structure that defines coordinate transformation bwtn roi and
%    frame for whole movie (created by defineSystem2d)

%% Parameters

alphaLevel = 0.5;

if ~isunix
    % Set up output video file
    vOut = VideoWriter([data_path filesep movie_file '.mp4'],'MPEG-4');
    vOut.Quality = 50;
else
    % Set up output video file
%     vOut = VideoWriter([data_path filesep movie_file '.avi'],'Uncompressed AVI');
    vOut = VideoWriter([data_path filesep movie_file '.avi'],'Motion JPEG AVI');
    vOut.Quality = 75;
end

open(vOut)


%% Parse inputs

if strcmp(opType,'blobs G&L')
    
    % Blob data structure
    B = varargin{1};
      
      j = 1;
    for i = 1:length(B)
        if ~isempty(B(i).fr_num)
            frames(j,1) = B(i).fr_num;
            j = j + 1;
        end
    end
    
    if nargin > 5
        imVis = varargin{2};
    else
        imVis = 1;
    end
    
elseif strcmp(opType,'blobs L simple')
   
    S          = varargin{1};
    tVal       = varargin{2};
    areaMin    = varargin{3};
    areaMax    = varargin{4};
    imVis      = varargin{5};
    imRoiMean  = varargin{6};
    dSample    = varargin{7};
    
    frames = S.frames;      
    
elseif strcmp(opType,'Centroid tracking')
       
    %S       = varargin{1};
    Body = varargin{1};
    frames = Body.frames;
    %frames = S.frames;    
    numroipts = length(Body.roi(1).theta);
    
    imVis = varargin{2};
    iC = varargin{3};
    
    fColor = [1 1 1];
    fPos = [1 1 1150 1340];
    
elseif strcmp(opType,'Feet')
    
    Body     = varargin{1}; 
    imVis    = varargin{2};
    B_ft     = varargin{3};
    
    % Start index
    j = 1;
    
    % Loop thru frames
    for i = 1:length(B_ft)
        
       % If there is foot blob data . . .
       if ~isempty(B_ft(i).frIdx)
           
           % Store frame number
           frames(j,1) = B_ft(i).fr_num;
           
           % Index in the Body data
           idx = find(Body.frames==frames(end),1,'first');
           
           % Store form blob data
           propsG{j} = B_ft(i).propsG;
           propsL{j} = B_ft(i).propsL;
           
           % Store from Body
           roi(j,1)     = Body.Rotation.roi(idx);
           xCntr(j,1)   = Body.xCntr(idx);
           yCntr(j,1)   = Body.yCntr(idx);
           tform(j,1)   = Body.Rotation.tform(idx);
           
           % Arm coordinates
           xArm(j,:)   = Body.xArmG(idx,:);
           yArm(j,:)   = Body.yArmG(idx,:);
           
           % Advance index
           j = j + 1;
       end
       
       
    end
      
    numroipts = length(Body.Rotation.roi(1).xPerimL);    

    fColor = 0.2.*[1 1 1];
    
    % Clear 
    clear Body j idx B_ft
    
elseif strcmp(opType,'Local feet')
    
     
    Body         = varargin{1}; 
    visSteps     = varargin{2}; 
    F            = varargin{3};
    iC           = varargin{4};
    tText        = varargin{6};
    
    % Listing of data 
    aData = dir([currDataPath filesep 'foot_blobs' filesep 'foot_blobs*.mat']);
    
    % Loop thru frames to store properties
    for i = 1:length(aData)
           % Store frame number
           frames(i,1) = str2num(aData(i).name(end-9:end-4));
    end
    
    % List frame numbers in order
    frames = sort(frames);
    
    % Loop thru frames
    for i = 1:length(frames)
        
        % Index in the Body data
        idx = find(Body.frames==frames(i),1,'first');
        
        % Find matching coordinates
        if ~isempty(idx)
            % Store from Body
            roi(i,1)     = Body.Rotation.roi(idx);
            xCntr(i,1)   = Body.xCntr(idx);
            yCntr(i,1)   = Body.yCntr(idx);
            tform(i,1)   = Body.Rotation.tform(idx);  
        else
            error('No match in Body data for current foot data');
        end
        
        % Index for foot coordinates
        n = 1;
        
        % Find matches in foot data
        for j = 1:length(F)
            
            % Index of mathcing frames
            idx2 = F(j).frames==frames(i);
            
            % Log foot cooordinates, if there
            if sum(idx2)==1
                xTmp(n,1) = F(j).xL(idx2);
                yTmp(n,1) = F(j).yL(idx2);
                
                n = n + 1;
                
            elseif sum(idx2)>1
                error('More than one matching frame, somehow')
            end
        end
        
        % Store foot results
        if n > 1
            xFt{i} = xTmp;
            yFt{i} = yTmp;
        else
           xFt{i} = nan;
           yFt{i} = nan;
        end
    end
      
    % Number of points in ROI
    numroipts = length(Body.Rotation.roi(1).xPerimL);    
    
    % Clear 
    clear Body j idx B_ft
    
    % Number of rows and columns in each fig
    nRow = 2;
    nCol = 1;
    
    % Number of panels for each frame
    pNumAdvance = 1;
    
    % Figure color
    fColor = 0.2.*[1 1 1];


elseif strcmp(opType,'Individual feet, local')
    
    Body     = varargin{1}; 
    imVis    = varargin{2};
    F        = varargin{3};
    iFrames  = varargin{4};
    iC       = varargin{5};

    % Listing of raw data for feet 
    aData = dir([data_path filesep 'foot_blobs' filesep 'foot_blobs*.mat']);
    
    % Loop thru filenames to get frame numbers
    for i = 1:length(aData)
           % Store frame number
           frames(i,1) = str2num(aData(i).name(end-9:end-4));
    end
    
    % List frame numbers in order
    [frames,iFile] = sort(frames);
    
    % Loop thru frames
    for i = 1:length(frames)
        
        % Number of feet in frame
        n = 0;
        
        % Index for body data
        iMatchBody = Body.frames==frames(i);
        
        % Arm coordinates
%         xArm   = Body.xArmL;
%         yArm   = Body.yArmL;
        
        % Get roi and transformation values
        roi(i,1)     = Body.Rotation.roi(iMatchBody);
        tform(i,1)   = Body.Rotation.tform(iMatchBody); 
        
        % Body center
%         xCntr(i,:)   = Body.xCntr(iMatchBody);
%         yCntr(i,:)   = Body.yCntr(iMatchBody);
        
        % Loop thru feet
        for j = 1:length(F)
            
            % Index fo current frame to F data
            iMatch = F(j).frames==frames(i);
 
            % Check for repeats
            if sum(iMatch)>1
                error('Multiple instances of same frame for a foot');
            end
            
            % If there is a match, store coordinates
            if max(iMatch)
                % Advance index
                n = n + 1;
                
                % Store coordinates and colors
                x{i}(n,1)      = F(j).xL(iMatch);
                y{i}(n,1)      = F(j).yL(iMatch);
                clr{i}(n,:)    = F(j).clr(1,:);
            end  
        end
    end
    
    % Figure color and positon
    fColor = 0.2.*[1 1 1];
    fPos   = [ 1 530  2002  995];
    
    % Clear 
    clear Body j idx B_ft
    
  
elseif strcmp(opType,'Global feet')
    
    Body     = varargin{1}; 
    imVis    = varargin{2};
    B_ft     = varargin{3};
    
    % Start index
    j = 1;
    
    % Loop thru frames
    for i = 1:length(B_ft)
        
       % If there is foot blob data . . .
       if ~isempty(B_ft(i).frIdx)
           
           % Store frame number
           frames(j,1) = B_ft(i).fr_num;
           
           % Index in the Body data
           idx = find(Body.frames==frames(end),1,'first');
           
           % Store form blob data
           propsG{j} = B_ft(i).propsG;
           propsL{j} = B_ft(i).propsL;
           
           % Store from Body
           roi(j,1)     = Body.Rotation.roi(idx);
           xCntr(j,1)   = Body.xCntr(idx);
           yCntr(j,1)   = Body.yCntr(idx);
           tform(j,1)   = Body.Rotation.tform(idx);
           
           % Arm coordinates
           xArm(j,:)   = Body.xArmG(idx,:);
           yArm(j,:)   = Body.yArmG(idx,:);
           
           % Advance index
           j = j + 1;
       end
       
       
    end
      
    numroipts = length(Body.Rotation.roi(1).xPerimL);    

    % Figure color and positon
    fColor = 0.2.*[1 1 1];
    fPos = [ 1 530  2002  995];
    
    % Clear 
    clear Body j idx B_ft F
    
    
elseif strcmp(opType,'Individual feet, pretty')
    
    Body     = varargin{1}; 
    imVis    = varargin{2};
    F        = varargin{3};
    iFrames  = varargin{4};
    iC       = varargin{5};

    % Frame numebrs to analyze
    frames   = Body.frames(iFrames)';
    
    for i = 1:length(frames)
        
        % Number of feet in frame
        n = 0;
        
        % Index for body data
        iMatchBody = Body.frames==frames(i);
        
        % Arm coordinates
        xArm(i,:)   = Body.xArmG(iMatchBody,:);
        yArm(i,:)   = Body.yArmG(iMatchBody,:);
        
        % Body center
        xCntr(i,:)   = Body.xCntr(iMatchBody);
        yCntr(i,:)   = Body.yCntr(iMatchBody);
        
        % Loop thru feet
        for j = 1:length(F)
            
            % Index of current frame to F data
            iMatch = F(j).frames==frames(i);
 
            % Check for repeats
            if sum(iMatch)>1
                error('Multiple instances of same frame for a foot');
            end
            
            % If there is a match, store coordinates
            if max(iMatch)
                % Advance index
                n = n + 1;
                
                % Store coordinates and colors
                x{i}(n,1)      = F(j).xG(iMatch);
                y{i}(n,1)      = F(j).yG(iMatch);
                clr{i}(n,:)    = F(j).clr(1,:); 
            end
        end
        
        % Add nans, if no matching feet
        if n==0
            x{i}(1)      = nan;
            y{i}(1)      = nan;
            clr{i}(1,:)  = F(j).clr(1,:);
        end
    end
    
    % Figure color and positon
    fColor = 0.2.*[1 1 1];
    fPos   = [ 1 530  2002  995];
    
    % Clear 
    clear Body j idx B_ft varargin
    
    
elseif strcmp(opType,'Centroid & Rotation') || strcmp(opType,'no analysis')    
       
    % Inputs
    Body   = varargin{1}; 
    imVis  = varargin{2};
    iC     = varargin{3};
    
    % Define parameters
    frames = Body.frames;    
    numroipts = length(Body.Rotation.roi(1).xPerimL);    

    fColor = [1 1 1];
    %fPos = [1 1 1150 1340];
    fPos = [1 1 1000 700];
    
    % Offset border
        offVal = 2;
    
    %moviePath = varargin{3};

end


%% Initialize things

if ~imVis
    set(0,'DefaultFigureWindowStyle','normal')
    % Make figure
    f = figure('Position',fPos,'Color',fColor);
    set(f,'Visible','off')
else
    set(0,'DefaultFigureWindowStyle','docked')
    f = figure('Color',fColor);
end

if nargout>0
    % Initialize index
    idx = 1;
end

imRect = round([703  2868  3853  2168]);



%% Centroid only ('Centroid tracking')

if strcmp(opType,'Centroid tracking')
 
    % Loop thru data
    for i = 1:length(frames)
        
        % Get current frame
        im = getFrame(vid_path,v,frames(i),imInvert,'gray',[],iC.r);
        
        % Display frame
        imshow(im,'InitialMag','fit');
        hold on
        
        % Define on first run thru loop
        if i ==1
            totArea = size(im,1)*size(im,2);
            areaMin = totArea/10^4;
            areaMax = totArea/10^2;
        end
        
        % roi in global frame
        xG = Body.roi(i).xPerimG;
        yG = Body.roi(i).yPerimG;
        xC = Body.roi(i).xCntr;
        yC = Body.roi(i).yCntr;
        
        % Overlay blobs
        %[props,bw,areas,xB,yB] = findBlobs(im,iC.tVal,'area',areaMin,areaMax);
        [props,bw,areas,xB,yB] = findBlobs(im,iC.tVal,'coord',xC,yC);
        
        % Make a truecolor all-green image, make non-blobs invisible
        green = cat(3, zeros(size(im)), ones(size(im)), zeros(size(im)));
        h = imshow(green,'InitialMag','fit');
        set(h, 'AlphaData', bw*0.3)
        
        % Plot tracking
        h(1) = line(xG,yG,'Color','k','LineWidth',1);
        h(2) = plot(xC,yC,'k+');
        hold off
        
        % Add frame number
        text(40,100,['Frame ' num2str(frames(i))],'Color','k','FontSize',24);
        
        % Render and write frame
        drawnow
        imFrame = getframe(gca);
        writeVideo(vOut,imFrame);
        
        if imVis
            pause(0.0001);
        else
            disp(['aniData (' opType ') : ' num2str(i) ' of ' num2str(length(frames))])
        end
        
        clear h xB yB areas bw props
    end
end


%% Centroid and rotation movie ('Centroid & Rotation')

if strcmp(opType,'Centroid & Rotation')
    
    % Loop thru data
    for i = 1:length(frames)
        
        % Current whole frame
        im = getFrame(vid_path,v,frames(i),imInvert,'gray',[],iC.r);
        
        
        set(f,'Position',[1 1 500 500]);
        
        % Plot frame
%         subplot(1,2,1)
%         h = imshow(im,'InitialMag','fit');
%         hold on
        
        % Get current roi
        roi = Body.Rotation.roi(i);

        % roi in global frame
%         xG = roi.xPerimG;
%         yG = roi.yPerimG;
        
        % Body center
%         xC = Body.xCntr(i);
%         yC = Body.yCntr(i);
%         
%         % Other point (to show rotation)
%         xO = Body.xOther(i);
%         yO = Body.yOther(i);
%         
        % Local FOR border
        xL = [roi.xPerimL(1)+offVal roi.xPerimL(2)-offVal ...
            roi.xPerimL(3)-offVal roi.xPerimL(4)+offVal ...
            roi.xPerimL(5)+offVal];
        yL = [roi.yPerimL(1)+offVal roi.yPerimL(2)+offVal ...
            roi.yPerimL(3)-offVal roi.yPerimL(4)-offVal ...
            roi.yPerimL(5)+offVal];
        
        % Stabilized image
        imStable =  giveROI('stabilized',im,roi,0,Body.Rotation.tform(i));
        
        % Plot tracking
%         h(1) = line(xG,yG,'Color',[1 0 0 0.2],'LineWidth',3);
%         h(1) = line([xO xC],[yO yC],'Color',[1 0 0 0.2],'LineWidth',3);
        %h(2) = plot(xC,yC,'r+');
        
        % Plot roi
%         subplot(1,2,2)
        imshow(imStable,'InitialMag','fit')
        hold on
%         line(xL,yL,'Color',[1 0 0 0.2],'LineWidth',4);
        title(['Frame ' num2str(frames(i))])
        
        % Add frame number
        text(40,40,['Frame ' num2str(frames(i))],'Color','k','FontSize',18);
        
        % Render and write frame
        drawnow
        imFrame = getframe(gca);
        writeVideo(vOut,imFrame);
        hold off
        
        if imVis
            pause(0.001);
        else
            disp(['aniData (' opType ') : ' num2str(i) ' of ' num2str(length(frames))])
        end
        
        clear h im imStable roi xC yC x0 y0 xL yL 
    end
end


%% Movie showing the codes tracking of individual feet 
% ('Individual feet, local')

if strcmp(opType,'Individual feet, local')
    
    % Level of alpha transparency
    aLevel = 0.3;
    
    % Offset border
    offVal = 2;

    % Loop thru data
    for i = 1:length(frames)
        
        % Get current frame
        im = getFrame(vid_path,v,frames(i),imInvert,'gray',[],iC.r);
        
        % Stabilized image
        imStable =  giveROI('stabilized',im,roi(i),0,tform(i),[],0);
        
        % Display frame
        h = imshow(imStable,'InitialMag','fit');
        
        % Load B_ft data for current frame
        load([data_path filesep 'foot_blobs' filesep ...
            aData(iFile(i)).name])
        
        % Figure color
        set(f,'Color',0.2.*[1 1 1])
        
        % Plot positions of feet (after post-processing)
        for j = 1:length(x{i})
            h = scatter(x{i}(j),y{i}(j),...
                'MarkerEdgeColor','r','SizeData',1000,...
                'MarkerEdgeAlpha',0.5,'LineWidth',4);
        end
        
        % Create black local image
        bw_im = imStable==300;
        
        % Identify all blobs for candidate tube feet
        for j = 1:length(B_ft.propsL)
            % Put white pixels where blobs exist
            bw_im(B_ft.propsL(j).PixelIdxList) = 1;
        end
        
        % Make a truecolor all-green image, make non-blobs invisible
        green = cat(3, ones(size(imStable)), ones(size(imStable)), ...
            zeros(size(imStable)));
        h = imshow(green,'InitialMag','fit');
        set(h, 'AlphaData', bw_im.*aLevel)
        
        % Place frame number
        hold on
        text(20,30,['Frame ' num2str(frames(i))],'Color','w','FontSize',24);
        
        % Render and write frame
        drawnow
        imFrame = getframe(gca);
        writeVideo(vOut,imFrame);
        
        if imVis
            pause(0.001);
        else
            disp(['aniData (' opType ') : ' num2str(i) ' of ' num2str(length(frames))])
        end
        
        hold off
        
        clear h imStable B_ft green hTitle imFrame
    end
    
    clear aLevel offVal i
end


%% loop thru frames (not 'blobs G&L')

if ~strcmp(opType,'blobs G&L') && ~strcmp(opType,'Centroid & Rotation') && ...
    ~strcmp(opType,'Individual feet, local') && ~strcmp(opType,'Centroid tracking')


    % Loop thru data
    for i = 1:5 %length(frames)

        if strcmp(opType,'no analysis') 
            im = getFrame(vid_path,v,frames(i),0,'color');
            
        elseif strcmp(opType,'Individual feet, pretty')
            
            imRange = 0.002;
            
            im = getFrame(vid_path,v,frames(i),0,'color',[],iC.r);
            
            im(:,:,1) = adapthisteq(im(:,:,1),'ClipLimit',imRange);
            im(:,:,2) = adapthisteq(im(:,:,2),'ClipLimit',imRange);
            im(:,:,3) = adapthisteq(im(:,:,3),'ClipLimit',imRange);
            
            set(f,'Color',0.2.*[1 1 1])

            h = text(round(size(im,2)/7)+20,round(size(im,2)/6),...
                hTitle.String,'Color',0.8.*[1 1 1],'FontSize',18);
      
            % Plot positions of feet (after post-processing)
            for j = 1:length(x{i})
                h = scatter(x{i}(j),y{i}(j),...
                    'MarkerEdgeColor','y','SizeData',150,...
                    'MarkerEdgeAlpha',0.5,'LineWidth',1);
            end
            
            h = scatter(xArm(i,:),yArm(i,:),...
                    'MarkerEdgeColor','r','MarkerFaceColor','r',...
                    'SizeData',20,'MarkerEdgeAlpha',0.8);
           h = scatter(xCntr(i,:),yCntr(i,:),...
                    'MarkerEdgeColor','r','MarkerFaceColor','r',...
                    'SizeData',20,'MarkerEdgeAlpha',0.8);
                ttt = 3;
                
           imFrame = getframe(gca);
            
            
        else

            % Current whole frame
%             im = getFrame(vid_path,v,frames(i),imInvert,'gray');
        end
        
        if strcmp(opType,'blobs L simple')
            % Roi image, mean image subtracted
            im = giveROI('stabilized',im,S.roi(i),dSample, ...
                S.tform(:,:,i),imRoiMean);
        end
        

    
        if strcmp(opType,'Feet')
            
            subplot(3,1,1)
            
        end
        
        
        
        if strcmp(opType,'blobs G&L')
            
            if 1%~isnan(B(i).propsL)
                
                % Start with blank
                bwG  = logical(zeros(size(im)));
                
                % Score pixels with blobs
                for k = 1:length(B(i).propsG),
                    bwG(B(i).propsG(k).PixelIdxList) = 1;
                end
                
                % Make a truecolor all-green image, make non-blobs invisible
                green = cat(3, zeros(size(im)), ones(size(im)), zeros(size(im)));
                h = imshow(green,'InitialMag','fit');
                set(h, 'AlphaData', bwG.*alphaLevel)
            end
            
        elseif strcmp(opType,'blobs L simple')
            
            % Overlay blobs
            [props,bw,areas,xB,yB] = findBlobs(im,tVal,'area',areaMin,areaMax);
            
            % Make a truecolor all-green image, make non-blobs invisible
            green = cat(3, zeros(size(im)), ones(size(im)), zeros(size(im)));
            h = imshow(green,'InitialMag','fit');
            set(h, 'AlphaData', bw)
            
            clear props bw areas xB yB
         
         
        elseif strcmp(opType,'Feet')
            
            set(f,'Color',0.2.*[1 1 1])
            set(hTitle,'Color',0.8.*[1 1 1]);
            
             % Offset border
            offVal = 2;
            
            % roi in global frame
            xG = roi(i).xPerimG;
            yG = roi(i).yPerimG;
            
            % Body center
            xC = xCntr(i);
            yC = yCntr(i);
            
            % Local FOR border
            xL = [roi(i).xPerimL(1)+offVal roi(i).xPerimL(2)-offVal ...
                roi(i).xPerimL(3)-offVal roi(i).xPerimL(4)+offVal ...
                roi(i).xPerimL(5)+offVal];
            yL = [roi(i).yPerimL(1)+offVal roi(i).yPerimL(2)+offVal ...
                roi(i).yPerimL(3)-offVal roi(i).yPerimL(4)-offVal ...
                roi(i).yPerimL(5)+offVal];
            
            % Stabilized image
            imStable =  giveROI('stabilized',im,roi(i),0,tform(i),0);
            
            % Plot tracking
            h(1) = line(xG,yG,'Color',[1 1 0 0.5],'LineWidth',3);
            
            % Plot roi
            subplot(3,1,2:3)
            imshow(imStable,'InitialMag','fit')
            hold on
            line(xL,yL,'Color',[1 1 0 0.5],'LineWidth',6);
            
            % Plot centers of tube feet
            for j = 1:length(propsL{i})
                h = scatter(propsL{i}(j).Centroid(1),propsL{i}(j).Centroid(2),...
                    'MarkerEdgeColor',[1 1 0],'SizeData',300,...
                    'MarkerEdgeAlpha',0.5);
            end
       
            
        elseif strcmp(opType,'Global feet')
            
            set(f,'Color',0.2.*[1 1 1])
            set(hTitle,'Color',0.8.*[1 1 1]);
            
             % Offset border
            offVal = 2;
            
            % roi in global frame
            xG = roi(i).xPerimG;
            yG = roi(i).yPerimG;
            
            % Plot centers of tube feet
            for j = 1:length(propsG{i})
                h = scatter(propsG{i}(j).Centroid(1),propsG{i}(j).Centroid(2),...
                    'MarkerEdgeColor',[1 1 0],'SizeData',200,...
                    'MarkerEdgeAlpha',0.2);
            end
    
        elseif strcmp(opType,'Local feet')
            
            set(f,'Color',0.2.*[1 1 1])
            set(hTitle,'Color',0.8.*[1 1 1]);
            imFrame = getframe(gca);
            
        end
        
        % Run get frame
        if strcmp(opType,'Centroid & Rotation')
            drawnow
            imFrame = getframe(f);
        else
            hT = text(20,30,['Frame ' num2str(frames(i))],'Color','w','FontSize',24);
            drawnow
            imFrame = getframe(gca);
        end
        
%         pause(0.001)  
 
%         imRect = [314  1300-974 1744 974];

        
%         imFrame.cdata = imcrop(imFrame.cdata,imRect);
        writeVideo(vOut,imFrame);
        
        if imVis
            pause(0.001);
        else
            disp(['aniData (' opType ') : ' num2str(i) ' of ' num2str(length(frames))])
        end
        
        hold off
        
        clear imRange im h fColor fPos imFrame h hT
    end
end


%% loop thru frames ('blobs G&L')

if strcmp(opType,'blobs G&L')
    
    % Loop thru data
    for i = 1:length(B)
        
        if isfield(B(i).propsG,'Area')
            % Current whole frame
            im = getFrame(vid_path,v,B(i).fr_num,imInvert,'gray');
            
            % Display frame
            h = imshow(imcompliment(im),'InitialMag','fit');
            hold on

            % Start with blank
            bwG  = logical(zeros(size(im)));
            
            % Score pixels with blobs
            for k = 1:length(B(i).propsG),
                bwG(B(i).propsG(k).PixelIdxList) = 1;
            end
            
            % Make a truecolor all-green image, make non-blobs invisible
            green = cat(3, zeros(size(im)), ones(size(im)), zeros(size(im)));
            h = imshow(green,'InitialMag','fit');
            set(h, 'AlphaData', bwG.*alphaLevel)
            
            title(['Frame ' num2str(B(i).fr_num)]);
            
            if nargout>0
                % Capture frame
                M(idx) = getframe(gcf);
                
                % Advance index
                idx = idx + 1;
            end
            
            if imVis
                pause(0.001);
                hold off
                delete(h)
            else
                disp(['aniData (' opType ') : ' num2str(i) ' of ' num2str(length(frames))])
            end
            
            
        end
    end
end


%% Finish up

close(f)

% Output
%varargout{1} = M;
varargout{1} = [];

close(vOut)

set(0,'DefaultFigureWindowStyle','docked')



