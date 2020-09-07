function surveyData(vid_path,v,imInvert,opType,varargin)
% Based on aniData, it runs thru a movie and generates a bunch of movie
% stills to validate that an anlysis has worked well


%% Parameters

alphaLevel = 0.5;

% Default title text
tText = '';


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
    
    iC      = varargin{2};
    numVis  = varargin{3};
    
    if length(varargin)>3
        tText = varargin{4};
    end
    
    % Number of rows and columns in each fig
    nRow = 5;
    nCol = 2;
    
    % Number of panels for each frame
    pNumAdvance = 1;
    
    % Figure color
    fColor = [1 1 1];

elseif strcmp(opType,'Pred + prey')
    
    % Parse inputs
    S_py   = varargin{1};
    S_pd   = varargin{2};
    L      = varargin{3};
    imVis  = varargin{4};

    numroipts = 400;
    frames = S_pd.frames;    
    
elseif strcmp(opType,'Pred-prey cent track')
    
    % Parse inputs
    C_pd   = varargin{1};
    C_py   = varargin{2};
    r_pd   = varargin{3};
    r_py   = varargin{4};
    imVis  = varargin{5};

    numroipts = 400;
    frames = C_pd.frames;    
    
elseif strcmp(opType,'Centroid & Rotation') || strcmp(opType,'no analysis')    
       
    Body    = varargin{1}; 
    iC      = varargin{2};
    numVis  = varargin{3};
    
    if length(varargin)>3
        tText = varargin{4};
    end
    
    frames = Body.frames;    
    numroipts = length(Body.Rotation.roi(1).xPerimL);    

    % Number of rows and columns in each fig
    nRow = 4;
    nCol = 2;
    
    % Number of panels for each frame
    pNumAdvance = 2;
    
    %moviePath = varargin{3};
    
    % Figure color
    fColor = [1 1 1];

    
elseif strcmp(opType,'Feet')
    
    currDataPath    = varargin{1}; 
    Body            = varargin{2}; 
%     B_ft     = varargin{2};
    numVis          = varargin{3};
    iC              = varargin{4};
    tText           = varargin{5};
    
    % Path to foot data (B_ft files)
    footPath = [currDataPath filesep 'foot_blobs'];
    
    % List files of B_ft
    [aB_ft,frame_list] = listB_ft(footPath);
    
    % Start index
    j = 1;
    
    % Loop thru frames
    for i = 1:length(aB_ft)
        
        % load B_ft
        load([footPath filesep aB_ft(i).name])
        
       % If there is foot blob data . . .
       if ~isempty(B_ft.frIdx)           
           
           % Index in the Body data
           idx = find(Body.frames==frame_list(i),1,'first');
           
           % Store from blob data
           propsG{j} = B_ft.propsG;
           propsL{j} = B_ft.propsL;
           
           % Store from Body
           roi(j,1)     = Body.Rotation.roi(idx);
           xCntr(j,1)   = Body.xCntr(idx);
           yCntr(j,1)   = Body.yCntr(idx);
           tform(j,1)   = Body.Rotation.tform(idx);
           frames(j,1)  = frame_list(i);
           
           % Advance index
           j = j + 1;
       end
    end
      
    numroipts = length(Body.Rotation.roi(1).xPerimL);    

    % Clear 
    clear Body j idx B_ft
    
    % Number of rows and columns in each fig
    nRow = 3;
    nCol = 2;
    
    % Number of panels for each frame
    pNumAdvance = 2;
    
    % Figure color
    fColor = 0.2.*[1 1 1];   
    
elseif strcmp(opType,'Feet local')
    
    currDataPath = varargin{1};  
    Body         = varargin{2}; 
    numVis       = varargin{3};
    iC           = varargin{4};
    F            = varargin{5};
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

end


%% Initialize things

% Frames to analyze
anaFrames = round(linspace(min(frames),max(frames),numVis));

% anaFrames = round(linspace(round(max(frames)/2),max(frames),numVis));

set(0,'DefaultFigureWindowStyle','docked')
%set(0,'DefaultFigureWindowStyle','normal')

% Make figure
%f = figure('Position',[1 1 1150 1340]);
f = figure('Color',fColor,'Name',tText,'NumberTitle','off');

if nargout>0
    % Initialize index
    idx = 1;
end

% Current panel number
pNum = 1;


%% loop thru frames (not 'blobs G&L')

if ~strcmp(opType,'blobs G&L')
    
    % Loop thru data
    for i = 1:length(anaFrames)
        
        % Index in data for current frame
        iData = find(frames==anaFrames(i),1,'first');
        
        % Designate panel
        figure(f)
        subplot(nRow,nCol,pNum)

        % Current whole frame
        im = getFrame(vid_path,v,anaFrames(i),imInvert,'gray',[],iC.r);
        
        if strcmp(opType,'blobs L simple')
            % Roi image, mean image subtracted
            im = giveROI('stabilized',im,S.roi(iData),dSample, ...
                S.tform(:,:,iData),imRoiMean);
        end
        
        if strcmp(opType,'Centroid & Rotation')
            
            % Get current roi
            roi = Body.Rotation.roi(iData);
    
        elseif strcmp(opType,'Feet')
            
            %subplot(3,1,1)
            
        elseif strcmp(opType,'Centroid tracking')
            
            subplot(nRow,nCol,pNum)
            
            %roi = Body.roi(iData);
            if i ==1
                totArea = size(im,1)*size(im,2);
                areaMin = totArea/10^4;
                areaMax = totArea/10^2;
            end

        elseif strcmp(opType,'Pred + prey')
            subplot(4,2,[1:4])
        end
        
        if ~strcmp(opType,'Feet local')
            
            % Display frame
            h = imshow(im,'InitialMag','fit');
            hold on
            
            hTitle = title([tText ': Frame ' num2str(anaFrames(i))]);
        end
        
        if strcmp(opType,'blobs G&L')
            
            if 1%~isnan(B(iData).propsL)
                
                % Start with blank
                bwG  = logical(zeros(size(im)));
                
                % Score pixels with blobs
                for k = 1:length(B(iData).propsG),
                    bwG(B(iData).propsG(k).PixelIdxList) = 1;
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
         
            
        elseif strcmp(opType,'Centroid tracking')
            
            % roi in global frame
            xG = Body.roi(iData).xPerimG;
            yG = Body.roi(iData).yPerimG;
            xC = Body.roi(iData).xCntr;
            yC = Body.roi(iData).yCntr;
            
            % Overlay blobs
            %[props,bw,areas,xB,yB] = findBlobs(im,iC.tVal,'area',areaMin,areaMax);
            [props,bw,areas,xB,yB] = findBlobs(im,iC.tVal,'coord',xC,yC);
            
            % Make a truecolor all-green image, make non-blobs invisible
%             green = cat(3, zeros(size(im)), ones(size(im)), zeros(size(im)));
%             h = imshow(green,'InitialMag','fit');
%             set(h, 'AlphaData', bw*0.3)
            
            % Plot tracking
            h(1) = line(xG,yG,'Color','g','LineWidth',1);
            h(2) = plot(xC,yC,'g+');
            
            clear xB yB areas bw props
            
            
        elseif strcmp(opType,'Centroid & Rotation')
            
            % Offset border
            offVal = 2;
            
            % roi in global frame
            xG = roi.xPerimG;
            yG = roi.yPerimG;
            
            % Body center
            xC = Body.xCntr(iData);
            yC = Body.yCntr(iData);
            
            % Other point (to show rotation)
            xO = Body.xOther(iData);
            yO = Body.yOther(iData);
            
            % Local FOR border
            xL = [roi.xPerimL(1)+offVal roi.xPerimL(2)-offVal ...
                roi.xPerimL(3)-offVal roi.xPerimL(4)+offVal ...
                roi.xPerimL(5)+offVal];
            yL = [roi.yPerimL(1)+offVal roi.yPerimL(2)+offVal ...
                roi.yPerimL(3)-offVal roi.yPerimL(4)-offVal ...
                roi.yPerimL(5)+offVal];
            
            % Stabilized image
            imStable =  giveROI('stabilized',im,roi,0,...
                Body.Rotation.tform(iData));
            
            % Plot tracking
            h(1) = line(xG,yG,'Color',[0 1 0 0.2],'LineWidth',3);
            h(1) = line([xO xC],[yO yC],'Color',[0 1 0 0.2],'LineWidth',3);
            %h(2) = plot(xC,yC,'r+');
            
            % Plot roi
            subplot(nRow,nCol,pNum+1)
            imshow(imStable,'InitialMag','fit')
            hold on
            line(xL,yL,'Color',[0 1 0 0.2],'LineWidth',4);
            
            drawnow
            pause(0.001)
            
            
        elseif strcmp(opType,'Feet')
            
            set(f,'Color',0.2.*[1 1 1])
            set(hTitle,'Color',0.8.*[1 1 1]);
            
             % Offset border
            offVal = 2;
            
            % roi in global frame
            xG = roi(iData).xPerimG;
            yG = roi(iData).yPerimG;
            
            % Body center
            xC = xCntr(iData);
            yC = yCntr(iData);
            
            % Local FOR border
            xL = [roi(iData).xPerimL(1)+offVal roi(iData).xPerimL(2)-offVal ...
                roi(iData).xPerimL(3)-offVal roi(iData).xPerimL(4)+offVal ...
                roi(iData).xPerimL(5)+offVal];
            yL = [roi(iData).yPerimL(1)+offVal roi(iData).yPerimL(2)+offVal ...
                roi(iData).yPerimL(3)-offVal roi(iData).yPerimL(4)-offVal ...
                roi(iData).yPerimL(5)+offVal];
            
            % Stabilized image
            imStable =  giveROI('stabilized',im,roi(iData),0,tform(iData),[],0);
            
            % Plot tracking
            h(1) = line(xG,yG,'Color',[1 0 0 0.5],'LineWidth',3);
            plot(xC,yC,'g+')
            
            % Plot centers of tube feet
            for j = 1:length(propsG{iData})
                xFG(j,1) = propsG{iData}(j).Centroid(1);
                yFG(j,1) = propsG{iData}(j).Centroid(2);
            end
            
            h = scatter(xFG,yFG,'MarkerEdgeColor',[1 0 0],'SizeData',10,...
                    'MarkerEdgeAlpha',1);
            
            % Plot roi
            subplot(nRow,nCol,pNum+1)
            imshow(imStable,'InitialMag','fit')
            hold on
            line(xL,yL,'Color',[1 0 0 0.5],'LineWidth',6);
            
            % Plot centers of tube feet
            for j = 1:length(propsL{iData})
                h = scatter(propsL{iData}(j).Centroid(1),propsL{iData}(j).Centroid(2),...
                    'MarkerEdgeColor',[1 0 0],'SizeData',150,...
                    'MarkerEdgeAlpha',0.5);
            end
            
            drawnow
            pause(0.001)
        
            
        elseif strcmp(opType,'Feet local')
            
            % Level of alpha transparency
            aLevel = 0.5;
            
             % Offset border
            offVal = 2;
            
            % Load B_ft data for current frame
            load([currDataPath filesep 'foot_blobs' filesep ...
                  aData(iData).name])
           
            % roi in global frame
            xG = roi(iData).xPerimG;
            yG = roi(iData).yPerimG;
            
            % Body center
            xC = xCntr(iData);
            yC = yCntr(iData);
            
            % Local FOR border
            xL = [roi(iData).xPerimL(1)+offVal roi(iData).xPerimL(2)-offVal ...
                roi(iData).xPerimL(3)-offVal roi(iData).xPerimL(4)+offVal ...
                roi(iData).xPerimL(5)+offVal];
            yL = [roi(iData).yPerimL(1)+offVal roi(iData).yPerimL(2)+offVal ...
                roi(iData).yPerimL(3)-offVal roi(iData).yPerimL(4)-offVal ...
                roi(iData).yPerimL(5)+offVal];
            
            % Stabilized image
            imStable =  giveROI('stabilized',im,roi(iData),0,tform(iData),[],0);
            
            % Plot roi
            imshow(imStable,'InitialMag','fit')
            hold on
            %line(xL,yL,'Color',[1 1 0 0.5],'LineWidth',6);
            
            % Create black local image
            bw_im = imStable==300;
            
            % Identify all blobs for candidate tube feet
            for j = 1:length(B_ft.propsL)
                % Put white pixels where blobs exist
                bw_im(B_ft.propsL(j).PixelIdxList) = 1;
            end
            
             % Make a truecolor all-green image, make non-blobs invisible
            green = cat(3, zeros(size(imStable)), ones(size(imStable)), ...
                           zeros(size(imStable)));
            h = imshow(green,'InitialMag','fit');
            set(h, 'AlphaData', bw_im.*aLevel)
            
%             % Plot centers of tube feet
%             for j = 1:length(B_ft.propsL)
%                 h = scatter(B_ft.propsL(j).Centroid(1),B_ft.propsL(j).Centroid(2),...
%                     'MarkerEdgeColor',[1 1 0],'SizeData',150,...
%                     'MarkerEdgeAlpha',0.5);
%             end

             % Plot centers of tube feet
            for j = 1:length(xFt{iData})
                h = scatter(xFt{iData}(j),yFt{iData}(j),...
                    'MarkerEdgeColor',[1 1 0],'SizeData',150,...
                    'MarkerEdgeAlpha',0.5);
            end
            
            hTitle = title([tText ': Frame ' num2str(anaFrames(i))]);
            set(f,'Color',0.2.*[1 1 1])
            set(hTitle,'Color',0.8.*[1 1 1]);
            
            drawnow
                       
        end
        
        hold off
        
        % If panel number equal to the total number . . .
        if ((pNum+pNumAdvance)>nCol*nRow) && (i<length(anaFrames))
            
            % New figure window
            f = figure('Color',fColor,'Name',tText,'NumberTitle','off');
            
            % Reset panel number
            pNum = 1;
            
        % Otherwise . . .
        else
                     
            % Advance panel number
            pNum = pNum + pNumAdvance;

        end
    end
end


%% Additional image

if strcmp(opType,'Centroid & Rotation')
    
    f = figure;
    
    % Whole frames
    im1 = getFrame(vid_path,v,anaFrames(1),imInvert,'gray',[],iC.r);
    im2 = getFrame(vid_path,v,anaFrames(end),imInvert,'gray',[],iC.r);
    
    % Index in data for current frame
    iData1 = find(Body.frames==anaFrames(1),1,'first');
    iData2 = find(Body.frames==anaFrames(end),1,'first');
    
    % Get current roi
    roi1 = Body.Rotation.roi(iData1);
    roi2 = Body.Rotation.roi(iData2);
    
    % Stabilized images
    imStable1 =  giveROI('stabilized',im1,roi1,0,Body.Rotation.tform(iData1));
    imStable2 =  giveROI('stabilized',im2,roi2,0,Body.Rotation.tform(iData2));
    
    warning off
    imshowpair(imStable1,imStable2)
    warning on
    title(['Frames ' num2str(anaFrames(1)) ' and ' ...
        num2str(anaFrames(end)) ' compared'])
end


%% loop thru frames ('blobs G&L')

if strcmp(opType,'blobs G&L')
    
    f = figure;
    
    % Loop thru data
    for i = 1:length(B)
        
        if isfield(B(i).propsG,'Area')
            % Current whole frame
            im = getFrame(vid_path,v,B(i).fr_num,imInvert,'gray',[],iC.r);
            
            % Display frame
            h = imshow(imcompliment(im),'InitialMag','fit');
            hold on

            % Start with blank
            bwG  = logical(zeros(size(im)));
            
            % Score pixels with blobs
            for k = 1:length(B(i).propsG)
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




function [a,frNums] = listB_ft(footPath)

a = dir([footPath filesep 'foot_blobs*']);

frNums = [];

% Loop thru file entries
for i = 1:length(a)

     % Index of separator
    iSep = find(a(i).name=='_',1,'last');
    
    % Get frame number
    a(i).frNum = str2num(a(i).name((iSep+1):end-4));
    
    % Listing of frame numbers
    frNums = [frNums; a(i).frNum];
end



