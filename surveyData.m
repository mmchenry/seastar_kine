function surveyData(vid_path,v,imInvert,opType,varargin)
% Based on aniData, it runs thru a movie and generates a bunch of movie
% stills to validate that an anlysis has worked well


%% Parameters

alphaLevel = 0.5;


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
    
    % Number of rows and columns in each fig
    nRow = 2;
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
    
    frames = Body.frames;    
    numroipts = length(Body.Rotation.roi(1).xPerimL);    

    % Number of rows and columns in each fig
    nRow = 3;
    nCol = 2;
    
    % Number of panels for each frame
    pNumAdvance = 2;
    
    %moviePath = varargin{3};
    
    % Figure color
    fColor = [1 1 1];

    
elseif strcmp(opType,'Feet')
    
    Body     = varargin{1}; 
    B_ft     = varargin{2};
    numVis   = varargin{3};
    
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
           
           % Store from blob data
           propsG{j} = B_ft(i).propsG;
           propsL{j} = B_ft(i).propsL;
           
           % Store from Body
           roi(j,1)     = Body.Rotation.roi(idx);
           xCntr(j,1)   = Body.xCntr(idx);
           yCntr(j,1)   = Body.yCntr(idx);
           tform(j,1)   = Body.Rotation.tform(idx);
           
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
    
    Body     = varargin{1}; 
    B_ft     = varargin{2};
    numVis   = varargin{3};
    
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
           
           % Store from blob data
           propsG{j} = B_ft(i).propsG;
           propsL{j} = B_ft(i).propsL;
           
           % Store from Body
           roi(j,1)     = Body.Rotation.roi(idx);
           xCntr(j,1)   = Body.xCntr(idx);
           yCntr(j,1)   = Body.yCntr(idx);
           tform(j,1)   = Body.Rotation.tform(idx);
           
           % Advance index
           j = j + 1;
       end
       
       
    end
      
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

set(0,'DefaultFigureWindowStyle','docked')
%set(0,'DefaultFigureWindowStyle','normal')

% Make figure
%f = figure('Position',[1 1 1150 1340]);
f = figure('Color',fColor);

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
        im = getFrame(vid_path,v,anaFrames(i),imInvert,'gray');
        
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
            
            subplot(2,2,pNum)
            
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
            
            hTitle = title(['Frame ' num2str(anaFrames(i))]);
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
            green = cat(3, zeros(size(im)), ones(size(im)), zeros(size(im)));
            h = imshow(green,'InitialMag','fit');
            set(h, 'AlphaData', bw*0.3)
            
            % Plot tracking
            h(1) = line(xG,yG,'Color','k','LineWidth',1);
            h(2) = plot(xC,yC,'k+');
            
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
            imStable =  giveROI('stabilized',im,roi,0,Body.Rotation.tform(iData));
            
            % Plot tracking
            h(1) = line(xG,yG,'Color',[1 0 0 0.2],'LineWidth',3);
            h(1) = line([xO xC],[yO yC],'Color',[1 0 0 0.2],'LineWidth',3);
            %h(2) = plot(xC,yC,'r+');
            
            % Plot roi
            subplot(nRow,nCol,pNum+1)
            imshow(imStable,'InitialMag','fit')
            hold on
            line(xL,yL,'Color',[1 0 0 0.2],'LineWidth',4);
            
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
            imStable =  giveROI('stabilized',im,roi(iData),0,tform(iData),0);
            
            % Plot tracking
            h(1) = line(xG,yG,'Color',[1 1 0 0.5],'LineWidth',3);
            
            % Plot roi
            subplot(nRow,nCol,pNum+1)
            imshow(imStable,'InitialMag','fit')
            hold on
            line(xL,yL,'Color',[1 1 0 0.5],'LineWidth',6);
            
            % Plot centers of tube feet
            for j = 1:length(propsL{iData})
                h = scatter(propsL{iData}(j).Centroid(1),propsL{iData}(j).Centroid(2),...
                    'MarkerEdgeColor',[1 1 0],'SizeData',300,...
                    'MarkerEdgeAlpha',0.5);
            end
            
            drawnow
            pause(0.001)
        
            
        elseif strcmp(opType,'Feet local')

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
            imStable =  giveROI('stabilized',im,roi(iData),0,tform(iData),0);
            
            % Plot roi
            imshow(imStable,'InitialMag','fit')
            hold on
            %line(xL,yL,'Color',[1 1 0 0.5],'LineWidth',6);
            
            % Plot centers of tube feet
            for j = 1:length(propsL{iData})
                h = scatter(propsL{iData}(j).Centroid(1),propsL{iData}(j).Centroid(2),...
                    'MarkerEdgeColor',[1 1 0],'SizeData',300,...
                    'MarkerEdgeAlpha',0.5);
            end
            
            hTitle = title(['Frame ' num2str(anaFrames(i))]);
            set(f,'Color',0.2.*[1 1 1])
            set(hTitle,'Color',0.8.*[1 1 1]);
            
            drawnow
                       
        end
        
        hold off
        
        % If panel number equal to the total number . . .
        if ((pNum+pNumAdvance)>nCol*nRow) && (i<length(anaFrames))
            
            % New figure window
            f = figure('Color',fColor);
            
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

f = figure;

if strcmp(opType,'Centroid & Rotation')
    
        % Whole frames
        im1 = getFrame(vid_path,v,anaFrames(1),imInvert,'gray');
        im2 = getFrame(vid_path,v,anaFrames(end),imInvert,'gray');
        
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


