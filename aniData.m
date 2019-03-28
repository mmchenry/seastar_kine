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

% Set up output video file
vOut = VideoWriter([data_path filesep movie_file '.mp4'],'MPEG-4');
vOut.Quality = 50;

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
           
           % Advance index
           j = j + 1;
       end
       
       
    end
      
    numroipts = length(Body.Rotation.roi(1).xPerimL);    

    % Clear 
    clear Body j idx B_ft
    
elseif strcmp(opType,'Centroid & Rotation') || strcmp(opType,'no analysis')    
       
    Body = varargin{1}; 
    frames = Body.frames;    
    numroipts = length(Body.Rotation.roi(1).xPerimL);    
    if nargin > 1
        imVis = varargin{2};
    else
        imVis = 1;
    end
    
    %moviePath = varargin{3};

end



%% Initialize things

set(0,'DefaultFigureWindowStyle','normal')

% Make figure
f = figure('Position',[1 1 1150 1340]);

if ~imVis
    set(f,'Visible','off')
end

if nargout>0
    % Initialize index
    idx = 1;
end


%% loop thru frames (not 'blobs G&L')

if ~strcmp(opType,'blobs G&L')
    
    % Loop thru data
    for i = 1:length(frames)

        if strcmp(opType,'no analysis')
            im = getFrame(vid_path,v,frames(i),0,'color');
        else
            % Current whole frame
            im = getFrame(vid_path,v,frames(i),imInvert,'gray');
        end
        
        if strcmp(opType,'blobs L simple')
            % Roi image, mean image subtracted
            im = giveROI('stabilized',im,S.roi(i),dSample, ...
                S.tform(:,:,i),imRoiMean);
        end
        
        if strcmp(opType,'Centroid & Rotation')
            subplot(1,2,1)

            % Get current roi
            roi = Body.Rotation.roi(i);
    
        elseif strcmp(opType,'Feet')
            
            subplot(3,1,1)
            
        elseif strcmp(opType,'Centroid tracking')
 
            %roi = Body.roi(i);
            if i ==1
                totArea = size(im,1)*size(im,2);
                areaMin = totArea/10^4;
                areaMax = totArea/10^2;
            end

        elseif strcmp(opType,'Pred + prey')
            subplot(4,2,[1:4])
        end
        
        % Display frame
        h = imshow(im,'InitialMag','fit');
        hold on
        hTitle = title(['Frame ' num2str(frames(i))]);
        
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
         
            
        elseif strcmp(opType,'Centroid tracking')
            
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
            
            clear xB yB areas bw props
            

        elseif strcmp(opType,'Pred + prey')
            
            set(f,'WindowStyle','normal')
            set(f,'Position',[1 1 881 692])
            
            % Colors
            pyClr = [41 171 226]./255;
            pdClr = [241 90 36]./255;
            
            
            angPd = -L.ang_pd + S_pd.ang(i);
            angPy = -L.ang_py + S_py.ang(i);
            
            thd = linspace(0,2*pi,400);
            roiPd.x = (S_pd.roi(i).r-1).*cos(thd)+S_pd.roi(i).r+0.5;
            roiPd.y = (S_pd.roi(i).r-1).*sin(thd)+S_pd.roi(i).r+0.5;
            roiPy.x = (S_py.roi(i).r-1).*cos(thd)+S_py.roi(i).r+0.5;
            roiPy.y = (S_py.roi(i).r-1).*sin(thd)+S_py.roi(i).r+0.5;
            
            % Current rois                
            imPd = giveROI('stabilized',im,S_pd.roi(i),0,angPd);
            imPy = giveROI('stabilized',im,S_py.roi(i),0,angPy);
                 
            % roi of pred in global frame
            xGpd = S_pd.roi(i).xPerimG;
            yGpd = S_pd.roi(i).yPerimG;
            xCpd = [S_pd.roi(i).xCntr S_pd.roi(i).xCntr+S_pd.roi(i).r*cosd(-angPd)];
            yCpd = [S_pd.roi(i).yCntr S_pd.roi(i).yCntr+S_pd.roi(i).r*sind(-angPd)];
            
            % roi of prey in global frame
            xGpy = S_py.roi(i).xPerimG;
            yGpy = S_py.roi(i).yPerimG;
            xCpy = [S_py.roi(i).xCntr S_py.roi(i).xCntr+S_py.roi(i).r*cosd(-angPy)];
            yCpy = [S_py.roi(i).yCntr S_py.roi(i).yCntr+S_py.roi(i).r*sind(-angPy)];
            
            % Plot tracking
            h(1) = line(xGpd,yGpd,'Color',[pdClr 0.8],'LineWidth',2);
            h(2) = line(xCpd,yCpd,'Color',[pdClr 0.8],'LineWidth',1);
            h(3) = line(xGpy,yGpy,'Color',[pyClr 0.8],'LineWidth',1);
            h(4) = line(xCpy,yCpy,'Color',[pyClr 0.8],'LineWidth',0.5);
 
            subplot(4,2,[5 7])
            imshow(imPd,'InitialMag','fit');
            hold on
            line(roiPd.x,roiPd.y,...
                        'Color',[pdClr],'LineWidth',7);
            line([S_pd.roi(i).r 2*S_pd.roi(i).r],[S_pd.roi(i).r S_pd.roi(i).r],...
                        'Color',[pdClr 0.5],'LineWidth',2);      
            hold off
            
            subplot(4,2,[6 8])
            imshow(imPy,'InitialMag','fit');
            hold on
            line(roiPy.x,roiPy.y,...
                        'Color',[pyClr],'LineWidth',7);
            line([S_py.roi(i).r 2*S_py.roi(i).r],S_py.roi(i).r.*[1 1],...
                        'Color',[pyClr 0.5],'LineWidth',2);
            hold off
            
      
        elseif strcmp(opType,'Pred-prey cent track')
            
            % Current rois
            roiPd = giveROI('define','circular',numroipts,r_pd,...
                     C_pd.x(i),C_pd.y(i));
            roiPy = giveROI('define','circular',numroipts,r_py,...
                     C_py.x(i),C_py.y(i));
                 
            % Plot
            h(1) = line(roiPd.xPerimG,roiPd.yPerimG,'Color',[1 0 0 0.2],'LineWidth',3);
            h(2) = line(roiPy.xPerimG,roiPy.yPerimG,'Color',[0 0 1 0.2],'LineWidth',3);
            h(3) = plot(C_pd.x(i),C_pd.y(i),'r+');
            h(4) = plot(C_py.x(i),C_py.y(i),'b+');   
            
            
        elseif strcmp(opType,'Centroid & Rotation')
            
            % Offset border
            offVal = 2;
            
            % roi in global frame
            xG = roi.xPerimG;
            yG = roi.yPerimG;
            
            % Body center
            xC = Body.xCntr(i);
            yC = Body.yCntr(i);
            
            % Other point (to show rotation)
            xO = Body.xOther(i);
            yO = Body.yOther(i);
            
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
            h(1) = line(xG,yG,'Color',[1 0 0 0.2],'LineWidth',3);
            h(1) = line([xO xC],[yO yC],'Color',[1 0 0 0.2],'LineWidth',3);
            %h(2) = plot(xC,yC,'r+');
            
            % Plot roi
            subplot(1,2,2)
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
            
            drawnow
            pause(0.001)
            
                       
        end
 
        
%         if nargout>0
%             % Capture frame
%             M(idx) = getframe(gcf);
%             
%             % Advance index
%             idx = idx + 1;
%         end

        imFrame = getframe(f);
        writeVideo(vOut,imFrame);
        
        if imVis
            pause(0.001);
        else
            disp(['aniData (' opType ') : ' num2str(i) ' of ' num2str(length(frames))])
        end
        
        hold off
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


% function visTrack(im,x,y,r,theta,tform,t_txt)
% 
% % Do not downsample the roi image:
% dSample = 0;
% 
% % If rotation included . . 
% if nargin>5 && ~isempty(tform)
%     % Focus on roi
%     [im_roi,bw_mask,roi_rect,bw_roi,imStable] = giveROI('circular',im,x,y,r,theta,...
%         dSample,tform);
%   
%     subplot(1,2,2)
%     imshow(imStable,'InitialMagnification','fit');
%     
%     % Set up for main plot
%     subplot(1,2,1)
% end
% 
% % Circular coordinates for new roi
% xC    = r.*cos(theta) + x;
% yC    = r.*sin(theta) + y;
% 
% imshow(im,'InitialMagnification','fit');
% hold on
% h = line(xC,yC,'Color',[1 0 0 0.2],'LineWidth',3);
% title(t_txt)
% plot(x,y,'r+')
% hold off
% 


