function makeDataMovie(vid_path,v,S,mov_path,imVis,movType,varargin)


%% Parameters

% Video quality
vidQuality = 50;

% Invert image
imInvert = 0;

% Downsample roi
dSample = 0;

% Number of point for roi circle
numroipts = 300;

% Marker size in roi
mSize(1) = 2000;

% Marker size in whole frame
mSize(2) = 300;

% Color 
clrs{1} = [241 90 36]./255;

% Color around roi
clrs{2} = [41 171 226]./255;

% Color around tube feet
clrs{3} = [0 1 0];

% Number of frames for the 'two view skip' mode
num_frames = 50;


%% Intializing things


% Set up output video file
vOut = VideoWriter([mov_path '.mp4'],'MPEG-4');
vOut.Quality = 50;
open(vOut)

% Make figure
set(0,'DefaultFigureWindowStyle','normal')
f = figure('Visible',imVis);
fPos = get(f,'Position');
set(f,'Position',[fPos(1) fPos(2) 1020 1271])



%% Creat output

% Define frames
if strcmp(movType,'two view skip')
    frames = round(linspace(S.frames(1),S.frames(end),num_frames));
else
    frames = S.frames;
end

% Loop thru data
for i = 1:length(frames)
    
    % Index for frame in the data
    iFrame = find(S.frames==frames(i),1,'first');
    
    % Current whole frame
    im = getFrame(vid_path,v,S.frames(iFrame),imInvert,'rgb');
       
    if strcmp(movType,'two view') || strcmp(movType,'two view skip')
        
        % Coordinate origin
        %cOrigin = [S.xCntr(iFrame); S.yCntr(iFrame)];
        
        % Get roi image
        [im_roi,bw_mask] = giveROI('stabilized',im,S.roi(iFrame),...
            dSample,S.tform(iFrame));
        
        
        % Plot image
        %figure(f)
        
        % Top image
        subplot(2,1,1)
        h = imshow(im,'InitialMag','fit');
        %title(['Frame ' num2str(S.frames(iFrame))])
        a1 = gca;
        hold on      
        
        % Outline roi (global)
        line(S.roi(iFrame).xPerimG,S.roi(iFrame).yPerimG,'Color',...
            [clrs{2} 0.5],'LineWidth',1); 

        % BOTTOM IMAGE  ---------------------
        subplot(2,1,2)
        
        h = imshow(im_roi,'InitialMag','fit');
        a2 = gca;
        hold on
        
         % Set position of both images
        set(a1,'Position',[0 0.27 1 1])
        
        % Add frame number
        %axes(a2); 
        h = text(1,30,['Frame ' num2str(S.frames(iFrame))],...
                      'Color','k','FontSize',22,'FontWeight','bold','Parent',a2);
        
        %set(a2,'Position',[0.2 -.05 0.65 0.65])
        set(a2,'Position',[0.25 0 0.5 0.5])
        
        %axes(a2)
        
        % Loop thru tube feet 
        for j = 1:length(S.ft)
            
            % Top plot
            line([S.ft(j).xBase(iFrame) S.ft(j).xTip(iFrame)],...
                 [S.ft(j).yBase(iFrame) S.ft(j).yTip(iFrame)],...
                 'Color',[clrs{3} 0.3],'LineWidth',5,'Parent',a1);
            scatter(S.ft(j).xTip(iFrame),S.ft(j).yTip(iFrame),...
                    'MarkerFaceColor',clrs{3},...
                    'MarkerEdgeColor',clrs{3},'LineWidth',2,...
                    'MarkerEdgeAlpha',0.8,'MarkerFaceAlpha',0.8,...
                    'SizeData',100,'Parent',a1); 
               
            % Bottom plot
             
            % Tip plot in roi
            tipR = [S.ft(j).xTipL(iFrame) + size(im_roi,1)/2 ...
                    S.ft(j).yTipL(iFrame) + size(im_roi,2)/2];
                
           % Highlight foot
           scatter(tipR(1),tipR(2),'MarkerEdgeColor',clrs{3},...
               'LineWidth',2,'MarkerEdgeAlpha',0.5,'SizeData',mSize(1),...
               'Parent',a2);     
        end
        
%         scatter(tipROI(1,:),tipROI(2,:),'MarkerEdgeColor',clrs{3},...
%             'LineWidth',2,'MarkerEdgeAlpha',0.5,'SizeData',mSize(1));
        set(a1,'NextPlot','replace')
        set(a2,'NextPlot','replace')

        clear tipG tipROI
    end
    
    % Get image
    imFrame = getframe(gcf);
    
    % Write frame
    writeVideo(vOut,imFrame);
    
    if ~imVis
       disp(['DataMovie: Done ' num2str(i) ' of ' num2str(length(frames)) ' frames'])
    end
end

close(f)

% Output
%varargout{1} = M;
%varargout{1} = [];

close(vOut)


