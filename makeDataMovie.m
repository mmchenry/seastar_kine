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


%% Intializing things


% Set up output video file
vOut = VideoWriter([mov_path '.mp4'],'MPEG-4');
vOut.Quality = 50;
open(vOut)

% Make figure
f = figure('Visible',imVis);
fPos = get(f,'Position');
set(f,'Position',[fPos(1) fPos(2) 1020 1534])



%% Creat output


% Loop thru data
for i = 1:length(S.frames)
    
    % Current whole frame
    im = getFrame(vid_path,v,S.frames(i),imInvert,'rgb');
    
    
    if strcmp(movType,'two view')
        
        % Coordinate origin
        cOrigin = [S.xCntr(i); S.yCntr(i)];
        
        % Get roi image
        [im_roi,bw_mask] = giveROI('stabilized',im,S.roi(i),...
            dSample,S.tform(i));
        
        % Loop thru tube feet
        for j = 1:length(S.ft)
            
            % Global tip point
            tipG(:,j) = [S.ft(j).xTip(i); S.ft(j).yTip(i)];
                
        end
        
        % Local tip point
        tipL = G2L(cOrigin,-S.ang(i)/180*pi,tipG);
        
        % Tip point in region of interest
        tipROI(1,:) = tipL(1,:) + size(im_roi,1)/2;
        tipROI(2,:) = tipL(2,:) + size(im_roi,2)/2;
        
        % Plot image
        figure(f)
        
        % TOP IMAGE  ---------------------
        subplot(2,1,1)
        h = imshow(im,'InitialMag','fit');
        title(['Frame ' num2str(S.frames(i))])
        a1 = gca;

        hold on       
        % Outline roi (global)
%         line(S.roi(i).xCntr,S.roi(i).yCntr,'Color',...
%             [clrs{2} 0.5],'LineWidth',2);
        line(S.roi(i).xPerimG,S.roi(i).yPerimG,'Color',...
            [clrs{2} 0.5],'LineWidth',1);    
        scatter(tipG(1,:),tipG(2,:),'MarkerEdgeColor',clrs{3},...
            'LineWidth',1,'MarkerEdgeAlpha',0.5,'SizeData',mSize(2));
        hold off
        
        
        % BOTTOM IMAGE  ---------------------
        subplot(2,1,2)
        h = imshow(im_roi,'InitialMag','fit');
        hold on
        scatter(tipROI(1,:),tipROI(2,:),'MarkerEdgeColor',clrs{3},...
            'LineWidth',2,'MarkerEdgeAlpha',0.5,'SizeData',mSize(1));
        a2 = gca;
        hold off
        
        
        % Set position of both images
        set(a1,'Position',[0 0.3 1 1])
        set(a2,'Position',[0.1 -0.1 0.8 0.8])
        
        clear tipG tipL tipROI
    end
    
    % Get image
    imFrame = getframe(gcf);
    
    % Write frame
    writeVideo(vOut,imFrame);
end

close(f)

% Output
%varargout{1} = M;
%varargout{1} = [];

close(vOut)




function ptsL = G2L(cOrigin,theta,ptsG)
% Coordinate transformation from global to local coordinates
% cOrgin    - Coordinate origin (3x1)
% theta     - Azimuth angle of local system wrt global
% ptsG      - Coordinates in the global FOR (3xn)

if (size(cOrigin,2) > size(cOrigin,1)) || ...
   ((length(ptsG)>3) && (size(ptsG,1) > size(ptsG,2)))  
    error('All points should be arrange in column vectors');
end

% Remove nans
idx = ~isnan(ptsG(1,:));
ptsG = ptsG(:,idx);

% Add z-dimension
cOrigin = [cOrigin; 0];
ptsG    = [ptsG; zeros(1,size(ptsG,2))];

% Axes defined
xaxis = [cos(theta) sin(theta) 0]';
yaxis = [cos(theta+pi/2) sin(theta+pi/2) 0]';
zaxis = [0 0 1]';

% Rotation matrix
S = [xaxis yaxis zaxis];

% Local coordinate
ptsL = global2localcoord(ptsG,'rr',cOrigin,S);

% Remove z-dimension
ptsL = ptsL(1:2,:);