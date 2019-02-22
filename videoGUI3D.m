function H = videoGUI3D(vid_path_s,v_s,vid_path_c,v_c,S,cam,acqMode,varargin)
% Interactive acquisition of coordinate points from a video 

%% Parse inputs


if strcmp(acqMode,'arms')
    
    % Path for data
    if length(varargin)>0
        savePath = varargin{1};
    end
    
    H  = varargin{2};
    
    % Marker color
    mClr = [0 1 0];
    
    imInvert = 0;
    
    % Frame intervals for data collection 
    frInterval = 15;

else
    error('Do not recognize acqMode');
end

% Frame numbers for each camera
frames_c = cam.c.frames;
frames_s = cam.s.frames;


%% Default parameters

% Container of handle data
%H = [];

% Adjust window docking
set(0,'DefaultFigureWindowStyle','normal')

 % Initial frame
im0_c = getFrame(vid_path_c,v_c,frames_c(1),imInvert,'gray');   
im0_s = getFrame(vid_path_s,v_s,frames_s(1),imInvert,'gray');   

% Set starting roi
%roi = calcROI(cam,frames_s(1),frames_c(1));
%roi = [0 0 size(im0_s,2) size(im0_s,1)];


% Display options
disp(' ')
disp('COMMANDS  ------------------------------------------------ ')
disp('    left click  : select position in current frame')
disp('    delete      : delete current tubefoot data')
% disp('    right arrow : advance one frame')
% disp('    left arrow  : go back one frame')
disp('    up arrow    : advance multiple frames')
disp('    down arrow  : go back multiple frames')
% disp('    j           : set number of frames to jump for up/down arrows')
disp('    1 - 9       : jump to relative position in video segment')
disp('    f           : Jump to frame number')
disp('    t           : Jump to time')
disp('    z           : Toggle zoom mode')
disp('    s           : Save current tubefoot coordinates')
disp('    m           : Mode toggle: eyes/feet')
disp('    q           : quit interaction mode')
disp(' ')
 
% Create figure window
[hFig, hAxes] = createFigureAndAxes;

% Add data to figure
[hFig, hAxes] = putData(hFig, hAxes, H);

% Wait for completion of interactive mode
waitfor(hFig)


    function [hFig, hAxes] = createFigureAndAxes  
        
        % Close figure opened by last run
        figTag = 'CVST_VideoOnAxis_9804532';
        close(findobj('tag',figTag));
        
         AR_c = size(im0_c,2)/size(im0_c,1);
         AR_s = size(im0_s,2)/size(im0_s,1);
%         %AR = range(iC.yTank)/range(iC.xTank);
%         AR_roi = roi(4)/roi(3);
        
        % Relative width of figure window
        rel_width = 0.4;
        rel_corr = 0.2;
        
        % Create new figure
        hFig = figure('numbertitle', 'off', ...
            'name', 'Analysis GUI', ...
            'menubar','none', ...
            'toolbar','none', ...
            'resize', 'on', ...
            'tag',figTag, ...
            'renderer','painters',...
            'Units','pixels');
        
        mPos = get(0,'MonitorPositions');
        
        if size(mPos,1)==3
            scrsize = mPos(3,:);
        else       
            scrsize = get(0,'screensize');          
        end
        
        set(hFig,'Position',scrsize,'Color','k');
        
        %set(hFig,'Units','normalized');
        
        % Size of width and height of 2 windows
        win_size_c = [scrsize(3)*(rel_width-rel_corr) scrsize(3)*(rel_width-rel_corr)/AR_c];
        win_size_s = [scrsize(3)*(rel_width+rel_corr) scrsize(3)*(rel_width+rel_corr)/AR_s];
        
        % Spacing between windows
        spacer    = 10;
        
        % Relative position of main panel
        Lrect(3) = win_size_c(1);
        Lrect(4) = win_size_c(2);
        Lrect(1) = scrsize(3)/4 - Lrect(3) - spacer/2;
        Lrect(2) = scrsize(4)/2 - Lrect(4)/2;
        
        Ltitle(3) = Lrect(3);
        Ltitle(4) = 20;
        Ltitle(1) = Lrect(1);
        Ltitle(2) = Lrect(2) + Lrect(4);    
        
        % Relative position of small panel
        Rrect(3) = win_size_s(1);
        Rrect(4) = win_size_s(1)/AR_s;
        Rrect(1) = scrsize(3)/4 + spacer/2;
        Rrect(2) = scrsize(4)/2 - Rrect(4)/2;
        
        Rtitle(3) = Rrect(3);
        Rtitle(4) = 20;
        Rtitle(1) = Rrect(1);
        Rtitle(2) = Rrect(2) + Rrect(4) ;
        
        % Get figure position
        %fPos = get(hFig,'Position');
        
        % Create axes and titles
        hAxes.axis1 = createPanelAxisTitle(hFig,...
            Lrect, Ltitle,'Full frame', 'title1'); % [X Y W H]
        
        hAxes.axis2 = createPanelAxisTitle(hFig, ...
            Rrect, Rtitle, 'ROI', 'title2');

        % Keystrokes
        set(hFig, 'WindowButtonDownFcn', {@butDown, hFig, hAxes});
        
        % Set callback for button press
        set(hFig, 'WindowKeyPressFcn', {@keyPress, hFig, hAxes});

    end


%% Deposit data into figure

    function [hFig, hAxes] = putData(hFig, hAxes,H)

        
        if isempty(H)
            % Store inital values into axis1
            %H.roi = roi;
            H.iFrame   = 1;
            H.frames_c = min(frames_c):frInterval:max(frames_c);
            H.frames_s = min(frames_s):frInterval:max(frames_s);
            H.vid_path_c = vid_path_c;
            H.vid_path_s = vid_path_s;
            H.savePath = savePath;
            
            H.clr = mClr;
            %         H.x = xStart;
            %         H.y = yStart;
            H.imInvert = imInvert;
            H.adInterval = 20;
            H.currArm = [];
            H.S = S;
            H.acqMode = 'eyes';
            H.cam = cam;
            H.zoomOn = 1;
            H.zoomFactor = 1.7;
            H.yCntr = [];
            H.frInterval = frInterval;
        end
        
        H.v_c = v_c;
        H.v_s = v_s;
  
        % If there is no eye field, start at 1
        if ~isfield(H.S,'eye')
            
            H = defineStartingData(H);
           
        end
        
        clear numArm
                
        % Store coordinate data
        guidata(hAxes.axis1, H);
        
        % Render images in figure
        update_fig(hFig, hAxes);
        
        %acqModes: 'Choose arm', 'Track eye', 'Track feet
    end


%% Create Axis and Title
% Axis is created on uipanel container object. This allows more control
% over the layout of the GUI. Video title is created using uicontrol.
    function hAxis = createPanelAxisTitle(hFig, posWin, posTitle, axisTitle, textTag)

        % Create panel
        hPanel = uipanel('parent',hFig,'Units','pixels','Visible','on');
        
        % Set position
        hPanel.Position = posWin;
        hPanel.Units = 'pixels';
        
        % Create axis   
        hAxis = axes('Parent',hPanel,'Units','pixels'); 
        hAxis.Position = [0 0 posWin(3:4)];
        
        hAxis.XTick = [];
        hAxis.YTick = [];
        hAxis.XColor = [1 1 1];
        hAxis.YColor = [1 1 1];
        
        % Revert units to pixels
        hAxis.Units = 'pixels';
        
        % Set video title using uicontrol. uicontrol is used so that text
        % can be positioned in the context of the figure, not the axis.
        titlePos = posTitle;
        hUI = uicontrol('style','text',...
                        'String', axisTitle,...
                        'Units','pixels',...
                        'Parent',hFig,...
                        'Position', titlePos,...
                        'BackgroundColor',[1 1 1], ...
                        'HorizontalAlignment','left',...
                        'Tag',textTag);
            
            %'BackgroundColor',hFig.Color, ...
        
        % Revert units to pixels
        %hUI.Units = 'pixels';
    end
end



%% Small helper functions

function showFrameOnAxis(hAxis, frame, zoomlevel)
    % This helper function  displays a frame of video on a user-defined axis.

    frame = convertToUint8RGB(frame);

    try
        hChild = get(hAxis, 'Children');
    catch %#ok<CTCH>
        return; % hAxis does not exist; nothing to draw
    end

    isFirstTime = isempty(hChild);

    if isFirstTime
        hIm = displayImage(hAxis, frame);
        zoom(hAxis,zoomlevel)
        if 0
            addScrollPanel(hAxis, hIm);
        end
    else
        hIm = hChild(end);

        try
            set(hIm,'cdata',frame); 
            drawnow;
        catch  %#ok<CTCH>
            % figure closed
            return;
        end
    end
end

function frame = convertToUint8RGB(frame)
    % Convert input data type to uint8
    if ~isa(class(frame), 'uint8')
        frame = im2uint8(frame);
    end

    % If the input is grayscale, turn it into an RGB image
    if (size(frame,3) ~= 3) % must be 2d
        frame = cat(3,frame, frame, frame);
    end
end

function hIm = displayImage(hAxis, frame)
% Display image in the specified axis
frameSize = size(frame);
xdata = [1 frameSize(2)];
ydata = [1 frameSize(1)];

%xdata= [1 hAxis.Parent.Position(3)];
%ydata= [1 hAxis.Parent.Position(4)];
%xdata = [0 range(hAxis.XLim)];
%ydata = [0 range(hAxis.YLim)];
cdatamapping = 'direct';
%cdatamapping = 'scaled';


hIm = image(xdata,ydata,frame, ...
           'BusyAction', 'cancel', ...
           'Parent', hAxis, ...
           'CDataMapping', cdatamapping, ...
           'Interruptible', 'off');
set(hAxis, ...
    'YDir','reverse',...
    'TickDir', 'out', ...
    'XGrid', 'off', ...
    'YGrid', 'off', ...
    'PlotBoxAspectRatioMode', 'auto', ...
    'Visible', 'off');

ttt=3;
end




%% Key press callback

function keyPress(fig, key, hFig, hAxes)
            
    closefig = 0;

     % Load data       
     H = guidata(hAxes.axis2);
         
     
     % QUIT ('q')
     if strcmp(key.Key,'Q') || strcmp(key.Key,'q')
         
          % Prompt   
          bName = questdlg(['Do you want to save you work before you quit?' ...
              'Save data?'],'','Yes','No, Cancel','Yes');
          
          % If yes . . .
          if strcmp(bName,'Yes')
              saveData(H)
          end
          
          assignin('caller','H',H);
          closefig = 1;

     % TOGGLE MODE ('m')
     elseif strcmp(key.Key,'M') || strcmp(key.Key,'m')

        if strcmp(H.acqMode,'eyes')
            
            % Change to feet mode
            H.acqMode = 'feet';
            
            % If in foot acquisition mode, point is a nan and neightboring frame 
             % is not a nan . . .
             if  max(isnan(H.S.ft_s(H.iFrame).x)) && ...
                     (H.iFrame > 1) && max(~isnan(H.S.ft_s(H.iFrame-1).x))
                 % Make index one less
                 H.iFrame = H.iFrame-1;
             end

        else
            % Change to eye mode
            H.acqMode = 'eyes';
       
        end
        
        % Store data
        guidata(hAxes.axis1, H);
    
     % DELETE
     elseif strcmp(key.Key,'backspace')
         
         if strcmp(H.acqMode,'eyes')
             % Add nans
             H.S.eye(H.currArm).s.x(H.iFrame) = nan;
             H.S.eye(H.currArm).s.z(H.iFrame) = nan;
             
         elseif strcmp(H.acqMode,'feet')
             % Add nans
             H.S.ft_s(H.iFrame).x(end) = nan;
             H.S.ft_s(H.iFrame).z(end) = nan;
         end

       % Store coordinate data
       guidata(hAxes.axis1, H);
       
       % Update figure
       update_fig(hFig, hAxes)
     
     % (zoom on/off)
     elseif strcmp(key.Key,'z') || strcmp(key.Key,'Z')

         % Toggle zoom
         H.zoomOn = abs(H.zoomOn-1);
         
         % Store coordinate data
          guidata(hAxes.axis1, H);
    
          % Update figure
          update_fig(hFig, hAxes)

%       % J (set interval to jump over frames)
%       elseif strcmp(key.Key,'j') || strcmp(key.Key,'J')
%           
%           % Prompt for interval
%           answer = inputdlg({'Frame interval'},'Set frame interval',1,...
%                             {num2str(H.adInterval)});                   
%           if isempty(answer)
%               return
%           end
%           
%           % Advance interval
%           H.adInterval = str2num(answer{1});
%           
%           % Store coordinate data
%           guidata(hAxes.axis1, H);
%     
%           % Update figure
%           update_fig(hFig, hAxes)                  
          
    % S (save data)
      elseif strcmp(key.Key,'s') || strcmp(key.Key,'S')  
 
          saveData(H)
             
          
     % The following involve changing frame number ---------------------
     else
         newIndex = [];
         
         % Index of current frame
%         iFrame = find(H.cFrame==H.frames,1,'first');
         
         % RIGHT ARROW/SPACE
         if strcmp(key.Key,'rightarrow') || strcmp(key.Key,'space')
             
             if strcmp(H.acqMode,'eyes')
                 
                 if (H.iFrame+1)>length(H.frames_s)
                     beep
                 else
                     % Advance frame
                     newIndex = H.iFrame+1;
                 end
                 
             elseif strcmp(H.acqMode,'feet')
                 if (H.iFrame+2)>length(H.frames_s)
                     beep
                 else
                     % Advance frame
                     newIndex = H.iFrame+2;
                 end
             end
             
         % LEFT ARROW
         elseif strcmp(key.Key,'leftarrow')
             
             if strcmp(H.acqMode,'eyes')
                 
                 if (H.iFrame-1)<1
                     beep
                 else
                     % Reduce frame
                     newIndex = H.iFrame-1;
                 end
                 
             elseif strcmp(H.acqMode,'feet')
                 if (H.iFrame-2)<1
                     beep
                 else
                     % Reduce frame
                     newIndex = H.iFrame-2;
                 end
             end
             
%          % UP ARROW
%          if strcmp(key.Key,'uparrow') 
%              
%              % Advance frame
%              newIndex = min([H.iFrame+H.adInterval length(H.frames_c) ...
%                              length(H.frames_s)]);
%              
%          % DOWN ARROW
%          elseif strcmp(key.Key,'downarrow')
%              
%              % Reverse frame
%              newIndex = max([1 H.iFrame-H.adInterval]);
             
         % NUMBER
         elseif length(key.Key)==1 && (sum(key.Key==num2str([1:9]))==1)
             
             % Requested relative number
             req_num = str2num(key.Key);
             
             % Start 
             if req_num == 1
                 newIndex = 1;
                 
             % End
             elseif req_num==9
                 newIndex = min([length(H.frames_s) length(H.frames_c)]);
                 
             % Something between
             else
                 % Set new frame
                 newIndex = round((req_num/10) * ...
                                 min([length(H.frames_s) length(H.frames_c)]));
             end
               
         % Key in time (in s)
         elseif strcmp(key.Key,'t')|| strcmp(key.Key,'T')
             
             answer = inputdlg({'Time: min','Time: sec'},'',1,{'0',''});
             
             % If answer provided
             if ~isempty(answer)
                 
                % Convert time into seconds
                timeval = str2num(answer{1})*60 + str2num(answer{2});
                 
                % Translate into frames
                aFrame = ceil(timeval*H.v.FrameRate);
                
                if aFrame<min(H.frames_s)
                    warning('Requested time happens before analyzed interval');
                    newIndex = [];
                    
                elseif aFrame>max(H.frames_s)
                    warning('Requested time happens after analyzed interval');
                    newIndex = [];
                    
                else
                    
                    % Difference between desired and available frames
                    frameDiff = aFrame - H.frames_s;
                
                    % New index as closest to requested
                    newIndex = find(frameDiff==min(frameDiff),1,'first');
                    
                end
             
             % If no answer
             else
                 newIndex = [];
             end
             
         % Key in frame number
         elseif strcmp(key.Key,'f')|| strcmp(key.Key,'F')
             
             answer = inputdlg({'Time: frame'},'',1,{''});
             
             % If answer provided
             if ~isempty(answer)
                 
                % Translate into frames
                aFrame = round(str2num(answer{1}));
                
                if aFrame<min(H.frames_s)
                    warning('Requested time happens before analyzed interval');
                    newIndex = [];
                    
                elseif aFrame>max(H.frames_s)
                    warning('Requested time happens after analyzed interval');
                    newIndex = [];
                    
                else
                    newIndex = find(aFrame==H.frames_s,1,'first');
                end
             
             % If no answer
             else
                 newIndex = [];
             end
             
         end

         % If new newframe
         if ~isempty(newIndex)

             % If in foot acquisition mode, piint is a nan and neightboring frame 
             % is not a nan . . .
             if strcmp(H.acqMode,'feet') && max(isnan(H.S.ft_s(newIndex).x)) && ...
                     (newIndex > 1) && max(~isnan(H.S.ft_s(newIndex-1).x))
                 % Make index one less
                 newIndex = newIndex-1;
             end
             
             % Update frame number
             H.iFrame = newIndex;
                         
             % Store coordinate data
             guidata(hAxes.axis1, H);
         end
         
      end
     
    if closefig==1
         close(hFig)
    else
      % Update figure
       update_fig(hFig, hAxes)
    end
end


%% Button down callback
function butDown(fig, key, hFig, hAxes)
              
    % Load data       
    H = guidata(hAxes.axis2);   
     
    % Collect current coordinate
    C2 = get(hAxes.axis2, 'CurrentPoint');
    C1 = get(hAxes.axis1, 'CurrentPoint');
    
    %H = updateArm(H)
    hold on
    
    % If in left frame . . .
    if (C1(1,1)>=0) && (C1(1,1)<=hAxes.axis1.XLim(2)) && ...
            (C1(1,2)>=0) && (C1(1,2)<=hAxes.axis1.YLim(2))
        
        if strcmp(H.acqMode,'eyes')
            % Current coordinate
            xCurr = C1(1,1);
            yCurr = C1(1,2);
            
            % Coordinates for arms in current frame
            for i = 1:length(H.S.eye)
                xArms(i,1) = H.S.eye(i).c.x(H.iFrame);
                yArms(i,1) = H.S.eye(i).c.y(H.iFrame);
            end
            
            % Distance to arms
            dists = hypot(xArms-xCurr,yArms-yCurr);
            
            % Index of closest eye
            iArm = find(dists==min(dists),1,'first');
            
            
            % Insert current coordinate as closest
            H.S.eye(iArm).c.x(H.iFrame) = xArms(iArm);
            H.S.eye(iArm).c.y(H.iFrame) = yArms(iArm);
            
            % Store curret arm
            H.currArm = iArm;
            
            % Clear variables
            clear iArm iMin dists xArms yArms i
            
        % If in foot mode
        else
            warning('Foot mode only works in the Sony view');
        end
    
    % If in right frame . . .
    elseif (C2(1,1)>=0) && (C2(1,1)<=hAxes.axis2.XLim(2)) && ...
           (C2(1,2)>=0) && (C2(1,2)<=hAxes.axis2.YLim(2))     
    
       % Look for arm selection
       if strcmp(H.acqMode,'eyes') && isempty(H.currArm)
           warning('You first need to select an arm in the left frame');
           
       % If selecting eye coords
       elseif strcmp(H.acqMode,'eyes')
           
           % Store
           H.S.eye(H.currArm).s.x(H.iFrame) = C2(1,1);
           H.S.eye(H.currArm).s.z(H.iFrame) = C2(1,2);
           
           % y-Position for center of zoomed view
           H.yCntr = C2(1,2);
           
       % If selecting feet
       elseif strcmp(H.acqMode,'feet')
           
           % Index for new value
           if isnan(H.S.ft_s(H.iFrame).x(end))
               iAdd = 1;
           else
               iAdd = length(H.S.ft_s(H.iFrame).x) + 1;
           end
           
           % Add
           H.S.ft_s(H.iFrame).x(iAdd) = C2(1,1);
           H.S.ft_s(H.iFrame).z(iAdd) = C2(1,2);
           
           % y-Position for center of zoomed view
           H.yCntr = C2(1,2);
           
       end
       
    % If outside both camera views
    else
        set(gcf,'Pointer','arrow')
    end

    % Store coordinate data
    guidata(hAxes.axis1, H);
    
    % Update figure
    update_fig(hFig, hAxes);           
end


%% Update GUI
function update_fig(hFig, hAxes)

% Activate figure window
figure(hFig)

% Load data
H = guidata(hAxes.axis1);

% List mode in title 2
obj = findobj('tag','title2');
set(obj,'String','');


% CANON CAMERA --------------------------------------------------

% Define time variables
currTime_c = (H.frames_c(H.iFrame))./H.v_c.FrameRate;
currMin_c = floor(currTime_c/60);
currSec_c = floor(currTime_c - currMin_c*60);
currFrac_c = floor((currTime_c - currMin_c*60-currSec_c)*100);

% Convert to strings
currMin_c  = num2str(currMin_c);
currSec_c  = ['0' num2str(currSec_c)];
currFrac_c = num2str(currFrac_c);

% List frame number in title 1
obj = findobj('tag','title1');

if strcmp(H.acqMode,'eyes')
    % Set title text
    t_str = ['Canon :  Frame ' num2str(H.frames_c(H.iFrame)) ' (' ...
        currMin_c ':' currSec_c(end-1:end) ':' currFrac_c '). SELECT EYES.'];
elseif strcmp(H.acqMode,'feet')
    % Set title text
    t_str = ['Canon :  Frame ' num2str(H.frames_c(H.iFrame)) ' (' ...
        currMin_c ':' currSec_c(end-1:end) ':' currFrac_c ').'];
end


set(obj,'String',t_str);

% Current frame
cFrame_c = H.frames_c(H.iFrame);

% Read input video frame
frame_c = getFrame(H.vid_path_c,H.v_c,cFrame_c,H.imInvert,'rgb');

% Display full video frame
delete(hAxes.axis1.Children)
showFrameOnAxis(hAxes.axis1, frame_c, 0);

% Marker size
mSize = 100;

% If selecting eye coords
if strcmp(H.acqMode,'eyes')
    
    % Hold on Axis1
    set(hAxes.axis1,'NextPlot','Add')
    
    % Arm coordinates
    xArm_c = H.S.arm(H.iFrame).x;
    yArm_c = H.S.arm(H.iFrame).y;
    
    % Center coordinate
    xCntr_c = H.S.xCntr(H.iFrame);
    yCntr_c = H.S.yCntr(H.iFrame);
    
    % Color map
    %cmap = colormap('lines');
    
    % Arm orientations
    aTheta = atan2(yArm_c-yCntr_c,xArm_c-xCntr_c);
    
    % Arms directed toward camera
    idx = (aTheta > 0) | (aTheta < -160/180*pi) | (aTheta > -20/180*pi);
    %idx = (aTheta < 0) | (aTheta > 160/180*pi) | (aTheta < 20/180*pi);
    
    xArm_c = xArm_c(idx);
    yArm_c = yArm_c(idx);
    
    %cvals = cmap(1:length(xArm_c),:);
    
    % White circles on all arms
    h2 = scatter(xArm_c,yArm_c,'Parent',hAxes.axis1,...
        'MarkerfaceColor','none','SizeData',mSize,'MarkerEdgeColor','w');
    
    % If there is a current arm
    if  ~isempty(H.currArm)
        
        % Red circle on current arm
        h2 = scatter(H.S.eye(H.currArm).c.x(H.iFrame),...
            H.S.eye(H.currArm).c.y(H.iFrame),...
            'Parent',hAxes.axis1,...
            'MarkerFaceColor','none','SizeData',mSize,...
            'MarkerEdgeColor','r');
    end
    
    % Dot point
    plot(H.cam.c.xDot,H.cam.c.yDot,'+r')
    
    % Hold off Axis1
    set(hAxes.axis1,'NextPlot','Replace');
end


% SONY CAMERA -------------------------------------------------------

% Define time variables
currTime_s = (H.frames_s(H.iFrame))./H.v_s.FrameRate;
currMin_s = floor(currTime_s/60);
currSec_s = floor(currTime_s - currMin_s*60);
currFrac_s = floor((currTime_s - currMin_s*60-currSec_s)*100);

% Convert to strings
currMin_s  = num2str(currMin_s);
currSec_s  = ['0' num2str(currSec_s)];
currFrac_s = num2str(currFrac_s);

% List frame number in title 2
obj = findobj('tag','title2');

% Set title text
if strcmp(H.acqMode,'eyes')
    t_str = ['Sony :  Frame ' num2str(H.frames_s(H.iFrame)) ' (' ...
        currMin_s ':' currSec_s(end-1:end) ':' currFrac_s '). SELECT EYES.'];
elseif strcmp(H.acqMode,'feet')
    t_str = ['Sony :  Frame ' num2str(H.frames_s(H.iFrame)) ' (' ...
        currMin_s ':' currSec_s(end-1:end) ':' currFrac_s '). SELECT BASE OF FEET.'];
end

set(obj,'String',t_str);

% Current frame
cFrame_s = H.frames_s(H.iFrame);

% Read input video frame
frame_s = getFrame(H.vid_path_s,H.v_s,cFrame_s,H.imInvert,'rgb');

% Display image
delete(hAxes.axis2.Children);
showFrameOnAxis(hAxes.axis2, frame_s, 0);

% Adjust zoom
if H.zoomOn && ~isempty(H.yCntr)
    adjustZoom(H.cam,H.yCntr,cFrame_s,cFrame_c,hAxes.axis2,H.zoomFactor);
end

% Limits for x and y axes
xL = hAxes.axis2.XLim;
yL = hAxes.axis2.YLim;

% Hold on Axis 2
set(hAxes.axis2,'NextPlot','Add');

% Activate right panel
set(hAxes.axis2,'Selected','on')

% If selecting eye coords
if strcmp(H.acqMode,'eyes')
    
    % Common dot point
    plot(H.cam.s.xDot,H.cam.s.yDot,'+r')
    
    % Find estimated x-position of unselected points
    xArm_s = canon2sonyX(xArm_c,H.cam,cFrame_s,cFrame_c);
    
    % Plot unselected points
    h2 = line(xArm_s.*[1 1],yL,'LineWidth',3,'Color',[1 1 1 0.1]);
    
    % If there is a selected arm . . .
    if ~isempty(H.currArm)
        
        % Find sony cooridnates from canon coordinates
        xArm_s = canon2sonyX(H.S.eye(H.currArm).c.x(H.iFrame),H.cam,cFrame_s,cFrame_c);
        
        % Plot line for selected arm
        h3 = line(xArm_s.*[1 1],yL,'LineWidth',2,'Color',[1 0 0 0.15]);
        
        % Index for values
        idx = ~isnan(H.S.eye(H.currArm).s.z);
        
        if sum(idx)>1
            % Trajectory coordinates
            xS = H.S.eye(H.currArm).s.x(idx);
            yS = H.S.eye(H.currArm).s.z(idx);
            
            % Draw trajectory
            h2 = line(xS,yS,'LineWidth',6,'Color',[0.25 0.25 0 0.3]);
        end
        
        
        % If there is a y-point selected . . .
        if ~isnan(H.S.eye(H.currArm).s.z(H.iFrame))
            
            % Draw selected point
            h2 = scatter(H.S.eye(H.currArm).s.x(H.iFrame),...
                H.S.eye(H.currArm).s.z(H.iFrame),...
                'MarkerFaceColor','none',...
                'MarkerEdgeColor','r',...
                'SizeData',mSize*6*H.zoomFactor,'MarkerEdgeAlpha',0.2,...
                'LineWidth',4);
        end
        
        
    end
    
% If selecting feet . . .
else
    
    % If any points on this frame . . .
    if sum(~isnan(H.S.ft_s(H.iFrame).x))>0
        
        lineWidth = 75;
        
        % Loop thru points
        for i = 1:length(H.S.ft_s(H.iFrame).x )
            line([(H.S.ft_s(H.iFrame).x(i) - lineWidth/2) ...
                (H.S.ft_s(H.iFrame).x(i) + lineWidth/2)],...
                H.S.ft_s(H.iFrame).z(i).*[1 1], ...
                'LineWidth',3,'Color',[1 0 0 0.2])
        end
        
        clear i lineWidth
    end
end

% Hold off Axis 2
set(hAxes.axis2,'NextPlot','Replace');

% Store coordinate data
%guidata(hAxes.axis1, H);

% Activate figure
figure(hFig)


end

    

function xS = canon2sonyX(xC,cam,cFrame_s,cFrame_c)
% Determine x coordinates (in pix) for sony camera for canon points
% 
% % Scale factor (correct for parallax)
% scFactor = polyval(cam.relsize.coef,cFrame);
% 
% % Define dot point for canon (m)
% xDot_c = cam.c.xDot * cam.c.cal;

% Margins for canon
L_c = polyval(cam.c.mLcoef,cFrame_c);
R_c = polyval(cam.c.mRcoef,cFrame_c);

% Margins for sony
L_s = polyval(cam.s.mLcoef,cFrame_s);
R_s = polyval(cam.s.mRcoef,cFrame_s);

% Body width in pixels
bWidth_c = abs(R_c-L_c);
bWidth_s = abs(R_s-L_s);

% Body center in pixels
bCntr_c = mean([R_c L_c]);
bCntr_s = mean([R_s L_s]);

% Points relative to center (pix)
xC = xC - bCntr_c;

% Points normalized by body width in canon (BLs)
xC = xC ./ bWidth_c;

% Points in pix for sony (pix)
xS = xC .* bWidth_s; 

% Points translated into full frame roi (pix)
xS = bCntr_s - xS;


end


function adjustZoom(cam,yVal_s,cFrame_s,cFrame_c,axis2,scFctr)
% Determine zoom for sony view

% Width/height
AR = range(axis2.XLim)/range(axis2.YLim);

% Scaling factor for the width
%scFctr = 1.7;

% Margins for canon
L_c = polyval(cam.c.mLcoef,cFrame_c);
R_c = polyval(cam.c.mRcoef,cFrame_c);

% Margins for sony
L_s = polyval(cam.s.mLcoef,cFrame_s);
R_s = polyval(cam.s.mRcoef,cFrame_s);

% Body width in pixels
bWidth_c = abs(R_c-L_c);
bWidth_s = abs(R_s-L_s);

% Body center in pixels
bCntr_c = mean([R_c L_c]);
bCntr_s = mean([R_s L_s]);

% Width of roi
w = bWidth_s * scFctr;

% Height
h = w/AR;

% Adjust axes
axis2.XLim = [bCntr_s-w/2 bCntr_s+w/2];
axis2.YLim = [yVal_s-h/2  yVal_s+h/2];

end


function H = defineStartingData(H)

% Number of arms
numArm = length(H.S.arm(1).x);

% Number of frames
numFrame = length(H.frames_c);


% Loop thru eyes
for i = 1:numArm
    % Make slots for eye coordinates
    H.S.eye(i).c.x = nan(numFrame,1);
    H.S.eye(i).c.y = nan(numFrame,1);
    H.S.eye(i).s.x = nan(numFrame,1);
    H.S.eye(i).s.z = nan(numFrame,1);
    
end

% Loop thru frames for canon view
for i = 1:numFrame
    % Break, if points for all arms defined
    if sum(~isnan(H.S.arm(i).x))==numArm
        iStart = i;
        break
    end
end

% Make slots for feet coordinates
for i = 1:numFrame
    H.S.ft_s(i).x = nan;
    H.S.ft_s(i).z = nan;
end

% Coordinates of starting arms in canon
xArmStart = H.S.arm(iStart).x;
yArmStart = H.S.arm(iStart).y;

% Loop thru arms
for i = 1:numArm
    
    % Get starting position of current arm
    H.S.eye(i).c.x(iStart) = xArmStart(i);
    H.S.eye(i).c.y(iStart) = yArmStart(i);
    
    % Loop thru frames
    for j = (iStart+1):length(H.S.arm)
        
        % Coordinates for arms in current frame
        xArms = H.S.arm(j).x;
        yArms = H.S.arm(j).y;
        
        % Coord for prior
        xPrior = H.S.eye(i).c.x(j-1);
        yPrior = H.S.eye(i).c.y(j-1);
        
        % Distance of current arms to prior
        dists = hypot(xPrior-xArms,yPrior-yArms);
        
        % Index of closest eye/arm
        iMin = find(dists==min(dists),1,'first');
        
        % Insert current coordinate
        H.S.eye(i).c.x(j) = xArms(iMin);
        H.S.eye(i).c.y(j) = yArms(iMin);
        
        clear xArms yArms xPrior yPrior iMin
    end
end

end

function saveData(H)
% Strip out v structures and save

% Save data (without 'v' field)
v_s = H.v_s;
v_c = H.v_c;

% Remove fields
H = rmfield(H,'v_c');
H = rmfield(H,'v_s');

% Update status
disp(['Saving to ' H.savePath])

% Save data
save(H.savePath,'H')

end


