function A = videoGUI(vid_path,v,frames,imInvert,acqMode,varargin)
% Interactive acquisition of coordinate points from a video 
% 

%% Parse inputs

if nargin< 5
    acqMode = 'simple';
end

if strcmp(acqMode,'simple')
    
    % roi radius
    if length(varargin)>0
        r = varargin{1};
    else
        r = 0;
    end
    
    % Marker color
    if length(varargin)>1
        mClr = varargin{2};
    else
        mClr = [0 1 0];
    end
    
    A = [];
    
   A.savePath = varargin{3};
   A.numArms  = varargin{4};

else
    error('Do not recognize acqMode');
end


%% Default parameters

% Container of handle data
%A = [];

% Adjust window docking
set(0,'DefaultFigureWindowStyle','normal')

 % Initial frame
im0 = getFrame(vid_path,v,frames(1),imInvert,'gray');    

% Alter default
if r == 0
    r = min(size(im0))/5;
end

% Apply mask
%im0 = applyMask(im0,iC.xTank,iC.yTank);

% Set starting roi
% roi(1,1) = size(im0,2)/2-r(1)/2;
% roi(1,2) = size(im0,1)/2-r(1)/2;
roi(1,3) = min([3*r size(im0,1) size(im0,2)]);
roi(1,4) = roi(1,3);
roi(1,1) = size(im0,2)/2-roi(1,3)/2;
roi(1,2) = size(im0,1)/2-roi(1,3)/2;


% Display options
disp(' ')
disp('COMMANDS  ------------------------------------------------ ')
disp('    left click  : select position in current frame')
disp('    delete      : delete current coordinate')
disp('    right arrow : advance one frame')
disp('    left arrow  : go back one frame')
%disp('    1 - 9       : jump to relative position in video segment')
disp('    +           : zoom in')
disp('    -           : zoom out')
%disp('    i           : invert image')
%disp('    d           : delete point in current frame')
%disp('    a           : delete all points in interval')
disp('    n           : Choose next arm')
disp('    q           : quit interactive mode and save data')
disp(' ')
 
% Create figure window
[hFig, hAxes] = createFigureAndAxes;

% Add data to figure
[hFig, hAxes] = putData(hFig, hAxes);

% Wait for completion of interactive mode
waitfor(hFig)

% If current arm is incomplete . . .
if sum(~isnan(A.arm(end).x))==0 

    % And still on the first arm . . .
    if length(A.arm)==1
        return
       
    % If beyond first . . .
    else        
        % Delete data for current tube foot
        A.arm(end).x         = nan(length(A.arm(end).x),1);
        A.arm(end).y         = nan(length(A.arm(end).y),1);
    end
end

    function [hFig, hAxes] = createFigureAndAxes
        
        
        % Close figure opened by last run
        figTag = 'CVST_VideoOnAxis_9804532';
        close(findobj('tag',figTag));
        
        AR = size(im0,2)/size(im0,1);
        %AR = range(iC.yTank)/range(iC.xTank);
        AR_roi = roi(4)/roi(3);
        
        % Relative width of figure window
        rel_width = 0.4;
        
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
        win_size = [scrsize(3)*rel_width scrsize(3)*rel_width/AR];
        
        % Spacing between windows
        spacer    = 10;
        
        % Relative position of main panel
        Lrect(3) = win_size(1);
        Lrect(4) = win_size(2);
        Lrect(1) = scrsize(3)/2 - Lrect(3) - spacer/2;
        Lrect(2) = scrsize(4)/2 - Lrect(4)/2;
        
        Ltitle(3) = Lrect(3);
        Ltitle(4) = 20;
        Ltitle(1) = Lrect(1);
        Ltitle(2) = Lrect(2) + Lrect(4);    
        
        % Relative position of small panel
        Rrect(3) = win_size(1);
        Rrect(4) = win_size(1)/AR_roi;
        Rrect(1) = scrsize(3)/2 + spacer/2;
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

    function [hFig, hAxes] = putData(hFig, hAxes)
        
        % Initialize data repository
        A = makeNewArm(A,length(frames));
        
        % Store inital values into axis1
        A.roi = roi;
        A.iFrame = 1;
        A.frames = frames;
        A.vid_path = vid_path;
        A.v = v;
        A.clr = mClr;
%         A.x = xStart;
%         A.y = yStart;
        A.imInvert = imInvert;
        A.adInterval = 15;
        A.iMode = 'b';

        % Store coordinate data
        guidata(hAxes.axis1, A);
        
        % Render images in figure
        update_fig(hFig, hAxes);
    end


%% Create Axis and Title
% Axis is created on uipanel container object. This allows more control
% over the layout of the GUI. Video title is created using uicontrol.
    function hAxis = createPanelAxisTitle(hFig, posWin, posTitle, axisTitle, textTag)

        % Create panel
        %hPanel = uipanel('parent',hFig,'Position',pos,'Units','Normalized');
        hPanel = uipanel('parent',hFig,'Units','pixels','Visible','on');
        
        % Set position
        hPanel.Position = posWin;
        hPanel.Units = 'pixels';
        
        % Create axis   
        hAxis = axes('Parent',hPanel,'Units','pixels'); 
        hAxis.Position = [0 0 posWin(3:4)];
        %hAxis.Position = [-min(iC.xTank)/2 min(iC.yTank)/2 range(iC.xTank) range(iC.yTank)];
        
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
cdata = frame;
cdatamapping = 'direct';

hIm = image(xdata,ydata,cdata, ...
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
end

% function [x, y] = full_to_roi(xVal,yVal,roi)
% % Converts from full frame coordinates to ROI coords    
% 
%      % Normalize to size of rendered frame
%      xNorm = (xVal-rect(1))/roi(3);
%      yNorm = (yVal-rect(2))/roi(4);
% 
%      % Transform into frame coords
%      x = xNorm .* range(hAxes.axis2.XLim) + hAxes.axis2.XLim(1);
%      y = yNorm .* range(hAxes.axis2.YLim) + hAxes.axis2.YLim(1);
% end


%% Key press callback

function keyPress(fig, key, hFig, hAxes)
            
    closefig = 0;

     % Load data       
     A = guidata(hAxes.axis2);
     
     % QUIT ('q')
      if strcmp(key.Key,'Q') || strcmp(key.Key,'q')
         
         % If there's any unsaved data . . .
         if sum(isnan(A.arm(end).x))>0 
             
             bName = questdlg(['Are you sure you want to quit without ' ...
                 'saving the current tube foot data?'],'','Yes','No, Cancel','Yes');
             
             if strcmp(bName,'Yes')
                assignin('caller','A',A);
                closefig = 1;
             end
             
         % If no unsaveed data  
         else
             
             % If still on first arm, remove arm field   
             if length(A.arm)==1               
                 % Remove arm field
                 A = rmfield(A,'arm');
                 
             else
                 
                  % Save data (without 'v' field)
                  A = rmfield(A,'v');
                  save(A.savePath,'A')
                 
             end
             
             
             % Close figure
             assignin('caller','A',A);
             closefig = 1;
          
         end
         
         
     % DELETE
     elseif strcmp(key.Key,'backspace')
         
       % Delete data
       A.arm(A.iFrame).x         = nan(length(A.arm(end).x),1);
       A.arm(A.iFrame).y         = nan(length(A.arm(end).y),1);
       
       % Store coordinate data
       guidata(hAxes.axis1, A);
       
       % Update figure
       update_fig(hFig, hAxes)
     
     % PLUS (zoom in)
     elseif strcmp(key.Key,'hyphen') || strcmp(key.Key,'minus')
         
         zoomFactor = 1.5;
         
         xCntr = A.roi(1) + A.roi(3)/2;
         yCntr = A.roi(2) + A.roi(4)/2;
         
         % Enlarge roi
         if (zoomFactor*A.roi(3)>A.v.Height) || (zoomFactor*A.roi(4)>A.v.Width)
             warning('You cannot have an ROI larger than the video frame')
         else
             A.roi(3) = zoomFactor*A.roi(3);
             A.roi(4) = zoomFactor*A.roi(4);
         end
         %A.roi(1) = xCntr - A.roi(3)/2;
         %A.roi(2) = yCntr - A.roi(4)/2;
         
         % Update roi x-coordinate
         if (xCntr-A.roi(3)/2) < 0
             A.roi(1) = 1;
         elseif (xCntr+A.roi(3)/2) > A.v.Width
             A.roi(1) = A.v.Width-A.roi(3);
         else
             A.roi(1) = xCntr-A.roi(3)/2;
         end
         
         % Update roi y-coordinate
         if (yCntr-A.roi(4)/2) < 0
             A.roi(2) = 1;
         elseif (yCntr+A.roi(4)/2) > A.v.Height
             A.roi(2) = floor(A.v.Height - A.roi(4));
         else
             A.roi(2) = yCntr-A.roi(4)/2;
         end
         
         % Store coordinate data
         guidata(hAxes.axis1, A);
         
         % Update figure
         update_fig(hFig, hAxes)
       
     % PLUS (zoom out)
      elseif strcmp(key.Key,'equal') || strcmp(key.Key,'plus')
          
          
          xCntr = A.roi(1) + A.roi(3)/2;
          yCntr = A.roi(2) + A.roi(4)/2;
          
          % Enlarge roi
          A.roi(3) = A.roi(3)/1.5;
          A.roi(4) = A.roi(4)/1.5;
          %A.roi(1) = xCntr - A.roi(3)/2;
          %A.roi(2) = yCntr - A.roi(4)/2;
          
          % Update roi x-coordinate
          if (xCntr-A.roi(3)/2) < 0
              A.roi(1) = 1;
          elseif (xCntr+A.roi(3)/2) > A.v.Width
              A.roi(1) = A.v.Width-A.roi(3);
          else
              A.roi(1) = xCntr-A.roi(3)/2;
          end
          
          % Update roi y-coordinate
          if (yCntr-A.roi(4)/2) < 0
              A.roi(2) = 1;
          elseif (yCntr+A.roi(4)/2) > A.v.Height
              A.roi(2) = floor(A.v.Height - A.roi(4));
          else
              A.roi(2) = yCntr-A.roi(4)/2;
          end
          
          
          % Store coordinate data
          guidata(hAxes.axis1, A);
          
          % Update figure
          update_fig(hFig, hAxes)
       

%       % J (set interval to jump over frames)
%       elseif strcmp(key.Key,'j') || strcmp(key.Key,'J')
%           
%           % Prompt for interval
%           answer = inputdlg({'Frame interval'},'Set frame interval',1,...
%                             {num2str(A.adInterval)});                   
%           if isempty(answer)
%               return
%           end
%           
%           % Advance interval
%           A.adInterval = str2num(answer{1});
%           
%           % Store coordinate data
%           guidata(hAxes.axis1, A);
%     
%           % Update figure
%           update_fig(hFig, hAxes)
                  
          
    % n (next arm)
      elseif strcmp(key.Key,'n') || strcmp(key.Key,'N')  
          
          % If all frames acquired . . .
          if length(A.arm)==A.numArms
              
              warning('Choose "q" to save data and quit interactive mode');
              
          % If no nans in current arm
          elseif ~max(isnan(A.arm(end).x))
              
              % Make room for a new arm
              A = makeNewArm(A,length(A.arm(end).x));
              
              % Jump back to first frame index
              A.iFrame = 1;
     
              disp(' ');
              disp(' The next arm should be clockwise with respect to the ');
              disp(' one that you just selected.');

              
              % Store coordinate data
              guidata(hAxes.axis1, A);
              
              % Update figure
              update_fig(hFig, hAxes)             
                  
              
          % If some frames not acquired . . .
          else
              
              % Jump to first empty frame
              A.iFrame = find(isnan(A.arm(end).x),1,'first');

              % Issue warning
              warning(['You need to select a coordinate on all frames ' ...
                       'before proceeding to the next arm'])
              
              % Store A data
              guidata(hAxes.axis1, A);     
                   
              % Update figure
              update_fig(hFig, hAxes)
          end
          
     % The following involve changing frame number
     else
         newIndex = [];
         
         % Index of current frame
%          iFrame = find(A.cFrame==A.frames,1,'first');
         
         % RIGHT ARROW/SPACE
         if strcmp(key.Key,'rightarrow') || strcmp(key.Key,'space')
             
             if A.iFrame==length(A.frames)
                 beep
             else
                 % Advance frame
                 newIndex = A.iFrame+1;
             end
             
         % LEFT ARROW
         elseif strcmp(key.Key,'leftarrow')
             
             if A.iFrame==1
                 beep
             else
                 % Reverse frame
                 newIndex = A.iFrame-1;
             end
             
         % UP ARROW
         elseif strcmp(key.Key,'uparrow') 
             
             % Advance frame
             newIndex = min([A.iFrame+A.adInterval length(A.frames)]);
             
         % DOWN ARROW
         elseif strcmp(key.Key,'downarrow')
             

             % Reverse frame
             newIndex = max([1 A.iFrame-A.adInterval]);

             
         % NUMBER
         elseif length(key.Key)==1 && (sum(key.Key==num2str([1:9]))==1)
             
             % Requested relative number
             req_num = str2num(key.Key);
             
             % Start 
             if req_num == 1
                 newIndex = 1;
                 
             % End
             elseif req_num==9
                 newIndex = length(A.frames);
                 
             % Something between
             else
                 % Set new frame
                 newIndex = round((req_num/10)*length(A.frames));
             end
               
         % Key in time (in s)
         elseif strcmp(key.Key,'t')|| strcmp(key.Key,'T')
             
             answer = inputdlg({'Time: min','Time: sec'},'',1,{'0',''});
             
             % If answer provided
             if ~isempty(answer)
                 
                % Convert time into seconds
                timeval = str2num(answer{1})*60 + str2num(answer{2});
                 
                % Translate into frames
                aFrame = ceil(timeval*A.v.FrameRate);
                
                if aFrame<min(A.frames)
                    warning('Requested time happens before analyzed interval');
                    newIndex = [];
                    
                elseif aFrame>max(A.frames)
                    warning('Requested time happens after analyzed interval');
                    newIndex = [];
                    
                else
                    newIndex = find(aFrame==A.frames,1,'first');
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
                
                if aFrame<min(A.frames)
                    warning('Requested time happens before analyzed interval');
                    newIndex = [];
                    
                elseif aFrame>max(A.frames)
                    warning('Requested time happens after analyzed interval');
                    newIndex = [];
                    
                else
                    newIndex = find(aFrame==A.frames,1,'first');
                end
             
             % If no answer
             else
                 newIndex = [];
             end
             
         end

         % If new newframe
         if ~isempty(newIndex)

             % Update frame number
             A.iFrame = newIndex;
             
             % Store coordinate data
             guidata(hAxes.axis1, A);
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
    A = guidata(hAxes.axis2);   
     
    % Collect current coordinate
    C2 = get(hAxes.axis2, 'CurrentPoint');
    C1 = get(hAxes.axis1, 'CurrentPoint');
    
    hold on
    % If in full frame box . . .
    if (C1(1,1)>=0) && (C1(1,1)<=hAxes.axis1.XLim(2)) && ...
       (C1(1,2)>=0) && (C1(1,2)<=hAxes.axis1.YLim(2))     
        
       % Update roi x-coordinate      
       if (C1(1,1)-A.roi(3)/2) < 0
            A.roi(1) = 1;
       elseif (C1(1,1)+A.roi(3)/2) > A.v.Width
           A.roi(1) = A.v.Width-A.roi(3);
       else
           A.roi(1) = C1(1,1)-A.roi(3)/2;
       end
       
       % Update roi y-coordinate 
       if (C1(1,2)-A.roi(4)/2) < 0
           A.roi(2) = 1;
       elseif (C1(1,2)+A.roi(4)/2) > A.v.Height
           A.roi(2) = floor(A.v.Height - A.roi(4));
       else
           A.roi(2) = C1(1,2)-A.roi(4)/2;
       end
                
    % If in roi box . . .
    elseif (C2(1,1)>=0) && (C2(1,1)<=hAxes.axis2.XLim(2)) && ...
           (C2(1,2)>=0) && (C2(1,2)<=hAxes.axis2.YLim(2))     
    
       % Current coordinate
       xCurr = A.roi(1) + C2(1,1);
       yCurr = A.roi(2) + C2(1,2);
  
       % Contact point
       A.arm(end).x(A.iFrame) = xCurr;
       A.arm(end).y(A.iFrame) = yCurr;

    % If cursor not in ROI box . . .
    else            
       warning('Your cursor needs to point to the zoomed window')
      
    end

    set(gcf,'Pointer','arrow')

    % Store coordinate data
    guidata(hAxes.axis1, A);
    
    % Update figure
    update_fig(hFig, hAxes);           
end


%% Update GUI
function update_fig(hFig, hAxes)

    % Activate figure window
    figure(hFig)

    % Load data
    A = guidata(hAxes.axis1);

    % List frame number in title 1
    obj = findobj('tag','title1');
    
    % Define time variables
    currTime = (A.frames(A.iFrame))./A.v.FrameRate;
    currMin = floor(currTime/60);
    currSec = floor(currTime - currMin*60);
    currFrac = floor((currTime - currMin*60-currSec)*100);
    
    % Convert to strings
    currMin  = num2str(currMin);
    currSec  = ['0' num2str(currSec)];
    currFrac = num2str(currFrac);
    
    t_str = ['Full Frame :  Frame ' num2str(A.frames(A.iFrame)) ' (' ...
             currMin ':' currSec(end-1:end) ':' currFrac ')'];
    set(obj,'String',t_str);

    % List mode in title 2
    obj = findobj('tag','title2');
    set(obj,'String','ROI');

    % Current frame
    cFrame = A.frames(A.iFrame);
    
    % Get all data
    [allTip,currTip] = getArmData(A,A.iFrame);

    % Read input video frame
    frame = getFrame(A.vid_path,A.v,cFrame,A.imInvert,'rgb');

    % Center of ROI (Full frame FOR)
    xCntr = A.roi(1) + A.roi(3)/2;
    yCntr = A.roi(2) + A.roi(4)/2;

    
    % FULL FRAME WINDOW --------------------------------------------------
    
    % Display full video frame
    delete(hAxes.axis1.Children)
    showFrameOnAxis(hAxes.axis1, frame, 0);
    
    % Marker size
    mSize = 100;
    
    % Hold on Axis1
    set(hAxes.axis1,'NextPlot','Add')

    % Set A.roi square on full frame
    x1 = A.roi(1); x2 = A.roi(1)+A.roi(3);
    y1 = A.roi(2); y2 = A.roi(2)+A.roi(4);
    plot([x1 x2 x2 x1 x1], [y1 y1 y2 y2 y1], 'w-',...
        mean([x1 x2]).*[1 1],[y1 y2], 'w:', ...
        [x1 x2], mean([y1 y2]).*[1 1], 'w:',...
        'Parent',hAxes.axis1);
    
    % Plot prior tip points
    h2 = scatter(allTip.x,allTip.y,'Parent',hAxes.axis1,...
        'MarkerEdgeColor',[1 1 1],'MarkerfaceColor','none','SizeData',mSize);
    
    % Plot current tip point
    h2 = scatter(currTip.x,currTip.y,'Parent',hAxes.axis1,...
        'MarkerEdgeColor',[0 1 0],'MarkerfaceColor','none','SizeData',mSize);
    

    % Hold off Axis1
    set(hAxes.axis1,'NextPlot','Replace');

    
    % ZOOMED WINDOW -------------------------------------------------------
    
    % Get cropped image
    im2 = imcrop(frame, A.roi);

    if size(im2,2)/size(im2,1) ~= A.roi(3)./A.roi(4)
        beep;beep;
        warning(['DO NOT SELECT POINTS : Cropped frame does not have ' ...
                 'same dimensions as the roi. Choose another ROI.'])
    end
    
    % Display cropped image
    delete(hAxes.axis2.Children);
    showFrameOnAxis(hAxes.axis2, im2, 0);

    % Limits for x and y axes
    xL = hAxes.axis2.XLim;
    yL = hAxes.axis2.YLim;

    % Hold on Axis 2
    set(hAxes.axis2,'NextPlot','Add');
    
    % Plot prior tip points
    h2 = scatter(allTip.x-A.roi(1),allTip.y-A.roi(2),'Parent',hAxes.axis2,...
        'MarkerEdgeColor',[1 1 1],'MarkerfaceColor','none','SizeData',500);
    
    % Plot current tip point
    h2 = scatter(currTip.x-A.roi(1),currTip.y-A.roi(2),'Parent',hAxes.axis2,...
        'MarkerEdgeColor',[0 1 0],'MarkerfaceColor','none','SizeData',500);
    

    % Hold off Axis 2
    set(hAxes.axis2,'NextPlot','Replace');

    % Activate figure
    figure(hFig)

    % Activate ROI panel
    set(hAxes.axis2,'Selected','on')
    
end

function A = makeNewArm(A,dataLen)
% Makes new place for foot data

    % Set index for adding data
    if isempty(A) || ~isfield(A,'arm')
        iAdd = 1;
    else
        iAdd = length(A.arm) + 1;
    end

    % Column of nans
    nanCol = nan(dataLen,1);

    % Tip point
    A.arm(iAdd).x = nanCol;
    A.arm(iAdd).y = nanCol;
    
end

function [allTip,currTip] = getArmData(A,iFrame)
% Returns data for visualization

    % Initialize containers
    allTip.x   = [];   allTip.y   = [];

    % Loop thru feet, prior to current
    for i = 1:length(A.arm)-1
        
        
        % If there is a coordinate for current frame . . .
        if ~isnan(A.arm(i).x(iFrame))
            
            % Add coordinates
            allTip.x    = [allTip.x; A.arm(i).x(iFrame)];
            allTip.y    = [allTip.y; A.arm(i).y(iFrame)];
        end
    end
    
    % Current tip point
    currTip.x    = A.arm(end).x(iFrame);
    currTip.y    = A.arm(end).y(iFrame); 
end
