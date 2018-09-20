function H = videoGUI(vid_path,v,frames,imInvert,acqMode,varargin)
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
    
    if length(varargin)>2
        H = varargin{3};
    else
        H = [];
    end
    
%     % Starting values
%     if length(varargin)>3
%         xStart = varargin{3};
%         yStart = varargin{4};
%     else
%         xStart = nan(length(frames),1);
%         yStart = nan(length(frames),1);
%     end
else
    error('Do not recognize acqMode');
end


%% Default parameters

% Container of handle data
%H = [];

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
roi = [size(im0,2)/2-r(1)/2 size(im0,1)/2-r(1)/2 2*r 2*r];


% Display options
disp(' ')
disp('COMMANDS  ------------------------------------------------ ')
disp('    left click  : select position in current frame')
disp('    delete      : delete current tubefoot data')
disp('    right arrow : advance one frame')
disp('    left arrow  : go back one frame')
disp('    up arrow    : advance multiple frames')
disp('    down arrow  : go back multiple frames')
disp('    j           : set number of frames to jump for up/down arrows')
disp('    1 - 9       : jump to relative position in video segment')
disp('    f           : Jump to frame number')
disp('    t           : Jump to time')
disp('    +           : zoom in')
disp('    -           : zoom out')
%disp('    i           : invert image')
%disp('    d           : delete all points after current in interval')
%disp('    a           : delete all points in interval')
disp('    s           : Save current tubefoot coordinates')
disp('    q           : quit interaction mode')
disp(' ')
 
% Create figure window
[hFig, hAxes] = createFigureAndAxes;

% Add data to figure
[hFig, hAxes] = putData(hFig, hAxes);

% Wait for completion of interactive mode
waitfor(hFig)

% If current tube foot is incomplete . . .
if (sum(~isnan(H.ft(end).xBase))==0 || ...
    H.ft(end).choseContact==0 || ...
    H.ft(end).choseRelease==0)

    % And still on the first tube foot . . .
    if length(H.ft)==1
        return
       
    % If beyond first . . .
    else        
        % Delete data for current tube foot
        H.ft(end).xBase        = nan(length(H.ft(end).xBase),1);
        H.ft(end).yBase        = nan(length(H.ft(end).yBase),1);
        H.ft(end).xTip         = nan(length(H.ft(end).xTip),1);
        H.ft(end).yTip         = nan(length(H.ft(end).yTip),1);
        H.ft(end).choseContact = 0;
        H.ft(end).choseRelease = 0;
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


%% Deposite data into figure

    function [hFig, hAxes] = putData(hFig, hAxes)
        
        % Initialize data repository
        H = makeNewFoot(H,length(frames));
        
        % Store inital values into axis1
        H.roi = roi;
        H.iFrame = 1;
        H.frames = frames;
        H.vid_path = vid_path;
        H.v = v;
        H.clr = mClr;
%         H.x = xStart;
%         H.y = yStart;
        H.imInvert = imInvert;
        H.adInterval = 15;
        H.iMode = 'b';

        % Store coordinate data
        guidata(hAxes.axis1, H);
        
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
     H = guidata(hAxes.axis2);
     
     % QUIT ('q')
      if strcmp(key.Key,'Q') || strcmp(key.Key,'q')
         
         % If there's any unsaved data . . .
         if sum(~isnan(H.ft(end).xTip))>0 || sum(~isnan(H.ft(end).xBase))>0
             
             bName = questdlg(['Are you sure you want to quit without ' ...
                 'saving the current tube foot data?'],'','Yes','No, Cancel','Yes');
             
             if strcmp(bName,'Yes')
                assignin('caller','H',H);
                closefig = 1;
             end
             
         % If there's no data    
         else
             if length(H.ft)==1
                 H = rmfield(H,'ft');
             else
                 H.ft = H.ft(1:end-1);
             end
             
             assignin('caller','H',H);
             closefig = 1;
         end
         
         
     % DELETE
     elseif strcmp(key.Key,'backspace')
         
       % Delete data
       H.ft(end).xBase        = nan(length(H.ft(end).xBase),1);
       H.ft(end).yBase        = nan(length(H.ft(end).yBase),1);
       H.ft(end).xTip         = nan(length(H.ft(end).xTip),1);
       H.ft(end).yTip         = nan(length(H.ft(end).yTip),1);
       H.ft(end).choseContact = 0;
       H.ft(end).choseRelease = 0;
       
       % Store coordinate data
       guidata(hAxes.axis1, H);
       
       % Update figure
       update_fig(hFig, hAxes)
     
     % PLUS (zoom in)
     elseif strcmp(key.Key,'hyphen') || strcmp(key.Key,'minus')
         
         xCntr = H.roi(1) + H.roi(3)/2;
         yCntr = H.roi(2) + H.roi(4)/2; 
         
         % Enlarge roi
         H.roi(3) = 1.5*H.roi(3);
         H.roi(4) = 1.5*H.roi(4);
         H.roi(1) = xCntr - H.roi(3)/2;
         H.roi(2) = yCntr - H.roi(4)/2;
         
        % Store coordinate data
        guidata(hAxes.axis1, H);
       
       % Update figure
       update_fig(hFig, hAxes)
       
     % PLUS (zoom out)
      elseif strcmp(key.Key,'equal') || strcmp(key.Key,'plus')
          
         
         xCntr = H.roi(1) + H.roi(3)/2;
         yCntr = H.roi(2) + H.roi(4)/2; 
         
         % Enlarge roi
         H.roi(3) = H.roi(3)/1.5;
         H.roi(4) = H.roi(4)/1.5;
         H.roi(1) = xCntr - H.roi(3)/2;
         H.roi(2) = yCntr - H.roi(4)/2;
         
        % Store coordinate data
        guidata(hAxes.axis1, H);
       
       % Update figure
       update_fig(hFig, hAxes)
       
%      % D (delete data)
%      elseif strcmp(key.Key,'d') || strcmp(key.Key,'D')
%          
%          but = questdlg(['You want to delete all points after this frame?'],...
%              'ALERT!','Yes','No','Yes');
%          
%          if strcmp(but,'Yes')           
%              H.x((H.iFrame+1):end) = nan;
%              H.y((H.iFrame+1):end) = nan;
%          end
%          
%          % Store coordinate data
%          guidata(hAxes.axis1, H);
%     
%          % Update figure
%          update_fig(hFig, hAxes)
%    
%       % I (invert image)
%       elseif strcmp(key.Key,'i') || strcmp(key.Key,'I')
%           
%           % Invert image
%           H.imInvert = abs(H.imInvert - 1);
%           
%           % Store coordinate data
%           guidata(hAxes.axis1, H);
%           
%           % Update figure
%           update_fig(hFig, hAxes)

      % S (set interval)
      elseif strcmp(key.Key,'j') || strcmp(key.Key,'J')
          
          % Prompt for interval
          answer = inputdlg({'Frame interval'},'Set frame interval',1,...
                            {num2str(H.adInterval)});                   
          if isempty(answer)
              return
          end
          
          % Invert image
          H.adInterval = str2num(answer{1});
          
          % Store coordinate data
          guidata(hAxes.axis1, H);
    
          % Update figure
          update_fig(hFig, hAxes)
          
    % N (start new foot)
      elseif strcmp(key.Key,'s') || strcmp(key.Key,'S')  
          
          disp(' ');
          disp(' Tube feet are numbered 1, 2, 3 . . . with 1 in the first proximal position');
          disp(' and A or B, with A in the more clockwise position and B counter-clockwise');
          disp(' ');
          
          answer = inputdlg({'Tube foot number (1-12)','Tube foot letter (A or B)'},...
                   'Data for tube foot just finished',1,{'',''});
          
          if strcmp(answer{2},'A') || strcmp(answer{2},'a')
              
              H.ft(end).footLet = 'A';
              
          elseif strcmp(answer{2},'B') || strcmp(answer{2},'b')
              
              H.ft(end).footLet = 'B';
              
          end
          
          % Address for the tube foot
          H.ft(end).footNum = str2num(answer{1});
          
          if isnan(H.ft(end).footLet) || (H.ft(end).footNum<1) || (H.ft(end).footNum>20)             
              
              warning('Do not recognize tube foot address -- try again')
              
          else

              % If we have all three points . . .
              if H.ft(end).choseContact==1 && H.ft(end).choseRelease==1 && ...
                      sum(~isnan(H.ft(end).xBase))~=0
                  
                  % Make room for a new foot
                  H = makeNewFoot(H,length(H.ft(end).xBase));
                  
              else
                  warning(['Need to select contact, release, and base points ' ...
                      'before starting a new foot']);
              end
              
              % Store coordinate data
              guidata(hAxes.axis1, H);
              
              % Update figure
              update_fig(hFig, hAxes)
          end
         
     % The following involve changing frame number
     else
         newIndex = [];
         
         % Index of current frame
%          iFrame = find(H.cFrame==H.frames,1,'first');
         
         % RIGHT ARROW/SPACE
         if strcmp(key.Key,'rightarrow') || strcmp(key.Key,'space')
             
             if H.iFrame==length(H.frames)
                 beep
             else
                 % Advance frame
                 newIndex = H.iFrame+1;
             end
             
         % LEFT ARROW
         elseif strcmp(key.Key,'leftarrow')
             
             if H.iFrame==1
                 beep
             else
                 % Reverse frame
                 newIndex = H.iFrame-1;
             end
             
         % UP ARROW
         elseif strcmp(key.Key,'uparrow') 
             
             % Advance frame
             newIndex = min([H.iFrame+H.adInterval length(H.frames)]);
             
         % DOWN ARROW
         elseif strcmp(key.Key,'downarrow')
             

             % Reverse frame
             newIndex = max([1 H.iFrame-H.adInterval]);

             
         % NUMBER
         elseif length(key.Key)==1 && (sum(key.Key==num2str([1:9]))==1)
             
             % Requested relative number
             req_num = str2num(key.Key);
             
             % Start 
             if req_num == 1
                 newIndex = 1;
                 
             % End
             elseif req_num==9
                 newIndex = length(H.frames);
                 
             % Something between
             else
                 % Set new frame
                 newIndex = round((req_num/10)*length(H.frames));
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
                
                if aFrame<min(H.frames)
                    warning('Requested time happens before analyzed interval');
                    newIndex = [];
                    
                elseif aFrame>max(H.frames)
                    warning('Requested time happens after analyzed interval');
                    newIndex = [];
                    
                else
                    newIndex = find(aFrame==H.frames,1,'first');
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
                aFrame = round(answer{1});
                
                if aFrame<min(H.frames)
                    warning('Requested time happens before analyzed interval');
                    newIndex = [];
                    
                elseif aFrame>max(H.frames)
                    warning('Requested time happens after analyzed interval');
                    newIndex = [];
                    
                else
                    newIndex = find(aFrame==H.frames,1,'first');
                end
             
             % If no answer
             else
                 newIndex = [];
             end
             
         end

         % If new newframe
         if ~isempty(newIndex)

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
    
    hold on
     % If in full frame box . . .
    if (C1(1,1)>=0) && (C1(1,1)<=hAxes.axis1.XLim(2)) && ...
       (C1(1,2)>=0) && (C1(1,2)<=hAxes.axis1.YLim(2))     
        
       % Update roi
       H.roi(1) = C1(1,1)-H.roi(3)/2;
       H.roi(2) = C1(1,2)-H.roi(4)/2;
                
    % If in roi box . . .
    elseif (C2(1,1)>=0) && (C2(1,1)<=hAxes.axis2.XLim(2)) && ...
           (C2(1,2)>=0) && (C2(1,2)<=hAxes.axis2.YLim(2))     
    
       % Current coordinate
       xCurr = H.roi(1) + C2(1,1);
       yCurr = H.roi(2) + C2(1,2);
       
       % Prompt for type
       answer = questdlg('What type of point ?','','Base','Tip: first contact', ...
                         'Tip: release','Base');
           
       % Store base point
       if strcmp(answer,'Base')
           
           % Base point
           H.ft(end).xBase(H.iFrame) = xCurr;
           H.ft(end).yBase(H.iFrame) = yCurr;
           
       % Store tip at time of contact
       elseif strcmp(answer,'Tip: first contact')
           
           % Index of release 
           iRelease = find(~isnan(H.ft(end).xTip),1,'last');
                      
           % Determine contact coordinates and interval 
           if isempty(iRelease)
               idx       = H.iFrame:length(H.ft(end).xTip);
               xContact  = xCurr;
               yContact  = yCurr;
               iNan      = [];
           else
               idx        = H.iFrame:iRelease;
               xContact   = H.ft(end).xTip(iRelease);
               yContact   = H.ft(end).yTip(iRelease);
               iNan       = 1:(H.iFrame-1);
           end
           
           % Contact point
           H.ft(end).xTip(idx) = xContact;
           H.ft(end).yTip(idx) = yContact;
           
           % Fill after with nans
           H.ft(end).xTip(iNan) = nan;
           H.ft(end).yTip(iNan) = nan;  
           
           H.ft(end).choseContact = 1;
           
       % Store tip at time of release
       elseif strcmp(answer,'Tip: release')
           
           % Index of contact
           iContact = find(~isnan(H.ft(end).xTip),1,'first');
           
           % Determine contact coordinates and interval 
           if isempty(iContact)
               idx       = 1:H.iFrame;
               xContact  = xCurr;
               yContact  = yCurr;
               iNan      = [];
           else
               idx        = iContact:H.iFrame;
               xContact   = H.ft(end).xTip(iContact);
               yContact   = H.ft(end).yTip(iContact);
               iNan       = (H.iFrame+1):length(H.ft(end).xTip);
           end
           
           % Contact point
           H.ft(end).xTip(idx) = xContact;
           H.ft(end).yTip(idx) = yContact;  
           
           % Fill after with nans
           H.ft(end).xTip(iNan) = nan;
           H.ft(end).yTip(iNan) = nan;  
           
           H.ft(end).choseRelease = 1;
       end
       
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

    % List frame number in title 1
    obj = findobj('tag','title1');
    
    % Define time variables
    currTime = (H.frames(H.iFrame))./H.v.FrameRate;
    currMin = floor(currTime/60);
    currSec = floor(currTime - currMin*60);
    currFrac = floor((currTime - currMin*60-currSec)*100);
    
    % Convert to strings
    currMin  = num2str(currMin);
    currSec  = ['0' num2str(currSec)];
    currFrac = num2str(currFrac);
    
    t_str = ['Full Frame :  Frame ' num2str(H.frames(H.iFrame)) ' (' ...
        currMin ':' currSec(end-1:end) ':' currFrac ')'];
    set(obj,'String',t_str);

    % List mode in title 2
    obj = findobj('tag','title2');
    set(obj,'String','ROI');

    % Current frame
    cFrame = H.frames(H.iFrame);
    
    % Get all data
    [allBase,allTip,currBase,currTip] = getFootData(H,cFrame);

    % Read input video frame
    frame = getFrame(H.vid_path,H.v,cFrame,H.imInvert,'rgb');

    % Display full video frame
    delete(hAxes.axis1.Children)
    showFrameOnAxis(hAxes.axis1, frame, 0);

%         % Get window dimensions (full frame FOR)
%         winWidth1  = size(hAxes.axis1.Children(end).CData,2);
%         winHeight1 = size(hAxes.axis1.Children(end).CData,1);
%         
%         
    % Get coordinates for the ROI (Full frame FOR)
    %rect = get_val;

    % Center of ROI (Full frame FOR)
    xCntr = H.roi(1) + H.roi(3)/2;
    yCntr = H.roi(2) + H.roi(4)/2;


    % Hold on Axis1
    set(hAxes.axis1,'NextPlot','Add')

    % Set H.roi square on full frame
    x1 = H.roi(1); x2 = H.roi(1)+H.roi(3);
    y1 = H.roi(2); y2 = H.roi(2)+H.roi(4);
    plot([x1 x2 x2 x1 x1], [y1 y1 y2 y2 y1], 'w-',...
        mean([x1 x2]).*[1 1],[y1 y2], 'w:', ...
        [x1 x2], mean([y1 y2]).*[1 1], 'w:',...
        'Parent',hAxes.axis1);

    % Hold off Axis1
    set(hAxes.axis1,'NextPlot','Replace');

    % Get cropped image
    im2 = imcrop(frame, H.roi);

    % Display cropped image
    delete(hAxes.axis2.Children);
    showFrameOnAxis(hAxes.axis2, im2, 0);

    % Limits for x and y axes
    xL = hAxes.axis2.XLim;
    yL = hAxes.axis2.YLim;

    % Hold on Axis 2
    set(hAxes.axis2,'NextPlot','Add');

    % Plot prior base points
%     h2 = plot(allBase.x-H.roi(1),allBase.y-H.roi(2),'w+',...
%               'Parent',hAxes.axis2);
    h2 = scatter(allBase.x-H.roi(1),allBase.y-H.roi(2),'g+','Parent',hAxes.axis2,...
        'MarkerEdgeColor',[1 1 1],'MarkerfaceColor','none','SizeData',300);   
    
    % Plot prior tip points
    h2 = scatter(allTip.x-H.roi(1),allTip.y-H.roi(2),'Parent',hAxes.axis2,...
        'MarkerEdgeColor',[1 1 1],'MarkerfaceColor','none','SizeData',500);
    
    % Plot current tip point
    h2 = scatter(currTip.x-H.roi(1),currTip.y-H.roi(2),'Parent',hAxes.axis2,...
        'MarkerEdgeColor',[0 1 0],'MarkerfaceColor','none','SizeData',500);
    
     % Plot current base point
    h2 = scatter(currBase.x-H.roi(1),currBase.y-H.roi(2),'g+',...
              'Parent',hAxes.axis2,'MarkerEdgeColor',[0 1 0],...
              'MarkerfaceColor','none','SizeData',300); 

    % Hold off Axis 2
    set(hAxes.axis2,'NextPlot','Replace');

    % Activate figure
    figure(hFig)

    % Activate ROI panel
    set(hAxes.axis2,'Selected','on')
end

function H = makeNewFoot(H,dataLen)
% Makes new place for foot data

    % Set index for adding data
    if isempty(H)
        iAdd = 1;
    else
        iAdd = length(H.ft) + 1;
    end

    % Column of nans
    nanCol = nan(dataLen,1);

    % Base point
    H.ft(iAdd).xBase = nanCol;
    H.ft(iAdd).yBase = nanCol;

    % Tip point
    H.ft(iAdd).xTip = nanCol;
    H.ft(iAdd).yTip = nanCol;

    % Logical about whether contact/release yet chosen
    H.ft(iAdd).choseContact = 0;
    H.ft(iAdd).choseRelease = 0;
    
    % Address for the tube foot
    H.ft(iAdd).footNum = nan;
    H.ft(iAdd).footLet = nan;
end

function [allBase,allTip,currBase,currTip] = getFootData(H,cFrame)
% Returns data for visualization

    % Initialize containers
    allBase.x  = [];   allBase.y  = [];
    allTip.x   = [];   allTip.y   = [];

    % Loop thru feet, prior to current
    for i = 1:length(H.ft)-1
        
        % If there is a coordinate for current frame . . .
        if ~isnan(H.ft(i).xBase(cFrame))
            
            % Add coordinates
            allBase.x   = [allBase.x; H.ft(i).xBase(cFrame)];
            allBase.y   = [allBase.y; H.ft(i).yBase(cFrame)];
        end
        
        % If there is a coordinate for current frame . . .
        if ~isnan(H.ft(i).xTip(cFrame))
            
            % Add coordinates
            allTip.x    = [allTip.x; H.ft(i).xTip(cFrame)];
            allTip.y    = [allTip.y; H.ft(i).yTip(cFrame)];
        end
    end
    
    % Current base point
    currBase.x    = H.ft(end).xBase(cFrame);
    currBase.y    = H.ft(end).yBase(cFrame); 
    
    % Current tip point
    currTip.x    = H.ft(end).xTip(cFrame);
    currTip.y    = H.ft(end).yTip(cFrame); 
end
