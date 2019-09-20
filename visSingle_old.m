function visSingle_old
(Orientation,SSnum,seqnum,makeMovie)
% Visualize the results of a single sequence
% Orientation - indicates the wall orientation ('h','v','u')
% SSnum       - number of seastar
% seqnum      - sequence number



% Paths
paths = givePaths;

%Camera view (for 2D analysis)
camView = 'canon';


%% Batch mode

if nargin==0
    
    % Catalog list
    cList = catDatafiles(paths.data,camView,'Bundled Data');
    
    % Loop thru sequences
    for i = 1:length(cList.fName)
        
        % Update status
        disp(['Running visSingle on: ' cList.orient(i) ', ' ...
               num2str(cList.indiv(i)) ', ' num2str(cList.seq(i))])
        
        % Run code
        visSingle(cList.orient(i),cList.indiv(i),cList.seq(i),0);
    end
    
    return
end

%% Deal with inputs

if nargin >=3
    
    txt1 = ['0' num2str(SSnum)];
    txt2 = ['0' num2str(seqnum)];
    
    % Movie name
    movName = [Orientation 'SS' txt1(end-1:end) '_s' txt2(end-1:end) '_ana'];
    
    % Title name (for graphs)
    tName = [Orientation ' SS' txt1(end-1:end) ' s' txt2(end-1:end)];
    
    if strcmp(Orientation,'h')
        Orientation = 'Horizontal';
        
    elseif strcmp(Orientation,'v')
        Orientation = 'Vertical';
        
    elseif strcmp(Orientation,'u')    
        Orientation = 'Upside-down';
        
    else
        error('Do not recorgnize the orientation requested');
    end
    
    % Make string out of seastar number
    SSnum = ['0' num2str(SSnum)];
    SSnum = ['SS' SSnum(end-1:end)];
    
    % Make string out of sequence number
    seqnum = ['0' num2str(seqnum)];
    seqnum = seqnum(end-1:end);
    
    good.path  =  [Orientation filesep SSnum filesep 'canon']; 
    good.fName = ['s' seqnum];
    
    if nargin < 4
        makeMovie = 0;
    end
    
    
    clear Orientation SSnum seqnum txt1 txt2
    
else
    error('You need to request the sequence')
    
end

%% Code execution

% Makes an animation of the data
do.aniData = 0;

% Draw trajectory
do.drawTraj = 1;

% Draw gait diagram
do.gaitDiagram = 0;

% Map of tube feet
do.feetMap = 1;


%% Parameters

% Initial duration for identifying a trajectory-based coordinate system (s)
initialDur = 15;

% Dock figure windows
set(0,'DefaultFigureWindowStyle','docked')

% Color of seastar body
clr{1} = 0.7.*[1 1 1];

% Color of tube feet
clr{2} = [0 0.75  0.3];

% Center point of body
clr{3} = .3.*[1 1 1];

if makeMovie
    
    mov_path = [paths.data filesep 'Data movies' filesep movName];
    
    vOut = VideoWriter([mov_path '.mp4'],'MPEG-4');
    vOut.Quality = 50;
    open(vOut)
end


%% Load data

% Load bundled 2D data ('S')
load([paths.data filesep good.path filesep ...
      good.fName filesep 'Bundled Data.mat']);


%% Visualize trajectory

if do.drawTraj
    % Time interval for drawing body
    tInt = 30;
    
    % Times for drawing body
    %tDraw = S.t(1):tInt:S.t(end);
    tDraw = [S.t(1) S.t(end)];
    
    figure
    
    for i = 1:length(tDraw)
        
        % Arm points in trajectory FOR
        %armPntsT = transCoord2d('xax G2L',S.pStart,S.pEnd,[S.armT(i).x, S.armT(i).y]);
        
        % Time difference btwen current and all times
        tDiff = abs(tDraw(i)-S.t);
        
        % Index is closest time to desired
        idx = find(tDiff==min(tDiff),1,'first');
        
        % Points for whole body
        starPtsT = makeStar(S.armT(idx).x, S.armT(idx).y, S.xCntrT(idx),S.yCntrT(idx));
        
        % Draw body at interval
        h = fill(starPtsT(:,1),starPtsT(:,2),clr{1},'FaceAlpha',0.3,...
            'EdgeColor','none');
        
        hold on
        
        clear starPntsT idx tDiff
    end
    
    % Plot center points
    h = scatter(S.xCntrT,S.yCntrT,'SizeData',10,'MarkerEdgeColor','none', ...
        'MarkerFaceColor',0.5.*[1 1 1]);
    
    line([min(S.xCntrT) max(S.xCntrT)],[0 0],'Color','k')
    axis equal
    set(gca,'YColor','none','XColor','none')
    
    clear tInt
end

%% Gait diagram

if do.gaitDiagram
    
    % Make plot
    gaitPlot(S,tName)
    
end

%% Direction of tube feet

if do.feetMap
    
    feetMap(S,tName);
end



%% Data animation

if do.aniData
    figure
    
    
    % Diameter of sea star
    aLen = max(hypot(S.armL.x,S.armL.y));
    
    % Limits of x-axis
    xL = 1.4*[min(S.xCntrT)-aLen max(S.xCntrT)+aLen];
    
    % Step thru time
    for i = 1:length(S.t)
        
        h = scatter(S.xCntrT,S.yCntrT,'MarkerEdgeColor','none',...
            'MarkerFaceColor',clr{3},'Sizedata',10);
        a1 = gca;
        hold on
        axis equal
        
        % Arm points in trajectory FOR
        armPntsT = transCoord2d('xax G2L',pStart,pEnd,[S.arm(i).x, S.arm(i).y]);
        
        % Points for whole body
        starPtsT = makeStar(armPntsT(:,1), armPntsT(:,2), CntrPtsT(i,1), ...
                            CntrPtsT(i,2));
        
        % Cntr points in trajectory FOR
       % CntrPtsT = transCoord2d('xax G2L',pStart,pEnd,CntrPnts);
        
        h = fill(starPtsT(:,1),starPtsT(:,2),clr{1},'FaceAlpha',0.3,...
            'EdgeColor','none');
        
        % Center point
        scatter(CntrPtsT(i,1),CntrPtsT(i,2),'MarkerEdgeColor','r',...
               'SizeData',50')
        
        % Step thru tube feet
        for j = 1:length(S.ft)
            if ~isnan(S.ft(j).xBase(i))
                
                % Current coordinates
                tip  = [S.ft(j).xTip(i) S.ft(j).yTip(i)];
                base = [S.ft(j).xBase(i) S.ft(j).yBase(i)]; 
                
                % Transform in trajectory FOR
                tipT = transCoord2d('xax G2L',pStart,pEnd,tip);
                baseT = transCoord2d('xax G2L',pStart,pEnd,base);
                
                % Plot
                line([tipT(1) baseT(1)],[tipT(2) baseT(2)],...
                    'Color',[clr{2} 0.3],'LineWidth',5,'Parent',a1);
                scatter(tipT(1),tipT(2),...
                    'MarkerFaceColor',clr{2},...
                    'MarkerEdgeColor',clr{2},'LineWidth',2,...
                    'MarkerEdgeAlpha',0.8,'MarkerFaceAlpha',0.8,...
                    'SizeData',100,'Parent',a1);
            end
            
        end 

           
        % Arm points
%         scatter(armPntsT(:,1),armPntsT(:,2),'MarkerEdgeColor','g',...
%                'SizeData',50')
        
        
%         armsT = transCoord2d('xax G2L',pStart,pEnd,starPts);
%         line(armsT(:,1),armsT(:,2))
        
        title(['t = ' num2str(S.t(i))])
        
        hold off
        
        axis square
        xlim(xL);
        ylim([-range(xL)/2 range(xL)/2]);
        set(a1,'YColor','none','XColor','none')
        drawnow
        pause(0.01)
        
        
        if makeMovie
            % Get image
            imFrame = getframe(gcf);
            
            % Write frame
            writeVideo(vOut,imFrame);
        end

    end
    
    if makeMovie
       close(vOut) 
    end
    
end


function starPts = makeStar(xA,yA,xC,yC)
% Create periperal shape of seastar body

nArms = length(xA);

% xA = [xA;xA(1)];
% yA = [yA;yA(1)];

nPts = 500;

% Radial position of arms, sorted
[Atheta,idx] = sort(wrapTo2Pi(atan2(yA-yC,xA-xC)));

% Arm lengths
ALen = mean(hypot(xA(idx)-xC,yA(idx)-yC));

% Diameter of cental disc
D = mean(ALen)*0.5;

% Tack on starting arm
Atheta = unwrap([Atheta; Atheta(1)]);
%ALen   = [ALen; ALen(1)];

tVal0 = linspace(0,pi,round(nPts/nArms))';
tVal1 = linspace(pi,2*pi,round(nPts/nArms))';

xP = []; yP = [];

% Loop thru arms
for i = 1:nArms
    
    % Radial positions of current and next arms
    %Atheta0 = atan2(yA(i)-yC,xA(i)-xC);
    %Atheta1 = atan2(yA(i+1)-yC,xA(i+1)-xC);
    
    % Arm lengths of current and next arms
    %ALen0 = hypot(xA(i)-xC,yA(i)-yC);
    %ALen1 = hypot(xA(i+1)-xC,yA(i+1)-yC);
    
    % Radial position values wrapTo2Pi
    
    tmp1 = min([Atheta(i),mean([Atheta(i) Atheta(i+1)])]);
    tmp2 = max([Atheta(i),mean([Atheta(i) Atheta(i+1)])]);    
    theta0 = linspace(tmp1,tmp2,round(nPts/nArms))';
    
    
    tmp1 = min([Atheta(i+1),mean([Atheta(i) Atheta(i+1)])]);
    tmp2 = max([Atheta(i+1),mean([Atheta(i) Atheta(i+1)])]); 
    theta1 = linspace(tmp1,tmp2,round(nPts/nArms))';
      
    rVal0 = (ALen-D).*(cos(tVal0)+1)/2+D;
    rVal1 = (ALen-D).*(cos(tVal1)+1)/2+D;
    
    [x0,y0] = pol2cart(theta0,rVal0);
    [x1,y1] = pol2cart(theta1,rVal1);
    
    xP = [xP; x0; x1];
    yP = [yP; y0; y1];
    
    if 0
        
       subplot(2,2,1)
       plot(tVal0,rVal0,'-',tVal1,rVal1,'-')
       %plot(theta0,rVal0,'-',theta1,rVal1,'-')
       hold on
       
       subplot(2,2,3)
       plot(tVal0,theta0,'-',tVal1,theta1,'-')
       hold on
       
       subplot(2,2,[2 4])
       plot(x0,y0,'-',x1,y1,'-')
       axis equal
       hold on
       
       ttt=3;
    end
end

starPts = [xP+xC yP+yC];


function gaitPlot(S,ttext)

% Number of arms
numArms = 6;

% Height of all bars for one arm
hBar = 1;

% Maximum foot number used
maxFt = 6;

% Color of the seastar body
sClr = 0.8.*[1 1 1];

% Make figure window
f = figure('Units','normalized');

% Transparency of bars
fAlpha = 1;

% Initialize index
j = 1;

% Loop thru feet
for i = 1:length(S.ft)
     % If there are coodrinates . . .
     if sum(~isnan(S.ft(i).xBase))>0 && sum(~isnan(S.ft(i).xTip))>0

         % Get points for current foot in trajectory FOR
         xTip(:,j)       = S.ft(i).xTip;
         armNum(:,j)     = S.ft(i).armNum;
         ftNum(:,j)      = S.ft(i).footNum;
         
         % Advance index
         j = j + 1;     
     end 
end

% Width of body along x-axis
normLen = range(S.armLT.x);

% Mean x and y arm coordinates
meanArmx = 0.8*hBar .* S.armLT.x./normLen;
meanArmy = 0.8*hBar .* S.armLT.y./normLen;

% Coordinate for seastar body
starPts = makeStar(meanArmx, meanArmy, 0, 0);

% Check on foot number
if max(ftNum(:))>maxFt
    error('actual foot number greater than max assumed');
end

% Pull default colormap
cmap = colormap(f,'parula');
cval = linspace(0,1.2,size(cmap,1));

% Custom colormap by interpolation
for i = 1:maxFt
    cmap2(i,1) = interp1(cval,cmap(:,1),i./maxFt);
    cmap2(i,2) = interp1(cval,cmap(:,2),i./maxFt);
    cmap2(i,3) = interp1(cval,cmap(:,3),i./maxFt);
end

clear cmap cval

% Make legend for tube foot colors
if 0
    
    figure
    subplot(4,5,1)
    for i = 1:maxFt
        sSize = 0.25;
        
        yV = [i-sSize i+sSize i+sSize i-sSize i-sSize];
        xV = [-sSize -sSize sSize sSize -sSize];
        
        h = fill(xV, yV, cmap2(i,:));
        set(h,'EdgeColor','none');
        
        h = text(1.1*sSize,i,num2str(i));
        set(h,'FontSize',16,'HorizontalAlignment','left','VerticalAlignment','middle')
        
        hold on
    end
    set(gca,'XColor','none','YColor','none');
    %axis equal
    hold off
    title('Tube foot number')
    
    clear sSize h xV yV  
end

% Calculate speed
spd = hypot(diff(smooth(S.xCntr)),diff(smooth(S.yCntr)))./diff(S.t)';
spd = [spd(1); spd];

% Calculate heading
for i = 1:length(S.armT)
    theta(i,1) = 180/pi*atan2(S.armT(i).y(1)-S.yCntrT(i),S.armT(i).x(1)-S.xCntrT(i));
end

% Plot speed and heading
figure(f);
subplot(4,5,2:5)
[ax1,h(1),h(2)] = plotyy(S.t,spd,S.t,theta);
ylabel(ax1(1),'Speed (pix/s)')
ylabel(ax1(2),'Heading (deg)')
%xlabel('Time s)');
set(h(1),'Color','k','LineWidth',1.5)
set(ax1(1),'YColor','k')
set(h(2),'Color',0.7.*[1 1 1],'LineWidth',1.5)
set(ax1,'TickDir','out')
set(ax1(2),'YColor',0.5.*[1 1 1])
title(ttext)

% Get x-limits
xL = xlim;

% Create axes
h1 = subplot(4,5,6:5:16);
h2 = subplot(4,5,[7:10 12:15 17:20]);

% Loop thru arms and plot seastar bodies
for i = 1:numArms
    axes(h1)
    h = fill(starPts(:,1),starPts(:,2)+i-0.5,sClr,'EdgeColor','none');
    axis equal
    hold on
    
    h = line([0 meanArmx(i)],[i meanArmy(i)+i]-0.5,'Color',0.4.*[1 1 1],'LineWidth',2);

    drawnow
end
ylim([0 6])


axes(h2)

bClr = 0.95.*[1 1 1];

arms = 2:2:6;

% Loop thru every-other arm to create light bars in background
for i = 1:length(arms)
    
    % Current arm number
    n = arms(i);
    
    % Coords for bar
    xVals = [xL(1) xL(2) xL(2) xL(1) xL(1)];
    yVals = [n n n+hBar n+hBar n];

    % Make bar
    h = fill(xVals,yVals,bClr,'EdgeColor','none');
    hold on
end


% Loop thru each footprint
for i = 1:size(xTip,2)
    
    % Index for values
    idx = ~isnan(xTip(:,i));
    
    if sum(idx)>0
        % Check
        if min(xTip(idx,i))~=max(xTip(idx,i))
            error('Problem!');
        end
        
        % Current arm number
        n = armNum(i);
        
        % Current foot number
        nFt = ftNum(i)-1;
        
        % Start and end times
        tStart = S.t(find(~isnan(xTip(:,i)),1,'first'));
        tEnd   = S.t(find(~isnan(xTip(:,i)),1,'last'));
        
        % Micro bar height
        mhBar = hBar/maxFt;
        
        % Coords for bar
        xVals = [tStart tEnd tEnd tStart tStart];
        yVals = [n n n+mhBar n+mhBar n];
        
        % Offset by foot number
        yVals = yVals + nFt.*mhBar;
        
        % Plot box
        axes(h2)
        h = fill(xVals,yVals,cmap2(nFt+1,:));
        set(h,'FaceAlpha',fAlpha,'EdgeColor','none')
        hold on
    end
end

ylim([1 7])

% Position of top graph
pos3 = get(ax1(1),'Position');

% Position of gait graph
pos2 = get(h2,'Position');

% Position of body graphs
pos1 = get(h1,'Position');

% Adjust position of gait graph
% set(h2,'Position',[pos2(1)*0.8 pos2(2)*1.3 pos2(3) pos1(4)*.88],...
%     'YColor','none','TickDir','out');
 set(h2,'YColor','none','TickDir','out');

% Adjust position of top graph
set(h1,'XColor','none','YColor','none','Position',[pos1(1)*1.4 pos1(2:4)]);

set(ax1,'XColor','none')

% set(ax1,'Position',[pos2(1)*0.8 pos3(2) pos2(3) pos3(4)])
xlabel('Time (s)')



function feetMap(S,ttext)

maxFt = 6;

% Color of the seastar body
sClr = 0.8.*[1 1 1];

% Radius and angular position of body centers
rBodCntr =  1.1;
%angBodCntr = [30:60:360]./180*pi;
angBodCntr = S.armAngLT;

% Body length/diameter
bLen = mean([range(S.armLT.x) range(S.armLT.y)]);

% Coordinate for seastar body
starPts = makeStar(S.armLT.x./bLen, S.armLT.y./bLen, 0, 0);

f = figure;

% Pull default colormap
cmap = colormap(f,'parula');
cval = linspace(0,1.2,size(cmap,1));

% Custom colormap by interpolation
for i = 1:maxFt
    cmap2(i,1) = interp1(cval,cmap(:,1),i./maxFt);
    cmap2(i,2) = interp1(cval,cmap(:,2),i./maxFt);
    cmap2(i,3) = interp1(cval,cmap(:,3),i./maxFt);
end



for i = 1:6
    
    % Offset of points
    xOff = rBodCntr * cos(angBodCntr(i));
    yOff = rBodCntr * sin(angBodCntr(i));
    
    % Draw body
    h = fill(starPts(:,1)+xOff,starPts(:,2)+yOff,sClr,'EdgeColor','none');
    axis equal
    hold on
    
    for j = 1:length(S.ft)
        
        if i==S.ft(j).armNum
            
            k = S.ft(j).footNum;
            
            % Indices for start and end of contact
            idx    = ~isnan(S.ft(j).xBase);
            iStart = find(~isnan(S.ft(j).xBase),1,'first');
            iEnd  = find(~isnan(S.ft(j).xBase),1,'last');
            
            % Tube foot motion
%             xVals = S.ft(j).xTipLT([iStart iEnd]) ./ bLen + xOff;
%             yVals = S.ft(j).yTipLT([iStart iEnd]) ./ bLen + yOff;
            xVals = S.ft(j).xTipLT(idx) ./ bLen + xOff;
            yVals = S.ft(j).yTipLT(idx) ./ bLen + yOff;
            h = line(xVals,yVals,'Color',cmap2(k,:),'LineWidth',3);
            
            % Start point
            xVals = S.ft(j).xTipLT([iStart]) ./ bLen + xOff;
            yVals = S.ft(j).yTipLT([iStart]) ./ bLen + yOff;
            h = scatter(xVals,yVals,'MarkerFaceColor',cmap2(k,:),...
                        'MarkerEdgeColor','none','SizeData',50);
            
%             % Start line
%             xVals = [S.ft(j).xTipLT(iStart) S.ft(j).xBaseLT(iStart)] ./ bLen + xOff;
%             yVals = [S.ft(j).yTipLT(iStart) S.ft(j).yBaseLT(iStart)] ./ bLen + yOff;
%             h = line(xVals,yVals,'Color','k','LineWidth',1);
%             
%             % End line
%             xVals = [S.ft(j).xTipLT(iEnd) S.ft(j).xBaseLT(iEnd)] ./ bLen + xOff;
%             yVals = [S.ft(j).yTipLT(iEnd) S.ft(j).yBaseLT(iEnd)] ./ bLen + yOff;
%             h = line(xVals,yVals,'Color','k','LineWidth',1);
        end
        
        
       
    end
    
    ttt=3; 
end

set(gca,'XColor','none','YColor','none')
hold off
title(ttext)















