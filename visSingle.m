function visSingle(Orientation,SSnum,seqnum)
% Visualize the results of a single sequence
% Orientation - indicates the wall orientation ('h','v','u')
% SSnum       - number of seastar
% seqnum      - sequence number


%% Deal with inputs

if nargin >=3
    
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
    
    clear Orientation SSnum seqnum
    
else
    error('You need to request the sequence')
    
end


%% Parameters

% Paths
paths = givePaths;

%Camera view (for 2D analysis)
camView = 'canon';

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


%% Transfor in a trajectory coordinate system

% Load bundled 2D data ('S')
load([paths.data filesep good.path filesep ...
      good.fName filesep 'Bundled Data.mat']);

% Index of early points
iDur = S.t<=initialDur;

% Linear fit to trajectory
cX = polyfit(S.t(iDur),S.xCntr(iDur),1);
cY = polyfit(S.t(iDur),S.yCntr(iDur),1);

% Starting and end points
pStart = [polyval(cX,min(S.t(iDur))) polyval(cY,min(S.t(iDur)))];
pEnd   = [polyval(cX,max(S.t(iDur))) polyval(cY,max(S.t(iDur)))];

% Cntr points in global FOR
CntrPnts = [S.xCntr' S.yCntr'];

% Cntr points in tarjectory FOR
CntrPntsT = transCoord2d('xax G2L',pStart,pEnd,CntrPnts);



% Diameter of sea star
aLen = max(hypot(S.armL.x,S.armL.y));

% Limits of x-axis
xL = 1.4*[min(CntrPntsT(:,1))-aLen max(CntrPntsT(:,1))+aLen];


%% Visualize

if 1
    figure
    
    
    % Step thru time
    for i = 1:length(S.t)
        
        h = scatter(CntrPntsT(:,1),CntrPntsT(:,2),'MarkerEdgeColor','none',...
            'MarkerFaceColor',clr{3},'Sizedata',10);
        a1 = gca;
        hold on
        axis equal
        
        % Arm points in trajectory FOR
        armPntsT = transCoord2d('xax G2L',pStart,pEnd,[S.arm(i).x, S.arm(i).y]);
        
        % Points for whole body
        starPtsT = makeStar(armPntsT(:,1), armPntsT(:,2), CntrPntsT(i,1), ...
                            CntrPntsT(i,2),aLen);
        
        % Cntr points in tarjectory FOR
       % CntrPntsT = transCoord2d('xax G2L',pStart,pEnd,CntrPnts);
        
        h = fill(starPtsT(:,1),starPtsT(:,2),clr{1},'FaceAlpha',0.3,...
            'EdgeColor','none');
        
        % Center point
        scatter(CntrPntsT(i,1),CntrPntsT(i,2),'MarkerEdgeColor','r',...
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
        
        
        
        
        
    end
    
end


function starPts = makeStar(xA,yA,xC,yC,aLen)
% Create periperal shape of seastar body

nArms = length(xA);

% xA = [xA;xA(1)];
% yA = [yA;yA(1)];

nPts = 500;

% Diameter of cental disc
D = aLen*0.75;

% Radial position of arms, sorted
[Atheta,idx] = sort(wrapTo2Pi(atan2(yA-yC,xA-xC)));

% Arm lengths
ALen = hypot(xA(idx)-xC,yA(idx)-yC);

% Tack on starting arm
Atheta = unwrap([Atheta; Atheta(1)]);
ALen   = [ALen; ALen(1)];

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
      
    rVal0 = (ALen(i)-D).*(cos(tVal0)+1)/2+D;
    rVal1 = (ALen(i+1)-D).*(cos(tVal1)+1)/2+D;
    
    [x0,y0] = pol2cart(theta0,rVal0);
    [x1,y1] = pol2cart(theta1,rVal1);
    
    xP = [xP; x0; x1];
    yP = [yP; y0; y1];
    
    if 0
       subplot(2,2,1)
       plot(tVal0,rVal0,'-',tVal1,rVal1,'-')
       hold on
       
       subplot(2,2,3)
       plot(tVal0,theta0,'-',tVal1,theta1,'-')
       hold on
       
       subplot(2,2,[2 4])
       plot(xP,yP,'-')
       axis equal
       hold on
       
       ttt=3;
    end
end

starPts = [xP+xC yP+yC];








