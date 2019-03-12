function anaMaster(orientation,SSnum,seqnum)
% Runs analysis of kinematic data
% Direction       - 'h', 'u',  or 'v'
% SSnum             - Individual number of sea star.
% seqnum            - number of sequenc


%% Parameters

% Camera view (for 2D analysis)
camView = 'canon';


% Initial duration for identifying a trajectory-based coordinate system (s)
initialDur = 15;

% Dock figure windows
set(0,'DefaultFigureWindowStyle','docked')


%% Manage paths

paths = givePaths;


%% List of sequences to analyze

% Vist of all video files
cList0 = catVidfiles(paths.vid,camView);

% Define path relative to root

% Initialize index
k = 1;

% Make cList to only include finished sequences
cList = struct();
for j = 1:length(cList0.path)
    if ~isempty(dir([paths.data filesep cList0.path{j} filesep ...
                     cList0.fName{j} filesep 'Bundled Data.mat']))
        
        % Transfer all data to cList
        cList.age(k,1)      = cList0.age(j);
        cList.indiv(k,1)    = cList0.indiv(j);
        cList.orient(k,1)   = cList0.orient(j);
        cList.vidType{k,1}  = cList0.vidType{j};
        cList.path{k,1}     = cList0.path{j};
        cList.fName{k,1}    = cList0.fName{j};
        cList.ext{k,1}      = cList0.ext{j};
        cList.calPath{k,1}  = cList0.calPath{j};
        
        k = k + 1;
    end
end

if k==1
    error('No matching videos found in list')
end

clear cList0


%% Calculate local coordinates

% Loop thru sequences
for i = 1:length(cList.age)
    
    % Load bundled 2D data ('S')
    load([paths.data filesep cList.path{i} filesep ...
        cList.fName{i} filesep 'Bundled Data.mat']);  
     
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
    
    
    
    % Visualize trajectory transformation
    if 0
       figure;
       
       subplot(1,2,1)
       h = scatter(S.xCntr(~iDur),S.yCntr(~iDur),'MarkerEdgeColor','none',...
              'MarkerFaceColor','k','Sizedata',10);
       hold on
       scatter(S.xCntr(iDur),S.yCntr(iDur),'MarkerEdgeColor','r',...
               'MarkerfaceColor','r','Sizedata',10)
       plot([pStart(1) pEnd(1)],[pStart(2) pEnd(2)],'r-')
       hold off
       axis equal
       
       subplot(1,2,2)
       h = scatter(CntrPntsT(~iDur,1),CntrPntsT(~iDur,2),'MarkerEdgeColor','none',...
              'MarkerFaceColor','k','Sizedata',10);
       hold on
       scatter(CntrPntsT(iDur,1),CntrPntsT(iDur,2),'MarkerEdgeColor','r',...
               'MarkerfaceColor','r','Sizedata',10);
       axis equal
       hold off
    end
    
    
    ttt=3;
     %baseL = transCoord2d('ang G2L',cOrigin,-S.ang(j)/180*pi,baseG);
end


%% Animate 








%% Plot data


figure;

subplot(1,2,1)
plot(S.xCntr,S.yCntr,'-k')
axis equal
hold on

% Loop thru frames
for i = 1:length(S.xCntr)
    
    % Loop trhu feet
    for j = 1:length(S.ft)
        
        % if there's a base point . . .
        if sum(~isnan(S.ft(j).xBase))>0
            
            % If an attachment frame . . .
            if i==S.ft(j).iStart
                
                h = line([S.ft(j).xBase(i) S.ft(j).xTip(i)],[S.ft(j).yBase(i) S.ft(j).yTip(i)]);
                
            % If a release frame . . .
            elseif i==S.ft(j).iEnd
                
                h = line([S.ft(j).xBase(i) S.ft(j).xTip(i)],[S.ft(j).yBase(i) S.ft(j).yTip(i)]);
                
            end
        end
    end
end

ttt=3;





% 
% % Current paths
% currpaths.data = [paths.data filesep cList.path{i} filesep cList.fName{i}];
% currpaths.vid  = [paths.vid filesep cList.path{i} filesep cList.fName{i} cList.ext{i}];
% 
% % Load video info (v)
% v = defineVidObject(currpaths.vid,'MOV');
% 
% % Load data bundle ('S')
% load([currpaths.data filesep 'Bundled Data.mat'])
% 
% % Define frames vector
% %frames = v.UserData.FirstFrame:v.UserData.LastFrame;
% 
% % Load initial conditions (iC)
% load([currpaths.data filesep 'Initial conditions'])
% 
% % % Cue up coordinates
% % xBase = nan(length(H.ft(1).xBase),length(H.ft));
% % yBase = xBase;
% % xTip  = xBase;
% % yTip  = xBase;
% 
% % Create matrices of 
% for i = 1:(length(H.ft)-2)
%     xBase(:,i)  = H.ft(i).xBase;
%     yBase(:,i)  = H.ft(i).yBase;
%     xTip(:,i)   = H.ft(i).xTip;
%     yTip(:,i)   = H.ft(i).yTip;
% end
% 
% % Data frames
% tmp = ~isnan(xTip);
% dFrames = H.frames(sum(tmp,2)>0);
% clear tmp
% 
% 
% % Loop trhu frames
% for i = min(dFrames):max(dFrames)
%     
%     cFrame = H.frames(i);
%     
%     %  frame
%     im = getFrame(currpaths.vid,v,cFrame,imInvert,'gray');
%     
%     
%     imshow(im,'InitialMag','fit')
%     hold on
%     
%     for j = 1:length(H.ft)
%         plot(H.ft(j).xTip(cFrame),H.ft(j).yTip(cFrame),'or')
%     end
%     
%     title(['Frame ' num2str(cFrame)])
%     
%     hold off
%     
%     pause(0.1)
% end
% ttt= 3;



function ptsL = G2L(cOrigin,theta,ptsG)
% Coordinate transformation from global to local coordinates
% cOrgin    - Coordinate origin (3x1)
% theta     - Azimuth angle of local system wrt global
% ptsG      - Coordinates in the global FOR (3xn)

if (size(cOrigin,2) > size(cOrigin,1)) || ...
   ((length(ptsG)>3) && (size(ptsG,1) > size(ptsG,2)))  
    error('All points should be arrange in column vectors');
end

% Output nans, if nans are input
if sum(~isnan(ptsG))==0
    ptsL = ptsG;
    
% Otherwise . . .
else
    
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
end

function ptsG = L2G(cOrigin,theta,ptsL)
% Coordinate transformation from global to local coordinates
% cOrgin    - Coordinate origin (3x1)
% theta     - Azimuth angle of local system wrt global
% ptsG      - Coordinates in the global FOR (3xn)

if (size(cOrigin,2) > size(cOrigin,1)) || ...
   ((length(ptsG)>3) && (size(ptsG,1) > size(ptsG,2)))  
    error('All points should be arrange in column vectors');
end

% Output nans, if nans are input
if sum(~isnan(ptsG))==0
    ptsL = ptsG;
    
% Otherwise . . .
else
    
    % Remove nans
    idx = ~isnan(ptsG(1,:));
    ptsL = ptsL(:,idx);
    
    % Add z-dimension
    cOrigin = [cOrigin; 0];
    ptsL    = [ptsL; zeros(1,size(ptsL,2))];
    
    % Axes defined
    xaxis = [cos(theta) sin(theta) 0]';
    yaxis = [cos(theta+pi/2) sin(theta+pi/2) 0]';
    zaxis = [0 0 1]';
    
    % Rotation matrix
    S = [xaxis yaxis zaxis];
    
    % Local coordinate
    ptsG = local2globalcoord(ptsL,'rr',cOrigin,S);
    
    % Remove z-dimension
    ptsG = ptsG(1:2,:);
end
