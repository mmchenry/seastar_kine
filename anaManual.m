function F = anaManual(seq,S,H,aud,Body,iC)
% Analyzes manually collected data
% Running without inputs ('anaManual' at the command line) executes batchMode, 
% which attempts to analyze all
% sequences where ana_video==1 from Weights_experiments.xlsx.
%
% The analysis of individual sequences draws from the audio syncing data,
% bottom manual data, initial conditions, and body data (from automated
% tracking) to do a 3D reconstrcution of the kinematics of individual tube
% feet

% TODO: Finish 3D reconstruction code


%% Define mode

if nargin < 1
    batchMode = 1;
else
    % Check inputs
    if nargin<4
        error('You need to provide values for seq, S, H, and aud')
    end
    batchMode = 0;
end

% Get root paths
paths = givePaths;

set(0,'DefaultFigureWindowStyle','docked')

%% Execution control

% Animate the manually-selected bottom-view data on video frames
do.visBottom = 0;


%% Read data from catalog spreadsheet

if batchMode

% Extension for video file names
extVid = 'MOV';

% Get list of all side csv files
%     dlc_path = [paths.data filesep 'data' filesep 'side_view' filesep '*.csv'];
dlc_path = [paths.data filesep 'rawCSV' filesep '*.csv'];
dlc_files = dir(dlc_path);

% Check for dlc files
if isempty(dlc_files)
    error(['No DLC csv files found in ' dlc_path]);
end

% Throw error if spreadsheet not in place
if exist([paths.data filesep 'Weights_experiments.xlsx'],'file')==0
    error(['Missing spreadsheet: ' paths.data filesep ...
        'Weights_experiments.xlsx'])
end

% Read table
T = readtable([paths.data filesep 'Weights_experiments.xlsx']);

j = 1;seq = [];

% Step thru rows
for i = 1:length(T.date)

    % If using video at all and currently analyzing
    %     if T.use_video(i)==1 && T.ana_video(i)==1 && T.complete_video(i)==0
    if T.use_video(i)==1 && T.ana_video(i)==1

        seq(j).dateNum        = datenum(T.date(i));
        seq(j).ext            = extVid;
        seq(j).fName_side     = [T.side_filename{i}];
        seq(j).fName_bot      = [T.bottom_filename{i}];
        seq(j).fName_calSide  = [T.side_cal_filename{i}];
        seq(j).fName_calBot   = [T.bot_cal_filename{i}];
        seq(j).expType        = T.exp_type{i};
        seq(j).indiv          = T.indiv_num(i);
        seq(j).addMass        = T.added_mass(i);
        seq(j).floatNum       = T.float_num(i);
        seq(j).bodyMass       = T.body_mass(i);
        seq(j).calConst       = T.cal_side(i);
%         seq(j).calConst_bot   = T.cal_bottom(i); % Note: do not use this cal constant
        seq(j).fps_side       = T.frame_rate_side(i);
        seq(j).fps_bot        = T.frame_rate_bottom(i);
        seq(j).SW_percent     = T.percent_sw(i);

        % Check experiment type
        if strcmp(seq(j).expType,'c')
            seq(j).dirName = 'control';
        elseif strcmp(seq(j).expType,'w')
            seq(j).dirName = 'weights';
        elseif strcmp(seq(j).expType,'f')
            seq(j).dirName = 'floats';
        else
            error(['Do not recognize experiment type: ' seq(i).expType])
        end

        j = j + 1;
    end
end

% Check for seq
if isempty(seq)
    error('No sequence data found from spreadsheet');
end

clear j T

end % Batch mode


%% Batch mode

if batchMode

% Load side kinematics data (S)
% load([paths.data filesep 'SideDataPooled.mat'])
load([paths.data filesep 'SideDataPooled_eventsManual2.mat']);

disp('Batch mode --------------')

% Loop thru sequences
for i = 1:length(seq)
    
    % Load audio sync data (aud)
    load([paths.data filesep seq(i).dirName filesep 'bottom' filesep ...
        'matlabData2021' filesep seq(i).fName_bot filesep 'audio_delay.mat']);

    % Load bottom manual data (H)
    load([paths.data filesep seq(i).dirName filesep 'bottom' filesep ...
        'matlabData2021' filesep seq(i).fName_bot filesep 'ManualFootData.mat'])
    
    % Load body data (Body)
    load([paths.data filesep seq(i).dirName filesep 'bottom' filesep ...
        'matlabData2021' filesep seq(i).fName_bot filesep 'Body.mat'])

    % Load frame rates (frRate)
    % load([paths.data filesep seq(i).dirName filesep 'bottom' filesep ...
    %       'matlabData2021' filesep seq(i).fName_bot filesep 'frame_rate.mat'])

    % Load initial conditions (iC)
     load([paths.data filesep seq(i).dirName filesep 'bottom' filesep ...
           'matlabData2021' filesep seq(i).fName_bot filesep 'Initial conditions.mat'])

    % Find matching data from compiled side data
    iMatch = nan;
    for j = 1:length(S)

        if strcmp(S(j).fName_bot,seq(i).fName_bot)
            iMatch = j;
            break
        end
    end

    % Check index
    if isnan(iMatch), error('No match in the compiled side view data'); end

    % Run present function
    F(i) = anaManual(seq(i),S(iMatch),H,aud,Body,iC);

    disp(['    Done ' num2str(i) ' of ' num2str(length(seq))])
end

% Save data 
save([paths.data filesep 'footData.mat'],'F')

F = [];

end % batchMode

%% Visualize bottom data

if do.visBottom
% Translate tip and base points into local and global FOR
H = addLocalFeet(H,Body.xCntr,Body.yCntr,Body.Rotation.rot_ang);

% Visualize feet onto video frames
visBotCoords(paths,seq,H,Body,iC)

end %do.visBottom


%% Calculate podium length and angle

% Run for single sequence
if ~batchMode

visCheck = 0;

% Number of points to trim from edges for range calibration
eTrim = 30;

% Time vectors, taking audio sync and kinematics cross-correlation into account
t_B = (H.frames'-min(H.frames)) ./ seq.fps_bot;
t_S = S.t-min(S.t) + aud.delay;

% Check for matching frame vectors
if length(Body.frames) ~= length(H.frames) || ...
          Body.frames(1)~=H.frames(1) || ...
          Body.frames(end)~=H.frames(end)       
    error('Frame vectors do not match up between body and tube feet')
end

% Adjust scaling of coordinates 
tMin        = max([min(t_B(2:end))+eTrim/seq.fps_bot ...
                   min(t_S)+eTrim/seq.fps_side]);
tMax        = min([max(t_B)-eTrim/seq.fps_bot ...
                   max(t_S)-eTrim/seq.fps_side]);
iB          = t_B>tMin & t_B<tMax;
iS          = t_S>tMin & t_S<tMax;
x_B         = Body.xCntr;
x_S         = S.xRaw;

% Range calibration 
rangeCal = range(x_S(iS)) / range(x_B(iB));

% Visual check on x-axis from two dimensions
if 0
    % Zero wrt minumum values
    x_B = smooth(x_B(iB)-min(x_B(iB)),50).*rangeCal;
    x_S = smooth(x_S(iS)-min(x_S(iS)),50);

    % Zero wrt minumum values
    s_B = diff(x_B)./diff(t_B(t_B>tMin & t_B<tMax));
    s_S = diff(x_S)./diff(t_S(t_S>tMin & t_S<tMax));
    ts_B    = t_B(iB);
    ts_S    = t_S(iS);
    ts_B    = ts_B(2:end);
    ts_S    = ts_S(2:end);

    figure
    subplot(2,1,1)
    plot(t_B(t_B>tMin & t_B<tMax),x_B,'-',...
         t_S(t_S>tMin & t_S<tMax),x_S,'-')
    title('raw')
    grid on
    subplot(2,1,2)
    plot(ts_B,s_B,'-',ts_S,s_S,'-')
    ylabel('Speeds')
    grid on
end

% Restate the timing of landmarks
tLand = S.tLand-min(S.t) + aud.delay;

% Side view coordinates
xSide  = S.xRaw ./ S.calConst;
zSide  = S.yRaw ./ S.calConst;

% Translate tip and base points into local and global FOR
H = addLocalFeet(H,Body.xCntr,Body.yCntr,Body.Rotation.rot_ang);

% Make figure windows
if visCheck
    f1 = figure;
    f2 = figure;
    f3 = figure;
end

% Extract the range of indicies for the power stroke among all tube feet
iStart = inf; iEnd = 0;
for i = 1:length(H.ft)  
    if H.ft(i).iStart<iStart
        iStart = H.ft(i).iStart;
    end
    if H.ft(i).iEnd>iEnd
        iEnd = H.ft(i).iEnd;
    end
end

% Index for side-view kinematics that span the power strokes
iS_s = t_S >= tLand(find(tLand<t_B(iStart),1,'last')) & ...
       t_S <= tLand(find(tLand>t_B(iEnd),1,'first'));

% Restate z-coordinates relative to current power stroke
zSide = zSide - nanmin(zSide(iS_s));

% Store in F
F            = H;
F.rangeCal   = rangeCal;
F.calConst_s = S.calConst;
F.seq        = seq;

% Loop thru feet, each with one pwr stroke
for i = 1:length(H.ft)

    % Flag if there are multiple power strokes
    iPwr = find(~isnan(H.ft(i).xTip));
    if max(diff(find(iPwr)))>1
        error('multiple power strokes for this foot')
    end

    % Index for power stroke in bottom view
    iB = H.ft(i).iStart:H.ft(i).iEnd;

    % Index for power stroke in side view
    iS = t_S>=t_B(H.ft(i).iStart) & t_S<=t_B(H.ft(i).iEnd);

    % Time values for current pwr stroke
    tPwr = t_B(iB);

    % Current foot coordinates
    xBase = H.ft(i).xBase(iB) .* rangeCal/S.calConst;
    yBase = H.ft(i).yBase(iB) .* rangeCal/S.calConst;
    xTip  = H.ft(i).xTip(iB)  .* rangeCal/S.calConst;
    yTip  = H.ft(i).yTip(iB)  .* rangeCal/S.calConst;

    % Interpolate wrt time to get z-coordinates for foot
    zTip  = 0.*yTip; 
    zBase = interp1(t_S,zSide,tPwr);

    % Restate coordinates, using tip coordinates as origin
    origin = [mean(xTip) mean(yTip)];

    % Transform coords relative to the plane of the power stroke
    c      = polyfit(xBase,yBase,1);
    xL     = max([range(xBase) range(yBase)])/2;
    xaxis  = [origin(1)+xL c(1)*(origin(1)+xL)+(origin(2)-c(1)*origin(1))];
    C      = defineSystem2d('x-axis',origin,xaxis);
    coordL = transCoord2d('G2L',C.tform,[xBase yBase]);
    xBaseL = coordL(:,1);
    yBaseL = coordL(:,2);
    
     % Podium length
    pod_len = hypot(xBaseL,yBaseL);

    % Podium angle
    theta = atan2(zBase,xBaseL);

    if visCheck

        figure(f1);

        subplot(3,2,1)
        plot(xBase,yBase,'ko',xBase,polyval(c,xBase),'r-',...
             origin(1),origin(2),'r+',xaxis(1),xaxis(2),'go')
        xlabel('xG');ylabel('yG')
        axis equal

        subplot(3,2,2)
        plot3(xBaseL,yBaseL,zBase,'ko')
        hold on
        for j = 1:length(xBaseL)
            plot3([0 xBaseL(j)],[0 yBaseL(j)],[0 zBase(j)],'k-')
        end
        axis equal
        xlabel('X_L');ylabel('Y_L');zlabel('Z_L')
        view([0 0])
        hold off

        subplot(3,2,[3:4])
        plot(tPwr,pod_len,'-')
        xlabel('t');ylabel('pod len');

        subplot(3,2,[5:6])
        plot(tPwr,theta.*180/pi,'-')
        xlabel('t');ylabel('\theta');

        figure(f2)
        subplot(2,1,1)
        plot3(xBase,yBase,zBase,'ko')
        hold on
        for j = 1:length(xBaseL)
            plot3([origin(1) xBase(j)],[origin(2) yBase(j)],[0 zBase(j)],'k-')
        end
        axis equal
        xlabel('X_L');ylabel('Y_L');zlabel('Z_L')
        title('raw-ish')
        grid on
        hold off

        subplot(2,1,2)
        plot3(xBaseL,yBaseL,zBase,'ko')
        hold on
        for j = 1:length(xBaseL)
            plot3([0 xBaseL(j)],[0 yBaseL(j)],[0 zBase(j)],'k-')
        end
        axis equal
        xlabel('X_L');ylabel('Y_L');zlabel('Z_L')
        title('transformed')
        grid on
        hold off

        figure(f3)
        subplot(4,1,1)
        plot(t_S,xSide,'k',t_S(iS),xSide(iS),'r-',...
             t_S(find(iS,1,'first')),xSide(find(iS,1,'first')),'r+',...
             t_S(find(iS,1,'last')),xSide(find(iS,1,'last')),'r+')
        xlabel('t');ylabel('X')

        subplot(4,1,2)
        plot(tPwr,xBase-min(xBase),'b',t_S(iS),xSide(iS)-min(xSide(iS)),'r')
        xlabel('t');ylabel('X')
        legend('manual','dlc')

        subplot(4,1,3)
        plot(t_S,zSide,'k',t_S(iS),zSide(iS),'r-',...
             t_S(find(iS,1,'first')),zSide(find(iS,1,'first')),'r+',...
             t_S(find(iS,1,'last')),zSide(find(iS,1,'last')),'r+')
        xlabel('t');ylabel('Z')

        subplot(4,1,4)
        plot(tPwr,zBase,'b',t_S(iS),zSide(iS),'r')
        xlabel('t');ylabel('Z')
    end

    clear coordL tform xaxis c C

    % Store results
    F.ftL(i).xBase      = xBaseL;
    F.ftL(i).yBase      = yBaseL;
    F.ftL(i).zBase      = zBase;
    F.ftL(i).theta      = theta;
    F.ftL(i).pod_len    = pod_len;
end

% Save new version of 'H'
% disp(' ')
% disp(['Saving tube foot length and angle data for ' seq.fName_bot])
% save([paths.data filesep seq.dirName filesep 'bottom' filesep ...
%         'matlabData2021' filesep seq.fName_bot filesep 'ManualFootData_len_ang.mat'],'H')

end %if ~batchMode


end %anaManual


function visBotCoords(paths,seq,H,Body,iC)
% Visual check on position of bottom-view coordinates

% video path
vid_path = [paths.vid filesep seq.dirName filesep 'bottom' filesep ...
                    seq.fName_bot '.' seq.ext];

imInvert = 0;

% Video object
v = VideoReader(vid_path);

% Frames
%     v.UserData.FirstFrame = clipInfo.startFrame;
%     v.UserData.LastFrame = clipInfo.endFrame;
%     frames = v.UserData.FirstFrame:v.UserData.LastFrame;

% Get starting and ending indicies
iStart = nan;iEnd = nan;
for i = 1:length(H.frames)
    for j = 1:length(H.ft)
        if ~isnan(H.ft(j).xBase(i))
            if isnan(iStart)
                iStart = H.frames(i);
            else
                iEnd = H.frames(i);
            end
        end
    end
end

% Check that indicies are okay
if isnan(iStart) || isnan(iEnd), error(' ');end

% Cue up the video
v.currentTime = H.frames(iStart)./seq.fps_bot;

f = figure;

for i = iStart:iEnd
    
    % Current image
%     im = getFrame(vid_path,v,H.frames(i),imInvert,'color',[]);


    im = readFrame(v);
     imshow(im,'InitialMag','fit')
    hold on

    plot(Body.xCntr(i)-iC.r,Body.yCntr(i)-iC.r,'w+')

    for j = 1:length(H.ft)
        h(j)  = plot(H.ft(j).xTip(i),H.ft(j).yTip(i),'og');
        h2(j) = plot([H.ft(j).xBase(i) H.ft(j).xTip(i)],...
            [H.ft(j).yBase(i) H.ft(j).yTip(i)],'g');
    end
    title(['Frame = ' num2str(H.frames(i))])
    pause(0.01)
    delete(h);delete(h2)
end
hold off

end %visBotCoords


function H = addLocalFeet(H,xCntr,yCntr,ang)
% Translate tip and base points into local FOR

% Loop thru tube feet
for k = 1:length(H.ft)

    % Fill out fields with nans
    H.ft(k).xBaseL = nan(length(H.frames),1);
    H.ft(k).yBaseL = nan(length(H.frames),1);
    H.ft(k).xTipL  = nan(length(H.frames),1);
    H.ft(k).yTipL  = nan(length(H.frames),1);

    % Index of non-nan frames
    iFrame = find(~isnan(H.ft(k).xBase));

    % Loop trhu frames for the base coordinates
    for j = 1:length(iFrame)

        % Current origin
        cOrigin = [xCntr(iFrame(j)) yCntr(iFrame(j))];

        % Current base point
        baseG = [H.ft(k).xBase(iFrame(j)) H.ft(k).yBase(iFrame(j))];

        %             % Base points in trajectory FOR
        %             baseT = transCoord2d('xax G2L',S.pStart,S.pEnd,baseG);

        % Local base point
        baseL = transCoord2d('ang G2L',cOrigin,-ang(iFrame(j))/180*pi,baseG);

        % Local tip point wrt traj
        %             baseLT = transCoord2d('ang G2L',[0 0],meanAng,baseL);

        % Store local coordinates for base
        H.ft(k).xBaseL(iFrame(j),1)   = baseL(1);
        H.ft(k).yBaseL(iFrame(j),1)   = baseL(2);
        %             S.ft(k).xBaseT(j,1)   = baseT(1);
        %             S.ft(k).yBaseT(j,1)   = baseT(2);
        %             S.ft(k).xBaseLT(j,1)  = baseLT(1);
        %             S.ft(k).yBaseLT(j,1)  = baseLT(2);

    end % j loop

    % Index of non-nan frames
    iFrame = find(~isnan(H.ft(k).xTip));

    % Loop trhu frames for the tip coordinates
    for j = 1:length(iFrame)

        %
        %         % Transfer data
        %         S.ft(k).xTip(j,1)    = H.ft(k).xTip(iH);
        %         S.ft(k).yTip(j,1)    = H.ft(k).yTip(iH);


        % Current origin
        cOrigin = [xCntr(iFrame(j)) yCntr(iFrame(j))];

        % Current tip point
        tipG = [H.ft(k).xTip(iFrame(j)) H.ft(k).yTip(iFrame(j))];

        % Local tip point
        %tipL = G2L(cOrigin,-S.ang(iFrame(j))/180*pi,tipG);
        tipL = transCoord2d('ang G2L',cOrigin,-ang(iFrame(j))/180*pi,tipG);

        %             % Tip points in trajectory FOR
        %             tipT = transCoord2d('xax G2L',S.pStart,S.pEnd,tipG);
        %
        %             % Local tip point wrt traj
        %             tipLT = transCoord2d('ang G2L',[0 0],meanAng,tipL);


        % Store local coordinates
        H.ft(k).xTipL(iFrame(j),1)    = tipL(1);
        H.ft(k).yTipL(iFrame(j),1)    = tipL(2);
        %             S.ft(k).xTipT(iFrame(j),1)    = tipT(1);
        %             S.ft(k).yTipT(iFrame(j),1)    = tipT(2);
        %             S.ft(k).xTipLT(iFrame(j),1)   = tipLT(1);
        %             S.ft(k).yTipLT(iFrame(j),1)   = tipLT(2);

        clear tipL tipG tipLT

    end % j loop
end

% FILL OUT BASE COORDINATES FROM LOCAL COORDINATES --------------------

% Loop thru tube feet
for j = 1:length(H.ft)

    % Loop trhu frames
%     for j = 1:length(H.frames)

        % Find indicies for contact and release
        H.ft(j).iStart = find(~isnan(H.ft(j).xTip),1,'first');
        H.ft(j).iEnd   = find(~isnan(H.ft(j).xTip),1,'last');

        % If no foot number . . .
        if ~isfield(H,'ft') || sum(~isnan(H.ft(j).footNum))==0
            warning(['No foot number given for foot ' num2str(j) ', assuming 2'])
            H.ft(j).footNum      = 2;
            H.ft(j).footLet      = 'A';
        else
            H.ft(j).footNum      = H.ft(j).footNum;
            H.ft(j).footLet      = H.ft(j).footLet;
        end

        % If there's a base point
        if sum(~isnan(H.ft(j).xBase))>0

            % Find local position of base
            idx     = find(~isnan(H.ft(j).xBase),1,'first');
            baseL   = [H.ft(j).xBaseL(idx) H.ft(j).yBaseL(idx)];
%             baseLT  = [H.ft(j).xBaseLT(idx) H.ft(j).yBaseLT(idx)];

            % Step thru frames to transform into global FOR
            for k = 1:length(xCntr)

                % Current origin
                cOrigin = [xCntr(k) yCntr(k)];

                % If there is a tip point
                if ~isnan(H.ft(j).xTip(k))

                    % Calcuate the global position of the base
                    %baseG = G2L(cOrigin,-H.ang(k)/180*pi,baseL);
                    baseG = transCoord2d('ang L2G',cOrigin,-ang(k)/180*pi,baseL);

                    % Overwrite the global position
                    H.ft(j).xBase(k,1)   = baseG(1);
                    H.ft(j).yBase(k,1)   = baseG(2);
                    H.ft(j).xBaseL(k,1)  = baseL(1);
                    H.ft(j).yBaseL(k,1)  = baseL(2);
                    %                 H.ft(j).xBaseLT(k,1) = baseLT(1);
                    %                 H.ft(j).yBaseLT(k,1) = baseLT(2);

                    clear baseG base L base LT
                end
            end

        % If no base point . . .
        elseif sum(~isnan(H.ft(j).xTip))>0
            warning(['    Foot ' num2str(j) ' in does not have a base, ' ...
                'but does have tip coordinates'])
        end


        clear iH
    end

% end

end %addLocalFeet








