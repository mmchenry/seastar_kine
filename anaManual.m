function F = anaManual(seq,S,H,aud,Body,iC)
% Analyzes manually collected data


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
load([paths.data filesep 'SideDataPooled.mat'])

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
end

end % batchMode


%% Load data

% Run for single sequence
if ~batchMode

% Time vectors, taking audio sync into account
t_B = (H.frames-min(H.frames)) ./ seq.fps_bot;
t_S = S.t-min(S.t) + aud.delay;

% Side view coordinates
xSide  = S.xRaw;
zSide  = S.yRaw;

%TODO: coordinate transform side coordinates 

% Check for matching frame vectors
if length(Body.frames) ~= length(H.frames) || ...
          Body.frames(1)~=H.frames(1) || ...
          Body.frames(end)~=H.frames(end)       
    error('Frame vectors do not match up')
end

H = addLocalFeet(H,Body.xCntr,Body.yCntr,Body.Rotation.rot_ang);

visBotCoords(paths,seq,H,Body,iC)



ttt=3;

return


% Loop thru feet
for i = 1:length(H.ft)
    
    % Current foot coordinates
    xBase = H.ft(i).xBase;
    yBase = H.ft(i).yBase;
    xTip  = H.ft(i).xTip;
    yTip  = H.ft(i).yTip;

    % Index for power stroke
    iPwr = ~isnan(xTip);

    % Identify multiple power strokes
    if max(diff(find(iPwr)))>1
        error('multiple power strokes for this foot')
    end

    ttt = 3;

    %TODO: Use linear interpolation to match the foot and side coordinate data 
    
    % Find body-coordinate transformation code to find base coordinates for
    % the whole power stroke

end

ttt = 3;

end %for


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

    % Loop trhu frames
    for j = 1:length(H.frames)

        %         % Index for current frame for manual data
        %         iH = find(H.frames==frames(j),1,'first');
        %
        %         % Transfer data
        %         S.ft(k).xBase(j,1)   = H.ft(k).xBase(iH);
        %         S.ft(k).yBase(j,1)   = H.ft(k).yBase(iH);
        %
        % If base point is a nan . . .
        if isnan(H.ft(k).xBase(j))

            H.ft(k).xBaseL(j,1)  = nan;
            H.ft(k).yBaseL(j,1)  = nan;
            %             S.ft(k).xBaseLT(j,1) = nan;
            %             S.ft(k).yBaseLT(j,1) = nan;
            %
            % If there is a tip point
        else
            %
            % Current origin
            cOrigin = [xCntr(j) yCntr(j)];

            % Current base point
            baseG = [H.ft(k).xBase(j) H.ft(k).yBase(j)];

            %             % Base points in trajectory FOR
            %             baseT = transCoord2d('xax G2L',S.pStart,S.pEnd,baseG);

            % Local base point
            baseL = transCoord2d('ang G2L',cOrigin,-ang(j)/180*pi,baseG);

            % Local tip point wrt traj
            %             baseLT = transCoord2d('ang G2L',[0 0],meanAng,baseL);

            % Store local coordinates
            H.ft(k).xBaseL(j,1)   = baseL(1);
            H.ft(k).yBaseL(j,1)   = baseL(2);
            %             S.ft(k).xBaseT(j,1)   = baseT(1);
            %             S.ft(k).yBaseT(j,1)   = baseT(2);
            %             S.ft(k).xBaseLT(j,1)  = baseLT(1);
            %             S.ft(k).yBaseLT(j,1)  = baseLT(2);
            %
            %             clear baseL baseLT
        end
        %
        %         % Transfer data
        %         S.ft(k).xTip(j,1)    = H.ft(k).xTip(iH);
        %         S.ft(k).yTip(j,1)    = H.ft(k).yTip(iH);

        % If tip point is a nan . . .
        if isnan(H.ft(k).xTip(j))
            %
            S.ft(k).xTipL(j,1)   = nan;
            S.ft(k).yTipL(j,1)   = nan;
            % %             S.ft(k).xTipLT(j,1)  = nan;
            % %             S.ft(k).yTipLT(j,1)  = nan;
            %
            % If there is a tip point
        elseif ~isnan(H.ft(k).xTip(j))

            % Current origin
            cOrigin = [xCntr(j) yCntr(j)];

            % Current tip point
            tipG = [H.ft(k).xTip(j) H.ft(k).yTip(j)];

            % Local tip point
            %tipL = G2L(cOrigin,-S.ang(j)/180*pi,tipG);
            tipL = transCoord2d('ang G2L',cOrigin,-ang(j)/180*pi,tipG);

            %             % Tip points in trajectory FOR
            %             tipT = transCoord2d('xax G2L',S.pStart,S.pEnd,tipG);
            %
            %             % Local tip point wrt traj
            %             tipLT = transCoord2d('ang G2L',[0 0],meanAng,tipL);


            % Store local coordinates
            H.ft(k).xTipL(j,1)    = tipL(1);
            H.ft(k).yTipL(j,1)    = tipL(2);
            %             S.ft(k).xTipT(j,1)    = tipT(1);
            %             S.ft(k).yTipT(j,1)    = tipT(2);
            %             S.ft(k).xTipLT(j,1)   = tipLT(1);
            %             S.ft(k).yTipLT(j,1)   = tipLT(2);

            clear tipL tipG tipLT
        end

    end
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
                    baseG = transCoord2d('ang L2G',cOrigin,-ang(j)/180*pi,baseL);

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








