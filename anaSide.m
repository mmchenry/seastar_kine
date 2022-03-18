function anaSide
% Runs analysis of the kinematics of the chip side view 

%% Execution control

% Import DLC data again
do.reImportData = 0;

% Analyze audio data to sync the timing of videos
do.anaAudioSync = 1;

% Annotate dataset manually
do.manAnnotate = 0;

% Visualize events for all sequences
do.visEvents = 0;

% Visualize details of each sequence
do.visSeqs = 0;

% Plot summary boxplot data
do.summaryPlots = 0;

% Plot data that keeps track of individuals
do.indivPlots = 0;

% Plot wrt trial number
do.trialPlot = 0;

% Run statistics
do.stats = 1;

% Marker and line colors
mClr = 0.5.*[1 1 1];
lClr = [0 0 0];


%% Parameters

% Get root paths
paths = givePaths;

% Extension for video file names
extVid = 'MOV';

% Get list of all side csv files
dlc_path = [paths.data filesep 'data' filesep 'side_view' filesep '*.csv'];
dlc_files = dir(dlc_path);

% Check
%if isempty(dlc_files)
    %error(['No DLC csv files found in ' dlc_path]);
%end


%% Read data from catalog spreadsheet

if do.reImportData || do.anaAudioSync
    
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
            
%             seq(j).dateNum        = datenum(T.date(i));
%             seq(j).ext            = extVid;
%             seq(j).fName_side     = [T.side_filename{i}];
%             seq(j).fName_bot      = [T.bottom_filename{i}];
%             seq(j).fName_calSide  = [T.side_cal_filename{i}];
%             seq(j).fName_calBot   = [T.bot_cal_filename{i}];
%             seq(j).expType        = T.exp_type{i};
%             seq(j).indiv          = T.indiv_num(i);
%             seq(j).addMass        = T.added_mass(i);
%             seq(j).floatNum       = T.float_num(i);
%             seq(j).bodyMass       = T.body_mass(i);
%             seq(j).calConst       = T.cal_side(i);
%             seq(j).fps            = T.frame_rate_side(i); 
%             seq(j).SW_percent     = T.percent_sw(i);
%             seq(j).trial          = T.trial_number(i);
%             seq(j).propBounce     = T.prop_bounce(i);
%             seq(j).SW_tot         = T.tot_sw(i);

            seq(j).dateNum         = datenum(T.date(i));
            seq(j).ext             = extVid;
            seq(j).fName_side      = [T.side_filename{i}];
            seq(j).fName_bot       = [T.bottom_filename{i}];
            seq(j).fName_calSide   = [T.side_cal_filename{i}];
            seq(j).fName_calBot    = [T.bot_cal_filename{i}];
            seq(j).expType         = T.exp_type{i};
            seq(j).indiv           = T.indiv_num(i);
            seq(j).addMass         = T.added_mass(i);
            seq(j).floatNum        = T.float_num(i);
            seq(j).bodyMass        = T.body_mass(i);
            seq(j).calConst        = T.cal_side(i);
            seq(j).fps             = T.frame_rate_side(i); 
            seq(j).SW_percent      = T.percent_sw(i);
            seq(j).trial           = T.trial_number(i);
            seq(j).numFtCon1       = T.prop_bounce1(i);
            seq(j).numFtCon2       = T.prop_bounce2(i);
            seq(j).numFtCon3       = T.prop_bounce3(i);
            seq(j).numFtCon4       = T.prop_bounce4(i);
            seq(j).numFtCon5       = T.prop_bounce5(i);
            seq(j).totFt           = T.tot_ft(i);
            seq(j).SW_tot          = T.tot_sw(i);
            
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
end


%% Import data

if  do.reImportData || ~isfile([paths.data filesep 'SideDataPooled.mat'])
    
    iSeq = 1;
    
    for j = 1:length(dlc_files)
        
        % Find match to current sequence --------------------------------------
        
        iCurrent = nan;

        % Loop thru sequences for a match
        for i = 1:length(seq)
            if strcmp(seq(i).fName_side,dlc_files(j).name(1:23))
                iCurrent = i;
            end
        end

        % Check if there a match
        if isnan(iCurrent)
            warning([dlc_files(j).name(1:23) ...
                ' has no match to the sequences in Weights_experiments.xlsx'])

        % If if there is a match . . .
        else

            % get info on current sequence
            cSeq = seq(iCurrent);

            % Compile data ------------------------------------------------

            % Read current csv data
            %         T = readtable([paths.data filesep 'data' filesep 'side_view' filesep ...
            %                        dlc_files(j).name],'HeaderLines',2);
            T = readtable([paths.data filesep 'rawCSV' filesep ...
                dlc_files(j).name],'HeaderLines',2);

            % Add time vector
            T.t = [0:(1/cSeq.fps):((length(T.x)-1)/cSeq.fps)]';

            % Assume that y-value of the
            T.y = T.y-min(T.y);

            % Y-displacement
            dY = abs(diff(T.y));

            % Index for problematic frames
            iMessedUp = T.likelihood(1:end-1)<0.95 | dY>10;

            % Trim up to messed up data
            if sum(iMessedUp)>0
                iEnd            = find(iMessedUp,1,'first');
                T(iEnd:end,:)   = [];
            end

            % Ignore short sequences, store the rest
            if max(T.t)>30

                % Sequence info
%                 S(iSeq).fName_side = cSeq.fName_side;
%                 S(iSeq).fName_bot  = cSeq.fName_bot;
%                 S(iSeq).expType    = cSeq.expType;
%                 S(iSeq).addMass    = cSeq.addMass;
%                 S(iSeq).floatNum   = cSeq.floatNum;
%                 S(iSeq).indiv      = cSeq.indiv;
%                 S(iSeq).bodyMass   = cSeq.bodyMass;
%                 S(iSeq).calConst   = cSeq.calConst;
%                 S(iSeq).SW_percent = cSeq.SW_percent;
%                 S(iSeq).fps        = cSeq.fps;
%                 S(iSeq).trial      = cSeq.trial;
%                 S(iSeq).propBounce = cSeq.propBounce;
%                 S(iSeq).SW_tot     = cSeq.SW_tot;

                S(iSeq).fName_side  = cSeq.fName_side;
                S(iSeq).fName_bot   = cSeq.fName_bot;
                S(iSeq).expType     = cSeq.expType;
                S(iSeq).addMass     = cSeq.addMass;
                S(iSeq).floatNum    = cSeq.floatNum;
                S(iSeq).indiv       = cSeq.indiv;
                S(iSeq).bodyMass    = cSeq.bodyMass;
                S(iSeq).calConst    = cSeq.calConst;
                S(iSeq).SW_percent  = cSeq.SW_percent;
                S(iSeq).fps         = cSeq.fps;
                S(iSeq).trial       = cSeq.trial;
                S(iSeq).numFtCon1   = cSeq.numFtCon1;
                S(iSeq).numFtCon2   = cSeq.numFtCon2;
                S(iSeq).numFtCon3   = cSeq.numFtCon3;
                S(iSeq).numFtCon4   = cSeq.numFtCon4;
                S(iSeq).numFtCon5   = cSeq.numFtCon5;
                S(iSeq).totFt       = cSeq.totFt;
                S(iSeq).SW_tot      = cSeq.SW_tot;

                % Copy over coordinates
                S(iSeq).t       = T.t;
                S(iSeq).xRaw    = T.x;
                S(iSeq).yRaw    = 1280 - T.y;

                % Advance index
                iSeq      = iSeq + 1;
            end

            if 0
                figure;
                subplot(4,1,1)
                plot(T.t,T.xRaw)
                xlabel('t (s)'); ylabel('X')
                grid on
                title(['Exp type = ' cSeq.expType, ', ' cSeq.fName_side])

                subplot(4,1,2)
                plot(T.t,T.yRaw)
                xlabel('t (s)'); ylabel('Y')
                grid on

                subplot(4,1,3)
                plot(T.t(2:end),abs(diff(T.yRaw)))
                xlabel('t (s)'); ylabel('dY')
                grid on

                subplot(4,1,4)
                plot(T.t,T.likelihood)
                xlabel('t (s)'); ylabel('Likelihood')
                grid on
            end
        end
    end
    
    % Save
    save([paths.data filesep 'SideDataPooled'],'S'); 
    
end



%% Analyze audio sync

if do.anaAudioSync
    
% Loop thru sequences to be analyzed
for i = 1:length(seq)
    
    % Video paths
    vPath_side = [paths.vid filesep seq(i).dirName filesep 'side' ...
                  filesep seq(i).fName_side '.' seq(i).ext];
    vPath_bot  = [paths.vid filesep seq(i).dirName filesep 'bottom' ...
                  filesep seq(i).fName_bot '.' seq(i).ext];
    
    % Define where to save the data
    dPath_save = [paths.data filesep seq(i).dirName filesep ...
        'bottom' filesep 'matlabData2021' filesep seq(i).fName_bot];

    

    % Update status
    disp(['Analyzing audio sync for ' seq(i).fName_side ' and ' seq(i).fName_bot])
    
    % Find delay
    [delay,info] = audio_sync(vPath_side,vPath_bot);

    % Store results
    aud.delay = delay;
    aud.info  = info;

    % Number of points to trim from edges
    eTrim = 30;

    % Load body data (Body)
    load([paths.data filesep seq(i).dirName filesep 'bottom' filesep ...
        'matlabData2021' filesep seq(i).fName_bot filesep 'Body.mat'])

    % Load side kinematics data (S)
    load([paths.data filesep 'SideDataPooled_eventsManual2.mat']);

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

    % Extract current S
    S = S(iMatch);

    %
    t_B = (Body.frames'-min(Body.frames)) ./ seq(i).fps_bot ;
    t_S = S.t-min(S.t) + aud.delay;

    % Extract x-values, smooth
    x_S = smooth(S.xRaw - min(S.xRaw),50);
    x_B = smooth(Body.xCntr - min(Body.xCntr),50);

    % Trim edges
    t_B   = t_B(eTrim:(end-eTrim));
    t_S   = t_S(eTrim:(end-eTrim));
    x_B   = x_B(eTrim:(end-eTrim));
    x_S   = x_S(eTrim:(end-eTrim));

    % Adjust starting position
    x_B  = x_B - x_B(1);
    x_S  = x_S - x_S(1);

    % Trim data wrt time
    tMin    = max([min(t_B) min(t_S)]);
    tMax    = min([max(t_B) max(t_S)]);
    iB      = t_B>=tMin & t_B<=tMax;
    iS      = t_S>=tMin & t_S<=tMax;
    t_B     = t_B(iB);
    t_S     = t_S(iS);
    x_B     = x_B(iB);
    x_S     = x_S(iS);

    % Resample data by interpolation
    if seq(i).fps > seq(i).fps_bot
        x_B = interp1(t_B,x_B,t_S);
        t = t_S;
        fps = seq(i).fps;
    elseif seq(i).fps < seq(i).fps_bot
        x_S = interp1(t_S,x_S,t_B);
        t = t_B;
        fps = seq(i).fps_bot;
    end
    clear t_B t_S tMin tMax iB iS

    % Check lengths
    if length(x_S)~=length(x_B)
        error('These should be equal here')
    end

    % Normalize by range
    x_S = x_S ./range(x_S);
    x_B = x_B ./ range(x_B);

    % Find derivative
    s_S = diff(x_S)./diff(t);
    s_B = diff(x_B)./diff(t);

    % Normalize speeds
    s_S = (s_S-min(s_S)) ./ range(s_S);
    s_B = (s_B-min(s_B)) ./ range(s_B);

    [C,lags] = xcorr(s_B,s_S);

    iMax = find(C==max(C),1,'first');
    lag = lags(iMax)./fps;

    % Add lag to delay
    aud.delay = delay + lag;

    if 0
        figure
        subplot(2,1,1)
        plot(t+lag,x_S,'-',t,x_B,'-')
        grid on
        subplot(2,1,2)
        plot(t(2:end)+lag,s_S,'-',t(2:end),s_B,'-')
        grid on

        ttt = 3;
    end

    % Report info on audio
    disp(['     ' info])

    % Video objects
    v_side = VideoReader(vPath_side);
    v_bot  = VideoReader(vPath_bot);

    % Store frame rates
    frRate.side = v_side.FrameRate;
    frRate.bot  = v_bot.FrameRate;

    % Save
    save([dPath_save filesep 'audio_delay'],'aud');   
    save([dPath_save filesep 'frame_rate'],'frRate');   
end

end %do.anaAudioSync


%% Event finding

if do.reImportData || ~isfile([paths.data filesep 'SideDataPooled_events.mat'])
    
    % Load 'S' structure
    load([paths.data filesep 'SideDataPooled.mat'])
    
    % Acceptable deviation from bounce amplitude and period
    yThresh  = 0.8;
    tThresh  = 0.8;
    
    % Loop trhu sequences
    for i = 1:length(S)

        % Rotate wrt central axis
        [xRot,yRot] = rotCoords(S(i).xRaw,S(i).yRaw);
        
        % Apply calibration constant
        xRot = xRot ./ S(i).calConst;
        yRot = yRot ./ S(i).calConst;

        % Find when the body hits the ground
        [tEvent,yEvent,xS,yS] = findPeaker(S(i).t,xRot,yRot);
        
%         yEvent   = yEvent  .* S(i).calConst;
%         
        % Forward speed
        xSpd = abs(diff(xS))./diff(S(i).t);
        xSpd = [xSpd; xSpd(end)];
        
        % Mean period
        tau = mean(diff(tEvent));
        
        % Individual step period
        sPeriod = diff([0;tEvent]);
        
        % Indicies for landings in the last 3-quarters
        idx = tEvent>range(S(i).t) * 0.25;
        
        % Loop trhu events
        for j = 1:length(tEvent)
            
            % Current interval, and y range
            if j==1
                idx = S(i).t<=tEvent(j);
            else
                idx = (S(i).t>tEvent(j-1)) & (S(i).t<=tEvent(j));
            end
            
            if sum(idx)>3
                % Range in y over bounce
                yRange(j,1) = max(S(i).yRaw(idx)) - min(S(i).yRaw(idx));
                
                % Mean speed over bounce
                meanSpd(j,1) = nanmean(xSpd(idx));
            else
                yRange(j,1) = 0;
                meanSpd(j,1) = nan;
            end
            
        end

        % Identify non-bouncing events
        iCrawl = (yRange < yThresh*mean(yRange)) | ...
                 (sPeriod < tau*tThresh) | (sPeriod > (2-tThresh)*tau);

        % Time to first bounce
        tFirstBounce = tEvent(find(~iCrawl,1,'first'));
        if isempty(tFirstBounce)
            tFirstBounce = nan;
        end

        % Store
        S(i).xS                 = xS;
        S(i).yS                 = yS;
        S(i).xSpd               = xSpd;
        S(i).tLand              = tEvent;
        S(i).yLand              = yEvent;
        S(i).iCrawl             = iCrawl;

        % Visual check
        if 0
            subplot(2,1,1)
            plot(S(i).t,S(i).yS,'k-');hold on
            plot(S(i).tLand(~S(i).iCrawl),S(i).yLand(~S(i).iCrawl),'ro')
            plot(S(i).tLand(S(i).iCrawl),S(i).yLand(S(i).iCrawl),'r+')
            hold off
            xlabel('t (s)'); ylabel('Y')
            grid on
        end

        clear tEvent yEvent xS yS iCrawl idx yRange meanSpd xSpd sPeriod
    end
    
    % Save
    save([paths.data filesep 'SideDataPooled_events'],'S');   
end


%% Manually annotate data

if do.manAnnotate

% Load data (S)
if isfile([paths.data filesep 'SideDataPooled_eventsManual.mat'])
    load([paths.data filesep 'SideDataPooled_eventsManual'])
else
    % Load pooled data with events identified automatically
    load([paths.data filesep 'SideDataPooled_events'])
end

f = figure;

disp(' ')
% disp('Select points, press RETURN (or double-click, or right-click) when done')
disp('Select a time value to modify, press RETURN when done with seq')
disp('   Press z to zoom (or unzoom)')
disp(' ')

xL = []; h = [];

% Loop trhu sequences
for i = 1:length(S)

    % Check of annotated
    if ~isfield(S(i),'annotated') || isempty(S(i).annotated)

        while true
        
            % Plot time series
            figure(f)
            ax = subplot(2,1,1);
            plot(S(i).t,S(i).yS,'k-');hold on
            h1 = plot(S(i).tLand(~S(i).iCrawl),S(i).yLand(~S(i).iCrawl),'ro');
            h2 = plot(S(i).tLand(S(i).iCrawl),S(i).yLand(S(i).iCrawl),'r+');
            title([S(i).fName_side])
            hold off

            if ~isempty(xL)
                set(ax,'XLim',xL);
            end

            % Prompt for coordinates
            [x,~,b] = ginput(1);

            if ~isempty(x) && b==122

                if isempty(h) 
                    h = zoom;
                    h.Motion = 'horizontal';
                    h.Enable = 'on';

                    disp('Unclick the magnifying glass when done')

                    waitfor(h,'Enable','off')

                    xL = get(ax,'XLim');
                else
                    xL = [];
                    h  = [];
                end

            elseif ~isempty(x)

                ButtonName = questdlg('What type of point?', ...
                         ' ', 'Switch type', 'New bounce peak', 'Nevermind', ...
                         'Switch type');

                switch ButtonName

                case 'Switch type'
                    
                    % Difference in time
                    tDiff = abs(x-S(i).tLand);
                    
                    % Closest peak
                    iClose = find(tDiff==min(tDiff),1,'first');

                    % Toggle crawl/bounce
                    S(i).iCrawl(iClose) = abs(S(i).iCrawl(iClose)-1);

                case 'New bounce peak'
                    
                    % Indicies around x
                    iBefore = find(S(i).tLand<x,1,'last');
                    iAfter  = find(S(i).tLand>x,1,'first');

                    % Interpolate for y
                    yLand = interp1(S(i).t,S(i).yS,x);

                    % Store new point
                    if isnan(yLand)
                        warning('You need to select a point within the time series')
                        
                    elseif isempty(iBefore)
                        S(i).tLand  = [x; S(i).tLand];
                        S(i).yLand  = [yLand; S(i).yLand];
                        S(i).iCrawl = [false; S(i).iCrawl];
                    elseif isempty(iAfter)
                        S(i).tLand   = [S(i).tLand; x];
                        S(i).yLand   = [S(i).yLand; yLand];
                        S(i).iCrawl  = [S(i).iCrawl; false];
                    else
                        S(i).tLand   = [S(i).tLand(1:iBefore); x; ...
                                        S(i).tLand(iAfter:end)];
                        S(i).yLand   = [S(i).yLand(1:iBefore); yLand; ...
                                        S(i).yLand(iAfter:end)];
                        S(i).iCrawl  = [S(i).iCrawl(1:iBefore); false; ...
                                        S(i).iCrawl(iAfter:end)];
                    end

                case 'Nevermind'
                    % Do nothing

                end % switch

            % If RETURN pressed . . .
            else
                ButtonName = questdlg('Save data?', ...
                         ' ', 'Yes', 'No, stop for now', 'Yes');
                switch ButtonName
                    case 'Yes'
                        S(i).annotated = 1;
                        save([paths.data filesep 'SideDataPooled_eventsManual'],'S'); 
                        break

                    case 'No, stop for now'
                        return
                end
            end
        end
    end
end
end %manAnnotate


%% Calculate bounce kinematics

% Run this code, if it hasn't been already
if ~isfile([paths.data filesep 'SideDataPooled_eventsManual2.mat'])

% Load S 
load([paths.data filesep 'SideDataPooled.mat'])
Sup = S;
clear S
   
% Load 'S' structure
load([paths.data filesep 'SideDataPooled_eventsManual.mat'])

% Loop trhu sequences
for i = 1:length(S)
    
    % Find match to Sup
    iMatch = [];
    for j = 1:length(Sup)
        if strcmp(Sup(j).fName_side,S(i).fName_side)
            iMatch = j;
            break
        end
    end
    if isempty(iMatch), error('No match here');end

    % Update fields in S
%     S(i).trial      = Sup(iMatch).trial;
%     S(i).t          = Sup(iMatch).t;
%     S(i).fps        = Sup(iMatch).fps;
%     S(i).propBounce = Sup(iMatch).propBounce;
%     S(i).SW_tot     = Sup(iMatch).SW_tot;
    
    S(i).trial       = Sup(iMatch).trial;
    S(i).t           = Sup(iMatch).t;
    S(i).fps         = Sup(iMatch).fps;
    S(i).numFtCon1   = Sup(iMatch).numFtCon1;
    S(i).numFtCon2   = Sup(iMatch).numFtCon2;
    S(i).numFtCon3   = Sup(iMatch).numFtCon3;
    S(i).numFtCon4   = Sup(iMatch).numFtCon4;
    S(i).numFtCon5   = Sup(iMatch).numFtCon5;
    S(i).totFt       = Sup(iMatch).totFt;
    S(i).SW_tot      = Sup(iMatch).SW_tot;
    

%     if isempty(S(i).propBounce)
%         S(i).propBounce = nan;
%     end
    
    
      if isempty(S(i).numFtCon1)
          S(i).numFtCon1 = nan;
      end
      
      if isempty(S(i).numFtCon2)
          S(i).numFtCon2 = nan;
      end
      
      if isempty(S(i).numFtCon3)
          S(i).numFtCon3 = nan;
      end
      
      if isempty(S(i).numFtCon4)
          S(i).numFtCon4 = nan;
      end

      if isempty(S(i).numFtCon5)
          S(i).numFtCon5 = nan;
      end

    % Forward speed
    xSpd = abs(diff(S(i).xS))./diff(S(i).t);
    xSpd = [xSpd; xSpd(end)];
    
    % Mean period
    tau = mean(diff(S(i).tLand));
    
    % Individual step period
    sPeriod = diff([0;S(i).tLand]);
    
    % 
    
    % Loop trhu events
    for j = 1:length(S(i).tLand)
        
        % Current interval, and y range
        if j==1
            idx = S(i).t<=S(i).tLand(j);
        else
            idx = (S(i).t>S(i).tLand(j-1)) & (S(i).t<=S(i).tLand(j));
        end
        
        % If the interval has sufficient duration . . .
        if sum(idx)>3
            % Range in y over bounce
            yRange(j,1) = max(S(i).yS(idx)) - min(S(i).yS(idx));
            
            % Mean speed over bounce
            meanSpd(j,1) = nanmean(xSpd(idx));

        % Skip, if too short
        else
            yRange(j,1)  = 0;
            meanSpd(j,1) = nan;
        end        
    end

    % Time to first bounce
    tFirstBounce = S(i).tLand(find(~S(i).iCrawl,1,'first'));
    if isempty(tFirstBounce)
        tFirstBounce = nan;
    end

    % Store
    S(i).crawlSpd           = nanmean(meanSpd(S(i).iCrawl));
    S(i).bounceSpd          = nanmean(meanSpd(~S(i).iCrawl));
    S(i).yBounceAmp         = mean(yRange(~S(i).iCrawl)); 
    S(i).stepPeriod         = mean(sPeriod(~S(i).iCrawl));   
    S(i).propTimeBounce     = sum(sPeriod(~S(i).iCrawl)) / range(S(i).t);
    S(i).meanBouncePeriod   = mean(sPeriod(~S(i).iCrawl));
    S(i).tFirstBounce       = tFirstBounce;

    % Visual check
    if 0
        subplot(2,1,1)
        plot(S(i).t,S(i).yS,'k-');hold on
        plot(S(i).tLand(~S(i).iCrawl),S(i).yLand(~S(i).iCrawl),'ro')
        plot(S(i).tLand(S(i).iCrawl),S(i).yLand(S(i).iCrawl),'r+')
        hold off
        xlabel('t (s)'); ylabel('Y')
        grid on
    end

    clear tEvent yEvent xS yS iCrawl idx yRange meanSpd xSpd sPeriod
end

% Save S
save([paths.data filesep 'SideDataPooled_eventsManual2.mat'],'S')

else

% Load S structure
load([paths.data filesep 'SideDataPooled_eventsManual2.mat'])

end


%% Visualize event finding

if do.visEvents
    
    makeEventPlots = 1;
    
    if ~exist('S','var')
        % Load 'S' structure
        load([paths.data filesep 'SideDataPooled_eventsManual2.mat'])
    end
    
    % Figure parameters
    nPanels        = 6;
    iPanel         = 1;
    
    % Acceptable deviation from bounce amplitude and period
    yThresh  = 0.8;
    tThresh  = 0.8;
    
    makeEventPlots = 1;
    
    if makeEventPlots
        f = figure;
    end
    
    % Loop trhu sequences
    for i = 1:length(S)
        
         figure(f)
% figure
        subplot(nPanels,1,iPanel)
%         plot(S(i).t,S(i).yRaw .* S(i).calConst); hold on
        plot(S(i).t,S(i).yS,'k-');hold on
        plot(S(i).tLand(~S(i).iCrawl),S(i).yLand(~S(i).iCrawl),'ro')
        plot(S(i).tLand(S(i).iCrawl),S(i).yLand(S(i).iCrawl),'r+')
        hold off
        xlabel('t (s)'); ylabel('Y')
        grid on
        title(['ExpType = ' S(i).expType ', ' S(i).fName_side])
        
        iPanel = iPanel + 1;
        
        if iPanel>nPanels && i<length(S)
            iPanel = 1;
            f = figure;
        end
        
        clear tEvent yEvent xS yS iCrawl idx yRange meanSpd xSpd sPeriod
    end
    
        for i = 1:length(S)
        
         figure(f)
% figure
        subplot(nPanels,1,iPanel)
%         plot(S(i).t,S(i).yRaw .* S(i).calConst); hold on
% smooth x and y ? knock off last value. and have 1/2 of each time pt
        spd        = sqrt(diff(S(i).xS).^2 + diff(S(i).yS)).^2./diff(S(i).t);
        spdSmooth  = smooth(spd,100);
        timeMinus1 = (S(i).t(1:end-1,:));
        plot(timeMinus1,spdSmooth,'b-');hold on
        %plot(S(i).tLand(~S(i).iCrawl),S(i).yLand(~S(i).iCrawl),'ro')
        %plot(S(i).tLand(S(i).iCrawl),S(i).yLand(S(i).iCrawl),'r+')
        hold off
        xlabel('t (s)'); ylabel('Velocity')
        grid on
        title(['ExpType = ' S(i).expType ', ' S(i).fName_side])
        
        iPanel = iPanel + 1;
        
        if iPanel>nPanels && i<length(S)
            iPanel = 1;
            f = figure;
        end
        
        clear tEvent yEvent xS spd spdSmooth timeMinus1 yS iCrawl idx yRange meanSpd xSpd sPeriod
        end
end


%% Visualize details on each sequence

if do.visSeqs
    
    if ~exist('S','var')
        % Load 'S' structure
        load([paths.data filesep 'SideDataPooled_eventsManual2.mat'])
    end

    % Acceptable deviation from bounce amplitude and period
    yThresh  = 0.8;
    tThresh  = 0.8;
    
    disp('Press any key to advance to next graph')

    f2 = figure;
    
    % Loop trhu sequences
    for i = 1:length(S)
        
        % Plots
        figure(f2)
        subplot(3,1,1)
        plot(S(i).t,S(i).yRaw .* S(i).calConst); hold on
        plot(S(i).t,S(i).yS,'k-');
        plot(S(i).tLand(~S(i).iCrawl),S(i).yLand(~S(i).iCrawl),'ro')
        plot(S(i).tLand(S(i).iCrawl),S(i).yLand(S(i).iCrawl),'r+')
        yL = ylim;
        for j = 1:length(S(i).tLand)
            plot(S(i).tLand(j).*[1 1],yL,'k-')
        end
        hold off
        xlabel('t (s)'); ylabel('Y')
        
        subplot(3,1,2)
        plot(S(i).t,S(i).xRaw.* S(i).calConst)
        hold on
        xlabel('t (s)'); ylabel('X (cm)')
        yL = ylim;
        for j = 1:length(S(i).tLand)
            plot(S(i).tLand(j).*[1 1],yL,'k-')
        end
        hold off
        
        subplot(3,1,3)
        plot(S(i).t,S(i).xSpd)
        xlabel('t (s)'); ylabel('x-spd (cm)')
        hold on
        yL = ylim;
        for j = 1:length(S(i).tLand)
            plot(S(i).tLand(j).*[1 1],yL,'k-')
        end
        hold off
        
        title(['ExpType = ' S(i).expType ', ' S(i).fName_side])
        
        pause
        
        clear tEvent yEvent xS yS iCrawl idx yRange meanSpd xSpd sPeriod
    end
end

%% Plot results for treatment within individuals

if do.indivPlots
    
if ~exist('S','var')
    % Load 'S' structure
    load([paths.data filesep 'SideDataPooled_eventsManual2.mat'])
end

% Loop thru sequences, collect data
% for i = 1:length(S)
% 
%     indiv(i,1)          = S(i).indiv;
%     meanBPeriod(i,1)    = S(i).meanBouncePeriod; 
%     SW_percent(i,1)     = S(i).SW_percent;
%     beat_spd(i,1)       = S(i).bounceSpd;
%     crawl_spd(i,1)      = S(i).crawlSpd;
%     tFirstB(i,1)        = S(i).tFirstBounce;
%     yAmpB(i,1)          = S(i).yBounceAmp;
%     propBounce(i,1)     = S(i).propBounce;
% end

 for i = 1:length(S)
 
     indiv(i,1)       = S(i).indiv;
     meanBPeriod(i,1) = S(i).meanBouncePeriod; % column vector
     SW_percent(i,1)  = S(i).SW_percent;
     beat_spd(i,1)    = S(i).bounceSpd;
     crawl_spd(i,1)   = S(i).crawlSpd;
     tFirstB(i,1)     = S(i).tFirstBounce;
     yAmpB(i,1)       = S(i).yBounceAmp;
     numFtCon(i,1)    = S(i).numFtCon1;% number of feet in contact in PS1
     numFtCon(i,2)    = S(i).numFtCon2;
     numFtCon(i,3)    = S(i).numFtCon3;
     numFtCon(i,4)    = S(i).numFtCon4;
     numFtCon(i,5)    = S(i).numFtCon5;% number of feet in contact in PS5
     totFt(i,1)       = S(i).totFt; % total number of feet of the sea star
   
 end
 
 
indNums = unique(indiv);


f = figure;

% Line colors
lClr = lines(length(indNums));

% Loop thru individuals
for i = 1:length(indNums)

    % Index of values for current individual
    iIdx = indiv==indNums(i);

    % Unique SW values
    SWs = unique(SW_percent(iIdx));

    % Loop thru SW values
    for j = 1:length(SWs) % each row for the matrix and std of that row

         % Index for current values
         idx = iIdx & SW_percent==SWs(j);
         
         % Log data
%          BPeriod_mean(j)    = nanmean(meanBPeriod(idx));
%          BPeriod_std(j)     = nanstd(meanBPeriod(idx));
%          Bspd_mean(j)       = nanmean(beat_spd(idx));
%          Bspd_std(j)        = nanstd(beat_spd(idx));
%          BtFirst_mean(j)    = nanmean(tFirstB(idx));
%          BtFirst_std(j)     = nanstd(tFirstB(idx));
%          Bamp_mean(j)       = nanmean(yAmpB(idx));
%          Bamp_std(j)        = nanstd(yAmpB(idx));
%          Bcrawlspd_mean(j)  = nanmean(crawl_spd(idx));
%          Bcrawlspd_std(j)   = nanstd(crawl_spd(idx));
%          Bprop_mean(j)      = nanmean(propBounce(idx));
%          Bprop_std(j)       = nanstd(propBounce(idx));
         
         BPeriod_mean(j)    = nanmean(meanBPeriod(idx));
         BPeriod_std(j)     = nanstd(meanBPeriod(idx));
         Bspd_mean(j)       = nanmean(beat_spd(idx));
         Bspd_std(j)        = nanstd(beat_spd(idx));
         BtFirst_mean(j)    = nanmean(tFirstB(idx));
         BtFirst_std(j)     = nanstd(tFirstB(idx));
         Bamp_mean(j)       = nanmean(yAmpB(idx));
         Bamp_std(j)        = nanstd(yAmpB(idx));
         Bcrawlspd_mean(j)  = nanmean(crawl_spd(idx));
         Bcrawlspd_std(j)   = nanstd(crawl_spd(idx));
         
         numFtCon_meanNaN   = nanmean(numFtCon(idx,1:5),2);% mean of the tf in contact in 5 PS
         BpropNaN           = numFtCon_meanNaN./totFt(idx);% with NaN's, proportion of feet
         Bprop(j)           = BpropNaN(~isnan(BpropNaN));
         
         
         numFtConNaN        = numFtCon(idx,1:5);
         numFtConCur        = numFtConNaN(~isnan(numFtConNaN));% current number of feet in contact
         Bprop_std(j)       = std(numFtConCur)./nanmean(totFt(idx));% std of numFtCon divided by total number of ft w/o NaN
    
         
         
         %S(idx).Bprop         = Bprop(j);
         %S(idx).Bprop_std     = Bprop_std(j);


         
         clear numFtCon_meanNaN BpropNaN propBNum numFtConNaN numFtConCur
       
    end

    ttt = 3;
    
    % Bounce period
    subplot(3,2,1)
    ePlot(SWs,BPeriod_mean,BPeriod_std,lClr(i,:),lClr(i,:))
    ylabel('Bounce period')
    xlabel('Percent submerged weight')
    axis square
    hold on
    xlim([40 175])

    % Bounce speed
    subplot(3,2,2)
    ePlot(SWs,Bspd_mean,Bspd_std,lClr(i,:),lClr(i,:))
    ylabel('Bounce speed (m/s)')
    xlabel('Percent submerged weight')
    axis square
    hold on
    xlim([40 175])

    % Time to first bounce
    subplot(3,2,3)
    ePlot(SWs,BtFirst_mean,BtFirst_std,lClr(i,:),lClr(i,:))
    ylabel('Time to first bounce (s)')
    xlabel('Percent submerged weight')
    axis square
    hold on
    xlim([40 175])

    % Boucne amplitude
    subplot(3,2,4)
    ePlot(SWs,Bamp_mean,Bamp_std,lClr(i,:),lClr(i,:))
    ylabel('Y Amplitude')
    xlabel('Percent submerged weight')
    axis square
    hold on
    xlim([40 175])

    % Crawl speed
    subplot(3,2,5)
    ePlot(SWs,Bcrawlspd_mean,Bcrawlspd_std,lClr(i,:),lClr(i,:))
    ylabel('Crawl speed')
    xlabel('Percent submerged weight')
    axis square
    hold on
    xlim([40 175])

    % Proportion of tube feet
    subplot(3,2,6)
    ePlot(SWs,Bprop,Bprop_std,lClr(i,:),lClr(i,:))
    %ePlot(SWs,Bprop_mean,Bprop_std,lClr(i,:),lClr(i,:))
    ylabel('Proportion of feet in pwr stroke')
    xlabel('Percent submerged weight')
    axis square
    hold on
    xlim([40 175])
    
    clear SWs B*
end

ttt = 3; % check std 

end %do.indivPlots

%% Plot effects of trial number

if do.trialPlot


if ~exist('S','var')
    % Load 'S' structure
    load([paths.data filesep 'SideDataPooled_eventsManual2.mat'])
end

% Loop thru sequences, collect data
for i = 1:length(S)
    indiv(i,1)          = S(i).indiv;
    trial(i,1)          = S(i).trial;
    meanBPeriod(i,1)    = S(i).meanBouncePeriod; 
    SW_percent(i,1)     = S(i).SW_percent;
    beat_spd(i,1)       = S(i).bounceSpd;
    crawl_spd(i,1)      = S(i).crawlSpd;
    tFirstB(i,1)        = S(i).tFirstBounce;
    yAmpB(i,1)          = S(i).yBounceAmp;
    expType(i,1)        = S(i).expType;
end

indNums = unique(indiv);

expts = ['f','c','w'];

f = figure;
lClr = lines(length(trial));


% Step trhu individuals
for i = 1:length(indNums)
    
    % Loop thru experiment types
    for k = 1:3

        % Index of values for current individual
        iIdx = indiv==indNums(i) & expType==expts(k);

        % Unique trial numbers
        trials = unique(trial(iIdx));

        % Loop thru trial numbers
        for j = 1:length(trials)

            % Index for current values
            idx = iIdx & trial==trials(j);

%             BPeriod_mean(j)  = nanmean(meanBPeriod(idx));
%             BPeriod_std(j)   = nanstd(meanBPeriod(idx));
%             Bspd_mean(j)     = nanmean(beat_spd(idx));
%             Bspd_std(j)      = nanstd(beat_spd(idx));
%             BtFirst_mean(j)  = nanmean(tFirstB(idx));
%             BtFirst_std(j)   = nanstd(tFirstB(idx));
%             Bamp_mean(j)     = nanmean(yAmpB(idx));
%             Bamp_std(j)      = nanstd(yAmpB(idx));
            
            BPeriod_mean(j)  = nanmean(meanBPeriod(idx));
            BPeriod_std(j)   = nanstd(meanBPeriod(idx));
            Bspd_mean(j)     = nanmean(beat_spd(idx));
            Bspd_std(j)      = nanstd(beat_spd(idx));
            BtFirst_mean(j)  = nanmean(tFirstB(idx));
            BtFirst_std(j)   = nanstd(tFirstB(idx));
            Bamp_mean(j)     = nanmean(yAmpB(idx));
            Bamp_std(j)      = nanstd(yAmpB(idx));
           
        end

        % Bounce period
        subplot(4,3,k)
        ePlot(trials,BPeriod_mean,BPeriod_std,lClr(i,:),lClr(i,:))
        ylabel('Bounce period')
        xlabel('Trial number')
        axis square
        hold on
        title([expts(k) ' type experiment'])
        xlim([0 200])
        ylim([2 10])

        % Bounce speed
        subplot(4,3,k+3)
        ePlot(trials,Bspd_mean,Bspd_std,lClr(i,:),lClr(i,:))
        ylabel('Bounce speed (m/s)')
        xlabel('Trial number')
        axis square
        hold on
        xlim([0 200])
        ylim([0 3.5e-3])

        % Time to first bounce
        subplot(4,3,k+6)
        ePlot(trials,BtFirst_mean,BtFirst_std,lClr(i,:),lClr(i,:))
        ylabel('Time to first bounce (s)')
        xlabel('Trial number')
        axis square
        hold on
        xlim([0 200])
        ylim([0 180])

        % Boucne amplitude
        subplot(4,3,k+9)
        ePlot(trials,Bamp_mean,Bamp_std,lClr(i,:),lClr(i,:))
        ylabel('Y Amplitude')
        xlabel('Trial number')
        axis square
        hold on
        xlim([0 200])
        ylim([0 2.5e-3])

        clear SWs B*
    end
end

if 1
    figure
%     subplot(2,2,1)
    idx = expType=='f';
    scatter(crawl_spd(idx),beat_spd(idx),50,'b','filled')
    hold on
    idx = expType=='c';
    scatter(crawl_spd(idx),beat_spd(idx),50,'r','filled')
    idx = expType=='w';
    scatter(crawl_spd(idx),beat_spd(idx),50,'g','filled')
    hold off
    axis equal
    xlabel('Crawl spd'); ylabel('Bounce Spd')
    legend('f','c','w')
end



end %do.trialPlot



%% Plots of summary data

if do.summaryPlots
    
    if ~exist('S','var')
        % Load 'S' structure
        load([paths.data filesep 'SideDataPooled_eventsManual2.mat'])
    end
    
    bouncePeriod = [];
    propBounce   = [];
    expType      = [];
    meanSpd      = [];
    bounceSpd    = [];
    crawlSpd     = [];
    SW           = [];
    
    for i = 1:length(S)
        
        % Boxplot categories
        if S(i).expType=='w' && S(i).addMass==2.09
            expType = [expType; 5];
            
        elseif S(i).expType=='w' && S(i).addMass==1.04
            expType = [expType; 4];
            
        elseif S(i).expType=='c'
            expType = [expType; 3];
            
        elseif S(i).expType=='f' && S(i).floatNum==1
            expType = [expType; 2];
            
        elseif S(i).expType=='f' && S(i).floatNum==2
            expType = [expType; 1];
            
        else
            error('Do not recognize experiment type');
        end
        
        % Store submerged weight for scatterplot
%          SW =  [SW; (S(i).SW_percent./100) * (S(i).bodyMass./1000.*9.81)];
        SW =  [SW; (S(i).SW_percent./100)];
        %TODO: This will probably need some adjustment when we run the
        %analysis for real
         
         
        % Store measurements
%         meanSpd         = [meanSpd;       mean(S(i).xSpd)];
%         bounceSpd       = [bounceSpd;     S(i).bounceSpd];
%         crawlSpd        = [crawlSpd;      S(i).crawlSpd];
%         
%         bouncePeriod    = [bouncePeriod;  S(i).meanBouncePeriod];
%         propBounce      = [propBounce;    S(i).propBounce];
        
        meanSpd         = [meanSpd;       mean(S(i).xSpd)];
        bounceSpd       = [bounceSpd;     S(i).bounceSpd];
        crawlSpd        = [crawlSpd;      S(i).crawlSpd];
        
        bouncePeriod    = [bouncePeriod;  S(i).meanBouncePeriod];
        propBounce      = [propBounce;    S(i).propBounce];
        
    end
    
    % Box plots
    figure;
    
    subplot(1,3,1)
    boxplot(meanSpd./10,expType)
    ylabel('Mean speed (cm/s)')
    set(gca,'TickDir','out')
     ylim([20 140])
     axis square
    
    subplot(1,3,2)
    boxplot(bounceSpd./10,expType)
    ylabel('Bounce speed (cm/s)')
    set(gca,'TickDir','out')
     ylim([20 140])
     axis square
    
    subplot(1,3,3)
    boxplot(crawlSpd./10,expType)
    ylabel('Crawl speed (cm/s)')
    set(gca,'TickDir','out')
     ylim([20 140])
     axis square
    
    
    figure;
    subplot(1,2,2)
    boxplot(bouncePeriod,expType)
    ylabel('Bounce Period')
    set(gca,'TickDir','out')
    axis square
    
    subplot(1,2,1)
    boxplot(propBounce,expType)
    ylabel('Proportion of time bouncing')
    set(gca,'TickDir','out')
    axis square
    
    
    % Scatterplots
    
    mSize = 75; % Marker size
    
    figure

    % Linear fit to bpunce period data
    cMean     = polyfit(SW,meanSpd./10,1);
    cBounce   = polyfit(SW,bounceSpd./10,1);
    idx       = ~isnan(crawlSpd);
    cCrawl    = polyfit(SW(idx),crawlSpd(idx)./10,1);
    
    subplot(1,3,1)
    h = scatter(SW,meanSpd./10,mSize, 'k','filled');
    hold on
    plot(SW,polyval(cMean,SW),'k-')
    hold off
    axis square
    xlabel('Submerged weight (N)')
    ylabel('Mean speed (cm/s)');
    set(gca,'TickDir','out')
    
    subplot(1,3,2)
    h = scatter(SW,bounceSpd./10,mSize, 'k','filled');
    hold on
    plot(SW,polyval(cBounce,SW),'k-')
    hold off
    axis square
    xlabel('Submerged weight (N)')
    ylabel('Bounce speed (cm/s)');
    set(gca,'TickDir','out')
    
    subplot(1,3,3)
    h = scatter(SW,crawlSpd./10,mSize, 'k','filled');
    hold on
    plot(SW,polyval(cCrawl,SW),'k-')
    hold off
    axis square
    xlabel('Submerged weight (N)')
    ylabel('Crawl speed (cm/s)');
    set(gca,'TickDir','out')
    
    figure

    % Linear fit to bpunce period data
    c = polyfit(SW,bouncePeriod,1);
    
    subplot(1,2,1)
    h = scatter(SW,bouncePeriod,mSize, 'k','filled');
    hold on
    plot(SW,polyval(c,SW),'k-')
    hold off
    axis square
    xlabel('Submerged weight (N)')
    ylabel('Bounce period (s)');
    set(gca,'TickDir','out')
    
    subplot(1,2,2)
    h = scatter(SW,propBounce,mSize, 'k','filled');
    axis square
    xlabel('Submerged weight (N)')
    ylabel('Proportion of time bouncing');
    set(gca,'TickDir','out')
    
end


%% Run statistics

if do.stats

% Put data into vectors
%for i = 1:length(S)

   % SW(i,1)             = S(i).SW_percent;
   % SW_tot(i,1)         = S(i).SW_tot;
   % indiv(i,1)          = S(i).indiv;
   % trial(i,1)          = S(i).trial;
   % propFtBounce(i,1)   = S(i).propBounce;
   % propTimeBounce(i,1) = S(i).propTimeBounce;
   % stepPeriod(i,1)     = S(i).stepPeriod;
   % bounceSpd(i,1)      = S(i).bounceSpd;
   % crawlSpd(i,1)       = S(i).crawlSpd;
   % tBounce(i,1)        = S(i).tFirstBounce;
   % bounceAmp(i,1)      = S(i).yBounceAmp;
%end

%Put data into vectors
for i = 1:length(S)

   SW(i,1)             = S(i).SW_percent;
   SW_tot(i,1)         = S(i).SW_tot;
   indiv(i,1)          = S(i).indiv;
   trial(i,1)          = S(i).trial;
   numFtCon(i,1)       = S(i).numFtCon1;% number of feet in contact in PS1
   numFtCon(i,2)       = S(i).numFtCon2;
   numFtCon(i,3)       = S(i).numFtCon3;
   numFtCon(i,4)       = S(i).numFtCon4;
   numFtCon(i,5)       = S(i).numFtCon5;% number of feet in contact in PS5
   totFt(i,1)          = S(i).totFt; % total number of feet of the sea star propTimeBounce(i,1) = S(i).propTimeBounce;
   stepPeriod(i,1)     = S(i).stepPeriod;
   bounceSpd(i,1)      = S(i).bounceSpd;
   crawlSpd(i,1)       = S(i).crawlSpd;
   tBounce(i,1)        = S(i).tFirstBounce;
   bounceAmp(i,1)      = S(i).yBounceAmp;
end


disp('-------------------------------------------------------------------')
disp('STEP PERIOD -------------------------------------------------------')
disp('-------------------------------------------------------------------')
idx = indiv>1;
glme = mixedModel(stepPeriod,SW,SW_tot,trial,indiv,'Step_period',...
                  'SW_percent','SW_tot');

disp('-------------------------------------------------------------------')
disp('BOUNCE SPEED ------------------------------------------------------')
disp('-------------------------------------------------------------------')
glme = mixedModel(bounceSpd,SW,SW_tot,trial,indiv,'Bounce_spd','SW_percent','SW_tot');

disp('-------------------------------------------------------------------')
disp('CRAWL SPEED -------------------------------------------------------')
disp('-------------------------------------------------------------------')
glme = mixedModel(crawlSpd,SW,SW_tot,trial,indiv,'Crawl_spd','SW_percent','SW_tot');

disp('-------------------------------------------------------------------')
disp('TIME TO FIRST BOUNCE ----------------------------------------------')
disp('-------------------------------------------------------------------')
glme = mixedModel(stepPeriod,SW,SW_tot,trial,indiv,'tFirst_bounce','SW_percent','SW_tot');

disp('-------------------------------------------------------------------')
disp('BOUNCE AMPLITUDE --------------------------------------------------')
disp('-------------------------------------------------------------------')
glme = mixedModel(bounceAmp,SW,SW_tot,trial,indiv,'Bounce_amp','SW_percent','SW_tot');

numFtConMean = nanmean(numFtCon,2); % mean of the 5 powerstrokes collected per indiv per Submerged weight
propFtBounce = numFtConMean./totFt; % Divide that mean by the total number of feet per sea star

disp('-------------------------------------------------------------------')
disp('PROPORTION OF FEET IN BOUNCES -------------------------------------')
disp('-------------------------------------------------------------------')
idx = ~isnan(propFtBounce);
glme = mixedModel(propFtBounce(idx),SW(idx),SW_tot(idx),trial(idx),indiv(idx),'Prop_feet',...
     'SW_percent','SW_tot');



% disp('-------------------------------------------------------------------')
% disp('PROPORTION OF FEET IN BOUNCES -------------------------------------')
% disp('-------------------------------------------------------------------')
% idx = ~isnan(propFtBounce);
% glme = mixedModel(propFtBounce(idx),SW(idx),SW_tot(idx),trial(idx),indiv(idx),'Prop_feet',...
%     'SW_percent','SW_tot');




end %stats

%% FUNCTIONS --------------------------------------


function glme = mixedModel(Y,X1,X2,trial,indiv_num,Yname,X1name,X2name)
% Runs a generalized linear mixed-effects model
% Y is the dependent variale
% X1 is a continuous independent variable
% X2-X3 are categorical variables
% indiv - individual number (random effect)

% Factors 

% Package school numbers into cells
for i = 1:length(indiv_num)
    indiv{i,1}   = num2str(indiv_num(i));
%     trial{i,1}   = num2str(tr_num(i));
end
    
% School number is a random variable
% trial = categorical(trial);
indiv = categorical(indiv);

% Table with data (2 independent continuous variables_
%tbl = table(Y,X1,X2,trial,indiv, 'VariableNames',{Yname,X1name,X2name,'trial','indiv'});

tbl = table(Y,X1,trial,indiv, 'VariableNames',{Yname,X1name,'trial','indiv'});


% Model in Wilkinson notation
% modelspec = [Yname ' ~ 1 + ' X1name '+' X2name ' + trial + (1|indiv)'];
modelspec = [Yname ' ~ 1 + ' X1name ' + trial + (1|indiv)'];

% Run generalized linear mixed-effects model
glme = fitglme(tbl,modelspec,'Distribution','normal');

disp(glme)





function ePlot(xVal,mVal,sVal,mClr,lClr)

eAlpha = 0.2;
lAlpha = 1;

for i = 1:length(xVal)
    h(i) = line(xVal(i).*[1 1],[mVal(i)-sVal(i) mVal(i)+sVal(i)],...
            'Color',[lClr eAlpha],'LineWidth',3);
    %shold on
    hold on
end
% h = errorbar(xVal,mVal,sVal,sVal,'Color',lClr);
% set(h,'LineStyle','none')
hold on
h = scatter(xVal,mVal,65,mClr,'MarkerFaceColor',mClr,...
    'MarkerEdgeColor','w');
hold off

c = polyfit(xVal,mVal,1);

line([min(xVal) max(xVal)],polyval(c,[min(xVal) max(xVal)]),...
    'Color',[lClr lAlpha],'LineWidth',2);

% Set x-axis
% xL = xlim;
% xlim([xL(1)-0.1*range(xL) xL(2)+0.1*range(xL)])

set(gca,'TickDir','out')






function [xRot,yRot] = rotCoords(xRaw,yRaw)

% Extract coordinates
x = xRaw - xRaw(1);
y = yRaw - yRaw(1);

% Fit trajectory
c = polyfit(x,y,1);

% Adjust coords
y = y - c(2);

% Refit trajectory
c = polyfit(x,y,1);

% x-axis point
xAx = [800 polyval(c,800)];

% Transform
coords = transCoord2d('xax G2L',[x(1) y(1)],xAx,[x y]);

% Record coords
xRot = coords(:,1);
yRot = coords(:,2);

if 0
    plot(x,y,'.')
end



function [tVal,yVal,xS,yS] = findPeaker(t,x,y)
% Finds valleys in the data with a variable cutoff frequency using a
% butterworth filter

% Whether to step trhu data plotting
makePlots = 0;

% Starting cutoff frequency, order of filter, sample rate
cutFreq       = 0.47;
cutStep      = 0.02;
filtOrd      = 4;
sampleRate   = 1/nanmean(diff(t));

% Threshold diference in time to stop iterations
tThresh = 1/sampleRate/1000;

% Low-pass filter the data, find peaks
[b, a] = butter(filtOrd, cutFreq./sampleRate,'low');
yS     = filtfilt(b, a, y);
[yEvents0,tEvents0] = findpeaks(yS,t);

% Mean period
meanT = mean(diff(tEvents0));

if makePlots 
    f = figure; 
end

% Loop thru cut-off frequencies
tEvents0_last = tEvents0;
nChange = 1;
stopLoop = 0;

while true

    % Advance cut-off frequency
    cutFreq = (cutFreq + cutStep);
    
    % Normalize cutoff
    Ws = cutFreq / sampleRate;
    if Ws >= 0.5
        break
    end
    
    % Low-pass filter the data, find peaks
    [b, a]   = butter(filtOrd, Ws,'low');
    yS       = filtfilt(b, a, y);
    xS       = filtfilt(b, a, x);
    [yE,tE]  = findpeaks(yS,t);
    
    % Store copy of tEvents0
    t_tEvents0 = tEvents0;
    t_yEvents0 = yEvents0;
    
    % Step trhu peaks
    for j = 1:length(tEvents0)

        % Index of peak in interval
        dVal  = abs(tEvents0(j)-tE);
        iPeak = find(dVal==min(dVal),1,'first');
        
        % Stop everything, if no values
        if isempty(iPeak)
            tEvents0 = t_tEvents0;
            yEvents0 = t_yEvents0;
            stopLoop = 1;
            break
        end
        
        % Store new peak, if there is only one
        yEvents0(j) = yE(iPeak);
        tEvents0(j) = tE(iPeak);
    end
    
    % Break the loop, or store new tEvents0
    if stopLoop
        break
    end
    
%     tEvents0 = t_tEvents0;
    
    % Mean change from last iteration
    dChange = mean(abs(tEvents0-tEvents0_last));
    
    % If time hasen't changed much, add to counter
    if dChange<tThresh
        nChange = nChange + 1;
    end
      
    % Break loop, if theres a lot of no changes
    if nChange>40
        break
    end
    
    % Store last eventtimes
    tEvents0_last = tEvents0;
    
    if makePlots
        figure(f)
        subplot(3,1,1)
        plot(t,y,'k-',t,yS,'b-',tEvents0,yEvents0,'ro')
        title(['cutFreq = ' num2str(cutFreq) ' Hz'])
        pause(0.1)
    end
end

% Store results
tVal  = tEvents0;
yVal  = yEvents0;






