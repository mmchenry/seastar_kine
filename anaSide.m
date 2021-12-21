function anaSide
% Runs analysis of the kinematics of the chip side view 

%% Execution control

% Import DLC data again
do.reImportData = 0;

% Analyze audio data to sync the timing of videos
do.anaAudioSync = 0;

% Visualize events for all sequences
do.visEvents = 0;

% Visualize details of each sequence
do.visSeqs = 0;

% Plot summary boxpolot data
do.summaryPlots = 0;

% Plot data that keeps track of individuals
do.indivPlots = 1;

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
            seq(j).fps            = T.frame_rate_side(i);
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
            T.t = [0:(1/seq(i).fps):((length(T.x)-1)/seq(i).fps)]';

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
                S(iSeq).fName_side = cSeq.fName_side;
                S(iSeq).fName_bot  = cSeq.fName_bot;
                S(iSeq).expType    = cSeq.expType;
                S(iSeq).addMass    = cSeq.addMass;
                S(iSeq).floatNum   = cSeq.floatNum;
                S(iSeq).indiv      = cSeq.indiv;
                S(iSeq).bodyMass   = cSeq.bodyMass;
                S(iSeq).calConst   = cSeq.calConst;
                S(iSeq).SW_percent = cSeq.SW_percent;
                S(iSeq).fps        = cSeq.fps;


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

if 1 %do.reImportData || ~isfile([paths.data filesep 'SideDataPooled_events.mat'])
    
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
        S(i).crawlSpd           = nanmean(meanSpd(iCrawl));
        S(i).bounceSpd          = nanmean(meanSpd(~iCrawl));
        S(i).tLand              = tEvent;
        S(i).yLand              = yEvent;
        S(i).yBounceAmp         = mean(yRange(~iCrawl)); 
        S(i).iCrawl             = iCrawl;
        S(i).stepPeriod         = sPeriod;
        S(i).propBounce         = sum(sPeriod(~iCrawl)) / range(S(i).t);
        S(i).meanBouncePeriod   = mean(sPeriod(~iCrawl));
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
    
    % Save
    save([paths.data filesep 'SideDataPooled_events'],'S');   
end


%% Visualize event finding

if do.visEvents
    
    makeEventPlots = 1;
    
    if ~exist('S','var')
        % Load 'S' structure
        load([paths.data filesep 'SideDataPooled_events.mat'])
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
end


%% Visualize details on each sequence

if do.visSeqs
    
    if ~exist('S','var')
        % Load 'S' structure
        load([paths.data filesep 'SideDataPooled_events.mat'])
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
    
% Load 'S' structure
load([paths.data filesep 'SideDataPooled_events.mat'])

% Loop thru sequences, collect data
for i = 1:length(S)

    indiv(i,1)          = S(i).indiv;
    meanBPeriod(i,1)    = S(i).meanBouncePeriod; 
    SW_percent(i,1)     = S(i).SW_percent;
    beat_spd(i,1)       = S(i).bounceSpd;
    tFirstB(i,1)        = S(i).tFirstBounce;
    yAmpB(i,1)          = S(i).yBounceAmp;
end

indNums = unique(indiv);


f = figure;

lClr = lines(length(indNums));

for i = 1:length(indNums)
    
    % Placeholder for data
%     nanPlace       = nan(sum(indiv==indNums(i)),1);
%     BPeriod_mean   = nanPlace;
%     BPeriod_std    = nanPlace;

    % Index of values for current individual
    iIdx = indiv==indNums(i);

    % Unique SW values
    SWs = unique(SW_percent(iIdx));

    
    % Loop thru SW values
    for j = 1:length(SWs)

         % Index for current values
         idx = iIdx & SW_percent==SWs(j);
         
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
    subplot(4,2,1)
    ePlot(SWs,BPeriod_mean,BPeriod_std,lClr(i,:),lClr(i,:))
    ylabel('Bounce period')
    xlabel('Percent submerged weight')
    axis square
    hold on
    xlim([40 175])

    % Bounce speed
    subplot(4,2,2)
    ePlot(SWs,Bspd_mean,Bspd_std,lClr(i,:),lClr(i,:))
    ylabel('Bounce speed (m/s)')
    xlabel('Percent submerged weight')
    axis square
    hold on
    xlim([40 175])

    % Time to first bounce
    subplot(4,2,3)
    ePlot(SWs,BtFirst_mean,BtFirst_std,lClr(i,:),lClr(i,:))
    ylabel('Time to first bounce (s)')
    xlabel('Percent submerged weight')
    axis square
    hold on
    xlim([40 175])

    % Time to first bounce
    subplot(4,2,4)
    ePlot(SWs,Bamp_mean,Bamp_std,lClr(i,:),lClr(i,:))
    ylabel('Y Amplitude')
    xlabel('Percent submerged weight')
    axis square
    hold on
    xlim([40 175])
    
    clear SWs B*
end


end %do.indivPlots


%% Plots of summary data

if do.summaryPlots
    
    if ~exist('S','var')
        % Load 'S' structure
        load([paths.data filesep 'SideDataPooled_events.mat'])
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


function ePlot(xVal,mVal,sVal,mClr,lClr)

errorbar(xVal,mVal,sVal,sVal,'Color',lClr)
hold on
h = scatter(xVal,mVal,50,mClr,'MarkerFaceColor',mClr,...
    'MarkerEdgeColor',lClr);
hold off

% Set x-axis
xL = xlim;
xlim([xL(1)-0.1*range(xL) xL(2)+0.1*range(xL)])

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






