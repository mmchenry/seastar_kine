function anaSide
% Runs analysis of the kinematics of the chip side view 

%% Execution control

% Import DLC data again
do.reImportData = 0;


%% Parameters

% Get root paths
paths = givePaths;

% Extension for viode file names
extVid = 'MOV';

% Frame rate 
%TODO: Read this from the spreadsheet
fps = 24;

% Get list of all side csv files
dlc_path = [paths.data filesep 'data' filesep 'side_view' filesep '*.csv'];
dlc_files = dir(dlc_path);

% Check
if isempty(dlc_files)
    error(['No DLC csv files found in ' dlc_path]);
end


%% Read data from catalog spreadsheet

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
    if T.use_video(i)==1 
       
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
        seq(j).calConst       = T.cal(i);
        
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


%% Import data

if do.reImportData || ~isfile([paths.data filesep 'SideDataPooled.mat'])
    
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
            error([sideFile ...
                ' has no match to the sequences in Weights_experiments.xlsx'])
        end
        
        % get info on current sequence
        cSeq = seq(iCurrent);
        
        % Compile data --------------------------------------------------------
        
        % Read current csv data
        T = readtable([paths.data filesep 'data' filesep 'side_view' filesep ...
                       dlc_files(j).name],'HeaderLines',2);
        
        % Add time vector
        T.t = [0:(1/fps):((length(T.x)-1)/fps)]';
        
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
        
        %TODO: read the calibration data
        
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
    
    % Save
    save([paths.data filesep 'SideDataPooled'],'S'); 
    
end


%% Event finding

if 1 %do.reImportData || ~isfile([paths.data filesep 'SideDataPooled.mat'])
    
    % Load 'S' structure
    load([paths.data filesep 'SideDataPooled.mat'])
    
    % Figure parameters
    makeEventPlots = 1;
    makeSeqPlots   = 0;
    nPanels        = 6;
    iPanel         = 1;
    
    % Acceptable deviation from bounce amplitude and period
    yThresh  = 0.8;
    tThresh  = 0.8;
    
    if makeEventPlots
        f = figure;
    end
    
    if makeSeqPlots
        f2 = figure;
    end
    
    % Loop trhu sequences
    for i = 1:length(S)
        
        % Find when the body hits the ground
        [tEvent,yEvent,xS,yS] = findPeaker(S(i).t,S(i).xRaw,S(i).yRaw);
        
        % Apply calibration constant
        xS       = xS      .* S(i).calConst;
        yS       = yS      .* S(i).calConst;
        yEvent   = yEvent  .* S(i).calConst;
        
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
                
                % Index for second half of power stroke
                %             idx2 = idx;
                %             idx2(1:round(length(idx2)/2)) = 0;
                
                yRange(j,1) = max(S(i).yRaw(idx)) - min(S(i).yRaw(idx));
                
                meanSpd(j,1) = nanmean(xSpd(idx));
            else
                yRange(j,1) = 0;
                meanSpd(j,1) = nan;
            end
            
        end
        
        % Identify non-bouncing events
        iCrawl = (yRange < yThresh*mean(yRange)) | ...
                 (sPeriod < tau*tThresh) | (sPeriod > (2-tThresh)*tau);
        
        % Store
        S(i).xS                 = xS;
        S(i).yS                 = yS;
        S(i).xSpd               = xSpd;
        S(i).crawlSpd           = nanmean(meanSpd(iCrawl));
        S(i).bounceSpd          = nanmean(meanSpd(~iCrawl));
        S(i).tLand              = tEvent;
        S(i).yLand              = yEvent;
        S(i).iCrawl             = iCrawl;
        S(i).stepPeriod         = sPeriod;
        S(i).propBounce         = sum(sPeriod(~iCrawl)) / range(S(i).t);
        S(i).meanBouncePeriod   = mean(sPeriod(~iCrawl));
        
        
        % Plots
        if makeSeqPlots
            figure(f2)
            subplot(3,1,1)
            plot(S(i).t,S(i).yRaw .* S(i).calConst); hold on
            plot(S(i).t,yS,'k-');
            plot(S(i).tLand(~iCrawl),S(i).yLand(~iCrawl),'ro')
            plot(S(i).tLand(iCrawl),S(i).yLand(iCrawl),'r+')
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
        end
        
        
        if makeEventPlots
            figure(f)
            subplot(nPanels,1,iPanel)
            plot(S(i).t,S(i).yRaw .* S(i).calConst); hold on
            plot(S(i).t,yS,'k-');
            plot(S(i).tLand(~iCrawl),S(i).yLand(~iCrawl),'ro')
            plot(S(i).tLand(iCrawl),S(i).yLand(iCrawl),'r+')
            hold off
            xlabel('t (s)'); ylabel('Y')
            grid on
            title(['ExpType = ' S(i).expType ', ' S(i).fName_side])
            
            iPanel = iPanel + 1;
            
            if iPanel>nPanels && i<length(S)
                iPanel = 1;
                f = figure;
            end
        end
        
        clear tEvent yEvent xS yS iCrawl idx yRange meanSpd xSpd sPeriod
    end
    % Save
    save([paths.data filesep 'SideDataPooled'],'S');
    
else
    
    % Load 'S' structure
    load([paths.data filesep 'SideDataPooled.mat'])
    
end


%% Summary box plots


bouncePeriod = []; 
propBounce   = []; 
expType      = [];
meanSpd      = [];
bounceSpd    = [];
crawlSpd     = [];

for i = 1:length(S)
    
    if S(i).expType=='w' && S(i).addMass==2.09
        expType = [expType; 1];
        
    elseif S(i).expType=='w' && S(i).addMass==1.04
        expType = [expType; 2];
    
    elseif S(i).expType=='c' 
        expType = [expType; 3];
        
    elseif S(i).expType=='f' && S(i).floatNum==1
        expType = [expType; 4];
        
    elseif S(i).expType=='f' && S(i).floatNum==2
        expType = [expType; 5];
        
    else
        error('Do not recognize experiment type');
    end
        

    meanSpd         = [meanSpd;       mean(S(i).xSpd)];
    bounceSpd       = [bounceSpd;     S(i).bounceSpd];
    crawlSpd        = [crawlSpd;      S(i).crawlSpd];
    
    bouncePeriod    = [bouncePeriod;  S(i).meanBouncePeriod];
    propBounce      = [propBounce;    S(i).propBounce];

end


figure;

subplot(1,3,1)
boxplot(meanSpd,expType)
ylabel('Mean speed (cm/s)')
set(gca,'TickDir','out')
ylim([0 3e-3])

subplot(1,3,2)
boxplot(bounceSpd,expType)
ylabel('Bounce speed (cm/s)')
set(gca,'TickDir','out')
ylim([0 3e-3])

subplot(1,3,3)
boxplot(crawlSpd,expType)
ylabel('Crawl speed (cm/s)')
set(gca,'TickDir','out')
ylim([0 3e-3])


figure;
subplot(1,2,2)
boxplot(bouncePeriod,expType)
ylabel('Bounce Period')
set(gca,'TickDir','out')

subplot(1,2,1)
boxplot(propBounce,expType)
ylabel('Proportion of time bouncing')
set(gca,'TickDir','out')



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






