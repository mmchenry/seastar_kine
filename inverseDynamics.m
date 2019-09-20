function inverseDynamics
% Perform an inverse dynamic analysis based on kinematics




%% Manage paths (need to modify for new PC)
   
if ~isempty(dir(['/Users/mmchenry/Documents/Matlab code']))

    % Path to root dir of video (CSULB project, external drive)
    vidPath = '/Users/mmchenry/Documents/Video/Chip sea stars/prelim video';
    
    % Path to root of data
    dataPath = '/Users/mmchenry/Documents/Projects/Chip sea stars/prelim data';

elseif isfolder('C:\Users\tpo\Documents\seastar_kine')
    
    % Path to root dir of video (CSULB project, external drive)
    vidPath = 'C:\Users\tpo\Documents\Video\Chip sea stars\prelim video';

    % Path to root of data
    dataPath = 'C:\Users\tpo\Documents\Chip sea star data\prelim data';
    
else
    
    error('Do not recognize computer')
    
end


%% Catalog sequence videos

%cList = catVidfiles(vidPath);

%cList.fName = 'S004_S001_T007';
if nargin<1
    %cList.fName = 'SS001_S001_T013';
    cList.fName = 'S005_S001_T011';
else  
    cList.fName = fileName;
end

cList.ext   = 'MOV';
cList.movtype = 'mov';
cList.path = '';

% Paths for current sequence
currDataPath = [dataPath filesep cList.path filesep cList.fName];
currVidPath  = [vidPath filesep cList.path filesep cList.fName '.' cList.ext];

% Check video path 
if ~isfile(currVidPath)
    error(['Video file does not exist at ' currVidPath]);
end

% Load video info (v)
v = defineVidObject(currVidPath);


%% Load data for sequence

% Load initial conditions (iC)
load([currDataPath filesep 'Initial conditions'])

% Load body kinematics (Body)
load([currDataPath filesep 'Body, post.mat'])

% Load F structure
%load([currDataPath filesep 'post- foot data'])
load([currDataPath filesep 'post- foot refined'])


%% Sequence-specific information

if strcmp(cList.fName,'S005_S001_T011')
    
    % Start and end of bouncing gait (s)
    Body.tBounceStart = 17;
    Body.tBounceEnd   = 35;
    
    % Body mass (kg)
    Body.mass = 10.27e-3;
    
    % Calibration constant (m/pix)
    Body.cal = 5.3238e-05;
end


%% Select times

if ~isfile([currDataPath filesep 'inverseDynTimes.mat'])
    
    % Index for duration of bouncing gait
    iBounce = (Body.t >= Body.tBounceStart) & (Body.t <= Body.tBounceEnd);
    
    % Calculate speed
    spd = Body.cal.*hypot(diff(smooth(Body.xCntr(iBounce))),...
        diff(smooth(Body.yCntr(iBounce))))./diff(Body.t(iBounce));
    spd = [spd(1); spd];
    
    f = figure;
    plot(Body.t(iBounce),spd*1000,'-k')
    grid on
    xlabel('t (s)')
    ylabel('Body speed (mm/s)')
    
    % Get times of peak velocity
    disp(' '); disp('Pick off peak velocity values (return when done)')
    T.tPeak = getTimes;
    
    % Get start times before velocity
    disp(' '); disp('Pick off start velocity values (return when done)')
    T.tStart = getTimes;
    
    % Check times
    if length(T.tPeak)~=length(T.tStart)
        error('Need to select same number of start and peak times');
        
    elseif max(T.tPeak<T.tStart)
        error('Start times need to preceed peak time');
    end
    
    close(f)
    
    % Save data
    save([currDataPath filesep 'inverseDynTimes'],'T');
    
else
    % Load T
    load([currDataPath filesep 'inverseDynTimes'])
end


%% Survey foot data

% Minimum number of feet to consider for force calculation
minFeet = 3;

% Start index
j = 1;

% Loop thru feet, saving properties
for i = 1:length(F)
    
     % If there are coodrinates (excluding feet that start at first frame)
     if max(~isnan(F(i).armNum)) %&& F(i).frames(1)>F(1).frames(1)

         % Mean position along arm
         armPos(j,1) = mean(F(i).xA);
        
         % Arm number
         armNum(j,1)     = F(i).armNum(1);
         
         % Index of current foot
         ftIdx(j,1) = i;
         
         % Start frame
         tStart(j,1) = F(i).frames(1)./v.FrameRate;
         
         % End frame
         tEnd(j,1) = F(i).frames(end)./v.FrameRate;   
         
         % Bounds of angle
         minAng(j,1) = min(F(i).ftAng);
         maxAng(j,1) = max(F(i).ftAng);
         
         % Advance index
         j = j + 1;     
     end 
end

% Calc normalized foot angle
for i = 1:length(armNum)
    angNorm{i} = (F(ftIdx(i)).ftAng-minAng(i))./(maxAng(i)-minAng(i));
end


%% Force calculation


f1 = figure;
f2 = figure;

% Loop thru power strokes
for i = 1:length(T.tStart)
    
    % Index for duration of power stroke
    idx = find((Body.t>=T.tStart(i)) & (Body.t<=T.tPeak(i)));

    % Solve for force
    t      = Body.t(idx);
    pos    = (Body.xT(idx)-Body.xT(idx(1))).*Body.cal;
    pos    = smooth(pos,10);
%     sp     = spaps(t,pos,1e-13);
%     f = fit(t,pos,'smoothingspline');
%     [vel,accel] = differentiate(f,t);
    warning off
    c = polyfit(t,pos,5);
    warning on

    accel = polyval(polyder(c,2),t);


    Force  = accel .* Body.mass;
    
    % Loop thru time values
    for j = 1:length(t)
        
        % Indicies of foot data include current frame
        iFeet = t(j)>=tStart & t(j)<=tEnd;
        
         % Indicies from F structure
        iF = ftIdx(iFeet);  
    
        % If there are enough feet at this time
        if ~isempty(iF) && length(iF)>minFeet
            
            % Loop thru each matching foot
            for k = 1:length(iF)
                
                % Time values for current foot
                tVal = F(iF(k)).frames./v.FrameRate ;
                
                %             cFtAng(k,1) = 2*pi*interp1(tVal,cAngVal,cTime);
                
                % Normalized position values in power stroke
                cAngVal = (F(iF(k)).ftAng-min(F(iF(k)).ftAng))./...
                                  range(F(iF(k)).ftAng);
                
                % Find current angle for this foot
                cFtAng(k,1) = interp1(tVal,cAngVal,t(j));

            end
        else
            cFtAng = nan;
        end
        
        ftPos(j,1) = mean(cFtAng);
        ftN(j,1)   = length(cFtAng);
        
        clear cFtAng
    end
    
    ftPos = smooth(ftPos);
    ftN = smooth(ftN);
    
    figure(f1);
    
    subplot(5,1,1)
    %plot(t,pos.*1000,'.r',t,f(t).*1000,'-k')
    plot(t,pos.*1000,'.r',t,polyval(c,t).*1000,'-k')
    hold on
    
    subplot(5,1,2)
    plot(t,Force.*1000,'-k')
    hold on
    
    subplot(5,1,3)
    plot(t,ftPos,'-k')
    hold on
    
    subplot(5,1,4)
    plot(t,ftN,'-k')
    hold on
    
    subplot(5,1,5)
    plot(t,Force.*1000./ftN,'-k')
    hold on
    
    figure(f2)
    
    plot(ftPos,Force.*1000./ftN,'o-')
    hold on
    
    clear t pos sp Force idx ftPos ftN
end

figure(f1)

subplot(5,1,1)
xlabel('Time (s)')
ylabel('Position (mm)')
grid on

subplot(5,1,2)
xlabel('Time (s)')
ylabel('Total force (mN)')
grid on

subplot(5,1,3)
xlabel('Time (s)')
ylabel('Mean foot position')
grid on

subplot(5,1,4)
xlabel('Time (s)')
ylabel('Number of feet')
grid on

subplot(5,1,5)
xlabel('Time (s)')
ylabel('Force per foot (mN)')
grid on


figure(f2)
xlabel('Relative position')
ylabel('Force per foot (mN)')


ttt=3;

%TODO: Fit SPAPS to position data, find second derivative and calculate
%force
% Then, determine the state of each foot, then plot force as a function
% of average stroke phase. (be sure to divide by number of feet)
    
    



function tPeak = getTimes

idx = 1;
h = [];
tPeak = [];
yL = ylim;

while true

    % Input
   [x,y,b] = ginput(1); 
    
   % If point selected
   if b==1
       
       % Store time values
       tPeak(idx) = x;
       idx = idx + 1;
       
   % If deleting point    
   elseif b==3 & isfield(T,'tPeak')
       
       % Subtract 
       T.tPeak = T.tPeak(1:end-1);
       
       % Reduce idx
       idx = max([1 idx-1]);
    
   % If return
   elseif isempty(b)
       break
   end
   
   hold on
   
   delete(h);
   h = [];
   
   % Display
   for i = 1:length(tPeak)   
       h(i) = plot(tPeak(i).*[1 1],yL,'r-');
   end
   hold off
   
end

if isempty(tPeak) 
    error('You need to select at least one peak')
end


    