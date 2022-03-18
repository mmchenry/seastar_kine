function visManual
%
% Visualize the manually-tracked data

%% Execution control

% Just plot everything
do.plotAll = 0;

% Collect and plot individual pwr strokes
do.poolIndividual = 1;

% Make representative figure
do.representative = 0;


%% Parameters and data loading

% Frame rate for foot data
fps = 23.97602398;

nSmooth = 15;

% Get root paths
paths = givePaths;

% Load fube foot data 'F'
load([paths.data filesep 'footData.mat'])

% Load side kinematics data (S)
% load([paths.data filesep 'SideDataPooled_eventsManual2.mat']);


%% Plot everything

if do.plotAll

% Loop thru sequences
for i = 1:length(F)

    % Find match in S
%     for j = 1:length(S)
%         
%         if strcmp(S(j).fName_side,F(i).seq.fName_side)
%             iCurr = j;
%             break
%         end
%     end

    figure
    fs = 1/mean(diff(F(i).tSide));
    zS = highpass(smooth(F(i).zSide,nSmooth),0.05,fs);
%     zS = bandpass(F(i).zSide,[0.3 3],fs));

    subplot(3,1,1)
    plot(F(i).tSide,zS,'k-')
    xlabel('t (s)'); ylabel('Z (cm)')

    % Loop thru feet, each with one pwr stroke
    for j = 1:length(F(i).ftL)
        
        % Extract current data
        t       = F(i).ftL(j).tBase;
        xBase   = smooth(F(i).ftL(j).xBase,nSmooth)   .* 100;
        yBase   = smooth(F(i).ftL(j).yBase,nSmooth)   .* 100;
        zBase   = smooth(F(i).ftL(j).zBase,nSmooth)   .* 100;
        theta   = smooth(F(i).ftL(j).theta,nSmooth)   .* 180/pi;
        pod_len = sqrt(xBase.^2 + yBase.^2 + zBase.^2);
%         pod_len = smooth(F(i).ftL(j).pod_len,nSmooth) .* 100;

        % Make time vector
%         t = F(i).ftL(j).frames ./ F(i).ftL(j).fps;

        subplot(3,1,2)
        plot(t,pod_len,'k-')
        xlabel('t (s)'); ylabel('pod. length (cm)')
        hold on

        subplot(3,1,3)
        plot(t,theta,'k-')
        xlabel('t (s)'); ylabel('\theta (deg)')
        hold on
    end

    xL = xlim;

    subplot(3,1,1)
    xlim(xL)

    % Plot peaks
    yL = ylim;
    for j = 1:length(F(i).tLand)
        hold on
        plot(F(i).tLand(j).*[1 1],yL,'b-')
    end
    hold off
    title(['Seq' num2str(i) ', ' F(i).seq.dirName])
    set(gca,'TickDir','out')

    subplot(3,1,2)
    yL = ylim;
    for j = 1:length(F(i).tLand)
        hold on
        plot(F(i).tLand(j).*[1 1],yL,'b-')
    end
    hold off
    xlim(xL)
    set(gca,'TickDir','out')

    subplot(3,1,3)
    yL = ylim;
    for j = 1:length(F(i).tLand)
        hold on
        plot(F(i).tLand(j).*[1 1],yL,'b-')
    end
    hold off
    xlim(xL)
    set(gca,'TickDir','out')

end

end %do.plotAll



%% Representative

if do.representative

% Sequence index
iSeq = 1;

% Index of feet to highlight
iFt = [6 22 31 60 70];

% Start time
tStart = 30;

% Make figures
f1 = figure;
f2 = figure;

% Filter side view 
fs = 1/mean(diff(F(iSeq).tSide));
zS = highpass(smooth(F(iSeq).zSide,nSmooth),0.05,fs);

% Plot side view
figure(f1)
subplot(3,1,1)
plot(F(iSeq).tSide-tStart,zS.*1000,'k-')
xlabel('t (s)'); ylabel('Z (mm)')

% Panel index
k = 1;

% Loop thru feet, each with one pwr stroke
for j = 1:length(F(iSeq).ftL)

    % Extract current data
    t       = F(iSeq).ftL(j).tBase - tStart;
    xBase   = smooth(F(iSeq).ftL(j).xBase,nSmooth)   .* 1000;
    yBase   = smooth(F(iSeq).ftL(j).yBase,nSmooth)   .* 1000;
    zBase   = smooth(F(iSeq).ftL(j).zBase,nSmooth)   .* 1000;
    theta   = smooth(F(iSeq).ftL(j).theta,nSmooth)   .* 180/pi;
    pod_len = sqrt(xBase.^2 + yBase.^2 + zBase.^2);

    figure(f1)
    subplot(3,1,2)
    h1 = plot(t,pod_len);
    xlabel('t (s)'); ylabel('pod. length (mm)')
    hold on

    subplot(3,1,3)
    h2 = plot(t,theta);
    xlabel('t (s)'); ylabel('\theta (deg)')
    hold on
    
    % Highlight iFt
    if max(j==iFt)
        set(h1,'Color',[0 0 0],'LineWidth',2);
        set(h2,'Color',[0 0 0],'LineWidth',2);
    else
        set(h1,'Color',0.5.*[1 1 1],'LineWidth',0.5);
        set(h2,'Color',0.5.*[1 1 1],'LineWidth',0.5);
    end

    if max(j==iFt)
        lClr = 0.5.*[1 1 1];
        figure(f2)

        idx = 1:4:length(xBase);

        subplot(3,2,k)
        for j = 1:length(idx)
            h1 = line([xBase(idx(j)) 0],[zBase(idx(j)) 0],'Color',lClr);
            hold on
%             h2 = scatter(xBase(j),zBase(j),7,'filled',...
%                 'MarkerFaceColor',lClr,'MarkerEdgeColor',lClr);
%             hold on
            axis square
            xlim([-6 6])
            ylim([-2 10])
        end

        k = k + 1;
    end
end

disp(['Time interval (s) = ' num2str(mean(diff(t(1:4:end))))])

figure(f1)
xL = [0 30];
subplot(3,1,1)
xlim(xL)

% Plot peaks
figure(f1)
yL = ylim;
for j = 1:length(F(iSeq).tLand)
    hold on
    plot((F(iSeq).tLand(j)-tStart).*[1 1],yL,'b-')
end
hold off
set(gca,'TickDir','out','XLim',[0 30])

subplot(3,1,2)
yL = ylim;
for j = 1:length(F(iSeq).tLand)
    hold on
    plot((F(iSeq).tLand(j)-tStart).*[1 1],yL,'b-')
end
hold off
xlim(xL)
set(gca,'YLim',[0 8],'YTick',[0:2:8])
set(gca,'TickDir','out','XLim',[0 30])

subplot(3,1,3)
yL = ylim;
for j = 1:length(F(iSeq).tLand)
    hold on
    plot((F(iSeq).tLand(j)-tStart).*[1 1],yL,'b-')
end
hold off
xlim(xL)
set(gca,'YLim',[0 180],'YTick',[0:45:180])
set(gca,'TickDir','out','XLim',[0 30])



end %do.representative


%% Extract individual sequences

if do.poolIndividual


showTraces = 0;

if showTraces
    f1 = figure;
end

f2 = figure;
f3 = figure;

k = 1;

% Loop thru sequences
for i = 1:length(F)

    % Loop thru feet, each with one pwr stroke
    for j = 1:length(F(i).ftL)
        
        % Extract current data
        t{j}       = F(i).ftL(j).tBase;
        xBase{j}   = smooth(F(i).ftL(j).xBase,nSmooth)   .* 100;
        yBase{j}   = smooth(F(i).ftL(j).yBase,nSmooth)   .* 100;
        zBase{j}   = smooth(F(i).ftL(j).zBase,nSmooth)   .* 100;
        theta{j}   = smooth(F(i).ftL(j).theta,nSmooth)   .* 180/pi;
        pod_len{j} = sqrt(xBase{j}.^2 + yBase{j}.^2 + zBase{j}.^2);

        % Index for the bounds of 
        iPwrStart = find(F(i).tLand<min(t{j}),1,'last');
        iPwrEnd   = find(F(i).tLand>max(t{j}),1,'first');

        % If kinematics span 2 oscillations, select middle peak
        if iPwrEnd-iPwrStart==2
            iPwr(j) = iPwrStart + 1;

        % Otherwise, choose closest one
        elseif iPwrEnd-iPwrStart==1
            
            if min(t{j})-F(i).tLand(iPwrStart) < F(i).tLand(iPwrEnd)-max(t{j}) 
                iPwr(j) = iPwrStart;
            else
                iPwr(j) = iPwrEnd;
            end

        else
            iPwr(j) = 0;
        end

        if iPwr(j) == 0
            tPwr(j) = 0;
            TPwr(j) = 0;
        else
            tPwr(j) = F(i).tLand(iPwr(j));
            TPwr(j) = (F(i).tLand(iPwr(j)-1) - F(i).tLand(iPwr(j)+1)) / 2;
        end
    end

    % Get panel number
    if F(i).seq.expType=='f'
        nPlt = 1;
    elseif F(i).seq.expType=='c'
        nPlt = 2;
    elseif F(i).seq.expType=='w'
        nPlt = 3;
    else
        error('huh?')
    end

    % Index for sequence
    if i==7
        iSeq = iPwr(2);
    else
        iSeq = max(iPwr);
    end
    
    % Extract body motion
    fs = 1/mean(diff(F(i).tSide));
    zS = highpass(smooth(F(i).zSide,nSmooth),0.05,fs);
    tS = (F(i).tSide - max(tPwr))./mean(TPwr);

    % Plot body motion
    if showTraces
        figure(f1)
        subplot(3,3,nPlt)
        plot(tS,(zS-mean(zS)).*100,'k')
        hold on
        set(gca,'TickDir','out','YLim',[-0.15 0.15],'YTick',[-0.15:0.05:0.15])
    end

    % Loop thru power strokes
    for j = 1:length(t)
    
        % Plot just iSeq
        if iPwr(j)==iSeq
            if showTraces
                figure(f1)
                subplot(3,3,nPlt+3)
                plot((t{j}-tPwr(j))./TPwr(j),pod_len{j},'k')
                set(gca,'TickDir','out','YLim',[0 1],'YTick',[0:0.2:1])
                hold on

                subplot(3,3,nPlt+6)
                plot((t{j}-tPwr(j))./TPwr(j),theta{j},'k')
                set(gca,'TickDir','out','YLim',[0 180],'YTick',[0:45:180])
                hold on
            end

            iMin = find(pod_len{j}==min(pod_len{j}),1,'first');

            D.lStart(k,1) = pod_len{j}(1);
            D.lMin(k,1)   = min(pod_len{j});
            D.lEnd(k,1)   = pod_len{j}(end);
            D.tStart(k,1) = (t{j}(1)-tPwr(j))     ./TPwr(j);
            D.tMin(k,1)   = (t{j}(iMin)-tPwr(j))  ./TPwr(j);
            D.tEnd(k,1)   = (t{j}(end)-tPwr(j))   ./TPwr(j);

            D.thMin(k,1)  = min(theta{j});
            D.thMax(k,1)  = max(theta{j});
            D.grp(k,1)    = F(i).seq.expType;

            k = k + 1;
        end
        end

    clear t *Base theta pod_len iPwr tPwr TPwr
end

% X-limits
if showTraces
    figure(f1)
    for i = 1:9
        subplot(3,3,i)
        xlim([-1 1])
        if i<4
            ylim([-0.15 0.15])
        end
    end

    % Titles
    figure(f1)
    subplot(3,3,1)
    title('floats')
    subplot(3,3,2)
    title('control')
    subplot(3,3,3)
    title('weights')
    set(f1,'Renderer','Painters')
end

% Boxplots of length
figure(f2)
subplot(2,3,1)
boxplot(D.lStart,D.grp,'Symbol','.k')
% axis square
set(gca,'TickDir','out','YLim',[0 1],'YTick',[0:.2:1])
ylabel('lStart')

subplot(2,3,2)
boxplot(D.lMin,D.grp,'Symbol','.k')
% axis square
set(gca,'TickDir','out','YLim',[0 1],'YTick',[0:.2:1])
ylabel('lMin')

subplot(2,3,3)
boxplot(D.lEnd,D.grp,'Symbol','.k')
% axis square
set(gca,'TickDir','out','YLim',[0 1],'YTick',[0:.2:1])
ylabel('lEnd')

subplot(2,3,4)
boxplot(D.tStart,D.grp,'Symbol','.k')
% axis square
% set(gca,'TickDir','out','YLim',[0 1],'YTick',[0:.2:1])
ylabel('tStart')
set(gca,'TickDir','out','YLim',[-1 1],'YTick',[-1:0.5:1])

subplot(2,3,5)
boxplot(D.tMin,D.grp,'Symbol','.k')
ylabel('tMin')
set(gca,'TickDir','out','YLim',[-1 1],'YTick',[-1:0.5:1])

subplot(2,3,6)
boxplot(D.tEnd,D.grp,'Symbol','.k')
ylabel('tEnd')
set(gca,'TickDir','out','YLim',[-1 1],'YTick',[-1:0.5:1])


% Boxplots of theta
figure(f3)
subplot(2,3,1)
boxplot(D.thMin,D.grp,'Symbol','.k')
% axis square
% set(gca,)
set(gca,'TickDir','out','YLim',[0 180],'YTick',[0:45:180])

subplot(2,3,2)
boxplot(D.thMax,D.grp,'Symbol','.k')
% axis square
set(gca,'TickDir','out','YLim',[0 180],'YTick',[0:45:180])
set(f2,'Renderer','Painters')
set(f3,'Renderer','Painters')

end %do.poolIndividual

