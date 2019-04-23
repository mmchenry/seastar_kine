function [a,rTick] = circ_plot_mjm(alpha,rTick,symbolSize,maxN)
% MJM's modification of the circular plot below
%

% Grid color
gClr = 0.6.*[1 1 1];

% Marker color
mClr = 0.2.*[1 1 1];

% Mean value color
muClr = gClr;

% Distribution color
vmClr = 0.5.*[1 1 1];

% Number of bins (for symbols) 
numbins = 180;

% Tick length
ticklen = 0.3;

% Inner and outer circle radii (graph units)
inR = 1;
outR = 3;

% If no tick range defined
if nargin < 2
    rTick = outR - inR;
end


if nargin < 3
    symbolSize = 150;
end


% Position of text around rings
textPos = -3*pi/4;

% Font size (for plot)
fSize = 12;

% Concentration parameter
kappa = circ_kappa(alpha);

% Circular mean
mu    = circ_mean(alpha);

% Bin intervals
binedge = linspace(-pi,pi,numbins);

% Angular values to evaluate angles
phi = linspace(-pi,pi,1000);

% Bin data
[N,tmp] = histcounts(alpha,binedge);

if nargin<4 
    maxN = max(N);
end


% Center of bins
cntrAng = binedge(1:end-1)+mean(diff(binedge))/2;

% Symbol size
symSpace = 0.9*rTick/maxN;
symSize = symSpace*symbolSize;


%% Preliminary plot

% Circles
h1 = line(inR.*cos(phi),inR.*sin(phi),'Color',gClr);
axis square
xlim([-1 1].*(outR+2*ticklen))
ylim([-1 1].*(outR+2*ticklen))
hold on
delete(h1)


%% Plot central line

rMu = [outR outR+ticklen*2];
xMu = rMu.*cos(mu);
yMu = rMu.*sin(mu);

line(xMu,yMu,'Color',muClr,'LineWidth',2)

if mu<=pi/4 && mu>-pi/4
    xT = xMu(2) + 2*ticklen;
    yT = yMu(2);
    algn = 'left';
    
elseif (mu > pi/4 && mu<=3*pi/4) || (mu < -pi/4 && mu>=-3*pi/4)
    algn = 'center';
    xT = xMu(2);
    yT = yMu(2)+ 2*ticklen;
    
else
    algn = 'right';
    xT = xMu(2) + 2*ticklen;
    yT = yMu(2);
end
    
% Text label
h = text(xT,yT,num2str(mu.*180/pi,'%.1f'),...
    'Color',muClr,'HorizontalAlignment',algn,...
    'fontSize',fSize);


%% Plot VM distribution

% PDF values for VM distribution
[pdfVM phiVM] = circ_vmpdf(phi, mu, kappa);

% Test for uniformity
[pval, v] = circ_vtest(alpha, mu);

% If significantly non-uniform
if pval<0.05
    %pdfVM = pdfVM/rTick.*(outR-inR) + inR;
    pdfVM2 = pdfVM/max(pdfVM).*(outR-inR) + inR;
    
    x = pdfVM2.*cos(phiVM);
    y = pdfVM2.*sin(phiVM);
    
    line(x,y,'Color',vmClr,'LineWidth',1)
end


%% Stat test for pure pursuit

% stat test
[h, mu2,ci1,ci2] = circ_mtest(alpha, 0);

if h==0
    warning('Distribution not significantly different from 0')
else
    disp(['Mean direction: ' num2str(mu2*180/pi) ' CIs: ' ...
        num2str(ci1*180/pi) ', ' num2str(ci2*180/pi)])
end


%% Plot observation circles

% Loop thru center angle position
for i = 1:length(cntrAng)
    if N(i)~=0
        rC = inR + [symSpace:symSpace:(N(i))*symSpace];
        xC = rC.*cos(cntrAng(i));
        yC = rC.*sin(cntrAng(i));
        
        % Points for current slice
        scatter(xC,yC,symSize,mClr,'MarkerEdgeColor','none',...
            'MarkerFaceColor',mClr)
    end
end


%% Plot gridlines

% Circles
line(inR.*cos(phi),inR.*sin(phi),'Color',gClr);
% axis square
% xlim([-1 1].*(outR+ticklen))
% ylim([-1 1].*(outR+ticklen))
% hold on
line(outR.*cos(phi),outR.*sin(phi),'Color',gClr);

% Clear graphics
box off
set(gca,'XTick',[],'YTick',[],'YColor',[1 1 1],'XColor',[1 1 1])

% Circle labels
text(inR*cos(textPos),inR*sin(textPos),'0','Color',gClr,...
    'BackgroundColor',[1 1 1],'HorizontalAlignment','center','fontSize',fSize);
text(outR*cos(textPos),outR*sin(textPos),num2str(rTick,'%.1f'),'Color',gClr,...
    'BackgroundColor',[1 1 1],'HorizontalAlignment','center','fontSize',fSize);

% Angular tick: 0 deg
line([outR outR+ticklen],[0 0],'Color',gClr);
text(outR+ticklen*1.5,0,'0','Color',gClr,...
    'HorizontalAlignment','left','fontSize',fSize);

% Angular tick: 90 deg
line([0 0],[outR outR+ticklen],'Color',gClr);
text(0,outR+ticklen*2.5,'90','Color',gClr,...
    'HorizontalAlignment','center','fontSize',fSize);

% Angular tick: 180 deg
line(-1*[outR outR+ticklen],[0 0],'Color',gClr);
text(-outR-ticklen*1.5,0,'180','Color',gClr,...
    'HorizontalAlignment','right','fontSize',fSize);

% Angular tick: -90 deg
line([0 0],-1*[outR outR+ticklen],'Color',gClr);
text(0,-outR-ticklen*2.5,'-90','Color',gClr,...
    'HorizontalAlignment','center','fontSize',fSize);


% h = plot(rTick.*cos(phi),rTick.*sin(phi),'k-');
% hold on
% 
% line([0 1.2.*rTick.*cos(mu)],[0 1.2.*rTick.*sin(mu)],'Color',.9.*[1 1 1]);
% text(1.25.*rTick.*cos(mu),1.25.*rTick.*sin(mu),num2str(mu.*180/pi,'%.2f'))
% 
% set(h,'Color',gClr)
% xlim(1.1.*[-rTick rTick])
% ylim(1.1.*[-rTick rTick])
% axis square
% 
% text(1.1*rTick,0,num2str(rTick))
% line([0 0],1.1.*[-rTick rTick],'Color',gClr)
% line(1.1.*[-rTick rTick],[0 0],'Color',gClr)
% h = scatter(rTick.*cos(alpha),rTick.*sin(alpha),mSize,0.5.*[1 1 1]);
% set(h,'MarkerEdgeColor','none','MarkerFaceColor',mClr)
% 
% h = fill(x,y,mClr);
% set(h,'EdgeColor','none')



% % Von Mises distribution
% p = polar(phiVM,pdfVM);
% h = findall(gca,'type','line');
% % remove the handle for the polar plot line from the array
% h(h == p) = [];
% % delete all other lines
% delete(h);
% %delete(findall(gcf,'type','text'))
% 
% hold on
% 
% % Get radius of plot
% %rTick = max(xlim);
% 
% % Plot measurements
% p = polar(alpha,rTick.*ones(length(alpha),1),'o');
%     set(p,'MarkerFaceColor',0.5.*[1 1 1],...
%         'MarkerEdgeColor','none','MarkerSize',8);
%        
    




% If you want, use the following lines to remove the text.
% (I prefer to leave it in place)
% find and remove all of the text objects in the polar plot
% delete(findall(gcf,'type','text'))    
    
hold off


% r = circ_r(alpha) * rTick;
% phi = circ_mean(alpha);
% hold on;
% zm = r*exp(1i*phi);
% plot([0 real(zm)], [0, imag(zm)],'-k')
% hold off;

a = gca;




