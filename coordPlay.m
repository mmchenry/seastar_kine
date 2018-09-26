function coordPlay
% Toying with matlab's bult-in coordinate transformations


%% Define local coordinate system


% Origin
or = [1 2 0]';

% Axes tilted at angle theta (30 deg)
theta = 30/180*pi;

% Axes defined
xaxis = [cos(theta) sin(theta) 0]';
yaxis = [cos(theta+pi/2) sin(theta+pi/2) 0]';
zaxis = [0 0 1]';

S = [xaxis yaxis zaxis];


%% Define point in global frame

% phi - angle of point in local FOR
phi = 45/180*pi;

% Dista from local FOR
r = 3; 

% Point in global FOR
pG = [ (or(1) + r*cos(theta+phi))  (or(2) + r*sin(theta+phi))  0]';

% Local coordinate
pL = global2localcoord(pG,'rr',or,S);

% Check value of phi
phiCheck = acos(pL(1)/r)/pi*180

% Check r
rCheck = hypot(pL(1),pL(2))


%figure;
subplot(1,2,1)
plot(0,0,'+k')
hold on
plot([xaxis(1)+or(1) or(1)],[xaxis(2)+or(2) or(2)],'k-')
plot([yaxis(1)+or(1) or(1)],[yaxis(2)+or(2) or(2)],'k-')
plot(pG(1),pG(2),'ro',[pG(1) or(1)],[pG(2) or(2)],'r-')
title('Global')
hold off
axis equal

subplot(1,2,2)
plot(0,0,'+k')
hold on
plot(pL(1),pL(2),'ro',[pL(1) 0],[pL(2) 0],'r-')
title('Local')
hold off
axis equal
