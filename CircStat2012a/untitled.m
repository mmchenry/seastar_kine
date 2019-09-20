



idx = index{i} & ~isnan(b.preyx2(:,2)) & ...
          spd1>spd2;
      
      
[m,l1,l2] = circ_mean(dir(idx));

m * (180/pi)
l1 * (180/pi)
l2 * (180/pi)       


sum(idx)



sum(index{i} & ~isnan(b.preyx2(:,2)))

sum(index{i} & ~isnan(b.preyx2(:,2)) & ~wrong)




%%

Dconv = 21.944; %px/mm
figure;

cd('/Volumes/RED/Ch4 Data');

[FILENAME, PATHNAME] = uigetfile('*.*');
    
cd(PATHNAME);

j = str2num(PATHNAME(end-3:end-1));


while 1
   
    who = dir(PATHNAME);
    
    FILENAME = who(10).name;
    
    i = imread([PATHNAME filesep FILENAME]);
    
  
    imshow(i);
    
    title(PATHNAME);
    
    [a,b] = ginput(1);
    
    zoomcenter(a,b);
    
    [x,y] = ginput(2);
    
    dx = diff(x);
    dy = diff(y);
    
    bl = (dx^2 + dy^2)^0.5;
    
    BL(j) = bl / 21.944;  %mm
    
    
    j = j+1;
    
    next = sprintf('%03d',j);
    
    PATHNAME(end-3:end-1) = next;
    
    cd(PATHNAME);
    
    
    
    
    
    
    
    
    clf;
    
    
    
    
    
    
end





dx = b.preyx(:,1) - b.preyx(:,3)
dy = b.preyy(:,1) - b.preyy(:,3)
bd = [dx dy]
bdp = cart2pol(bd(:,1),bd(:,2))


for j = 1: length(who)
    plot3(b.preyx(who(j),:),b.preyy(who(j),:),b.preyz(who(j),:));
    hold on;
    plot3(b.preyx(who(j),1),b.preyy(who(j),1),b.preyz(who(j),1),'o');
    
    
end

