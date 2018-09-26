function makeDataMovie(vid_path,v,S,mov_path,imVis,movType,varargin)


imInvert = 0;
dSample = 0;
numroipts = 300;
mSize(1) = 500;

clrs{1} = [241 90 36]./255;
clrs{2} = [41 171 226]./255;
clrs{3} = [0 1 0];

% Make figure
f = figure;



%% Creat output


% Loop thru data
for i = 1:length(S.frames)
    
    % Current whole frame
    im = getFrame(vid_path,v,S.frames(i),imInvert,'rgb');
    
    
    if strcmp(movType,'two view')
        
        % Get roi image
        [im_roi,bw_mask] = giveROI('stabilized',im,S.roi(i),...
            dSample,S.tform(:,:,i));
        
        % Loop thru tube feet
        for j = 1:length(S.ft)
            
            % Global tip point
            tipG(j,:) = [S.ft(j).xTip(i) S.ft(j).yTip(i)];
                
        end
        
        % Local tip point
        tipL = transCoord2d('G2L',S.tform(i),tipG);
        
        % Plot image
        figure(f)
        
        % LEFT IMAGE  ---------------------
        subplot(1,2,1)
        h = imshow(im,'InitialMag','fit');
        title(['Frame ' num2str(S.frames(i))])
        
        hold on
        
        % Outline roi (global)
%         line(S.roi(i).xCntr,S.roi(i).yCntr,'Color',...
%             [clrs{2} 0.5],'LineWidth',2);
        line(S.roi(i).xPerimG,S.roi(i).yPerimG,'Color',...
            [clrs{2} 0.5],'LineWidth',1);
        
        
        hold off
        
        % RIGHT IMAGE  ---------------------
        subplot(1,2,2)
        h = imshow(im_roi,'InitialMag','fit');
        hold on
        scatter(tipL(:,1),tipL(:,2),'MarkerEdgeColor',clrs{2},...
            'LineWidth',1,'MarkerEdgeAlpha',0.5,'SizeData',mSize(1));
        hold off
        
        clear tipG tipL
    end
    rrr=4;
end


