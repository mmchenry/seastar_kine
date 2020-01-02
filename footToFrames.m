function C = footToFrames(F,Body)
% Reorganized foot data wrt video frames

% Index for frames
iFrames = Body.frames>=min(F(1).frames) & Body.frames<=max(F(end).frames);

% Frame numebrs to analyze
C.frames   = Body.frames(iFrames)';

for i = 1:length(C.frames)
    
    % Number of feet in frame
    n = 0;
    
    % Index for body data
    iMatchBody = Body.frames==C.frames(i);
    
    % Arm coordinates
    C.xArm(i,:)   = Body.xArmG(iMatchBody,:);
    C.yArm(i,:)   = Body.yArmG(iMatchBody,:);
    
    % Body center
    C.xCntr(i,:)   = Body.xCntr(iMatchBody);
    C.yCntr(i,:)   = Body.yCntr(iMatchBody);
    
    % Loop thru feet
    for j = 1:length(F)
        
        % Index of current frame to F data
        iMatch = F(j).frames==C.frames(i);
        
        % Check for repeats
        if sum(iMatch)>1
            error('Multiple instances of same frame for a foot');
        end
        
        % If there is a match, store coordinates
        if max(iMatch)
            % Advance index
            n = n + 1;
            
            % Store coordinates and colors
            C.x{i}(n,1)      = F(j).xG(iMatch);
            C.y{i}(n,1)      = F(j).yG(iMatch);
            C.clr{i}(n,:)    = F(j).clr(1,:);
        end
    end
    
    % Add nans, if no matching feet
    if n==0
        x{i}(1)      = nan;
        y{i}(1)      = nan;
        clr{i}(1,:)  = F(j).clr(1,:);
    end
end

