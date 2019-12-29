function paths = givePaths
% Returns a structure of paths for running code


% Matt's computer
if ~isempty(dir(['/Users/mmchenry/Documents/Matlab code']))

    % Path to root dir of video (CSULB project, external drive)
%     vidPath = '/Users/mmchenry/Documents/Video/Chip sea stars/prelim video';
    paths.vid = '/Users/mmchenry/Documents/Video/Sea stars/SICB2020';
    
    % Location of video frames
    %vidFramePath = '/Users/mmchenry/Documents/Video/Chip sea stars/prelim video/video frames';
    
    % Path to root of data
%     paths.data = '/Users/mmchenry/Documents/Projects/Chip sea stars/prelim data';
    paths.data = '/Users/mmchenry/Documents/Projects/Seastars/SICB2020';

elseif isfolder('C:\Users\tpo\Documents\seastar_kine')
    
    % Path to root dir of video (CSULB project, external drive)
    paths.vid = 'C:\Users\tpo\Documents\Video\Chip sea stars\prelim video';

    % Path to root of data
    paths.data = 'C:\Users\tpo\Documents\Chip sea star data\prelim data';

elseif isfolder('C:\Users\McHenryLab\Documents\GitHub\seastar_kine')
    
    % Path to root dir of video (CSULB project, external drive)
    paths.vid = 'C:\Users\McHenryLab\Documents\Seastar\SICB2020\floats\videos';

    % Path to root of data
    paths.data = 'C:\Users\McHenryLab\Documents\Seastar\SICB2020\floats\data';

else
    
    error('Do not recognize computer')
    
end




% % Matt's laptop
% if ~isempty(dir(['/Users/mmchenry/Documents/Matlab code']))
%     
%     % Use flow, if connected
%     if ~isempty(dir('/Volumes/shared'))
%         paths.vid = '/Volumes/shared';
%     else
%         paths.vid = '/Users/mmchenry/Documents/Video/Sea stars';
%     end
%     
%     % Use Google Drive, if connected
%     if ~isempty(dir('/Volumes/GoogleDrive/My Drive'))
%         paths.data = '/Volumes/GoogleDrive/My Drive/Shared files/Projects/Andres sea stars/Kinematics';
%     else
%         paths.data = '/Users/mmchenry/Documents/Projects/Andres sea stars/Kinematics';
%      end
%     
% %Amberle's Laptop
%     elseif ~isempty(dir(['C:\Users\biolo_000\Desktop\Sea Star Videos']))
%     
%     paths.vid = 'C:\Users\biolo_000\Desktop\Sea Star Videos';
%     paths.data = 'G:\My Drive\Kinematics';
%     
%     
% elseif ~isempty(dir(['C:\Users\tpo\Documents']))
%     
%     paths.vid = 'C:\Users\tpo\Documents\Andres seastar video';
%     paths.data = 'G:\My Drive\Andres sea stars\Kinematics';
%     
% % Line to assign single vids    
% elseif ~isempty(dir(['G:\My Drive\Andres sea stars\Kinematics']))
%     
%     paths.vid = 'F:';
%     paths.data = 'G:\My Drive\Andres sea stars\Kinematics';
%     
% elseif ~isempty(dir(['/Users/andrescarrillo/seastar_kine']))
%     
%     paths.vid = '/Volumes/AChd2TB';
%     paths.data = '/Users/andrescarrillo/seastar_kine/dataPath';
%      
% 
%     
% else
%     error('Do not recognize computer')
% end

