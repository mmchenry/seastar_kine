function paths = givePaths
% Returns a structure of paths for running code

% Matt's laptop
if ~isempty(dir(['/Users/mmchenry/Documents/Matlab code']))
    
    % Use flow, if connected
    if ~isempty(dir('/Volumes/shared'))
        paths.vid = '/Volumes/shared';
    else
        paths.vid = '/Users/mmchenry/Documents/Video/Sea stars';
    end
    
    % Use Google Drive, if connected
    if ~isempty(dir('/Volumes/GoogleDrive/My Drive'))
        paths.data = '/Volumes/GoogleDrive/My Drive/Shared files/Projects/Andres sea stars/Kinematics';
    else
        paths.data = '/Users/mmchenry/Documents/Projects/Andres sea stars/Kinematics';
    end
    
elseif ~isempty(dir(['C:\Users\tpo\Documents']))
    
    paths.vid = 'C:\Users\tpo\Documents\Andres seastar video';
    paths.data = 'G:\My Drive\Andres sea stars\Kinematics';
    
% Line to assign single vids    
elseif ~isempty(dir(['G:\My Drive\Andres sea stars\Kinematics']))
    
    paths.vid = 'F:';
    paths.data = 'G:\My Drive\Andres sea stars\Kinematics';
    
elseif ~isempty(dir(['/Users/andrescarrillo/seastar_kine']))
    
    paths.vid = '/Volumes/AChd2TB';
    paths.data = '/Users/andrescarrillo/seastar_kine/dataPath';
     

    
else
    error('Do not recognize computer')
end

