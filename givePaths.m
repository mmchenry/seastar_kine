function paths = givePaths
% Returns a structure of paths for running code


if ~isempty(dir(['/Users/mmchenry/Documents/Matlab code']))
    
    paths.vid = '/Users/mmchenry/Documents/Video/Sea stars';
    %paths.vid = '/Volumes/shared/';
    paths.data = '/Users/mmchenry/Documents/Projects/Andres sea stars/Kinematics';
    
    
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

