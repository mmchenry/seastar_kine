function paths = givePaths
% Returns a structure of paths for running code


if ~isempty(dir(['/Users/mmchenry/Documents/Matlab code']))
    
    paths.vid = '/Users/mmchenry/Documents/Video/Sea stars';
    paths.data = '/Users/mmchenry/Documents/Projects/Andres sea stars/Kinematics';
    
% Line to assign single vids    
elseif ~isempty(dir(['C:\Program Files\MATLAB\R2016a']))
    
    paths.vid = 'C:\Users\andres\Documents\SS Assign';
    paths.data = 'C:\Users\andres\Documents\dataPath';
    
elseif ~isempty(dir(['/Users/andrescarrillo/seastar_kine']))
    
    paths.vid = '/Volumes/AChd2TB';
    paths.data = '/Users/andrescarrillo/seastar_kine/dataPath';
     
else
    error('Do not recognize computer')
end

