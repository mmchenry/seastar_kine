function batchMaster
% Handles the acquisition of data for all sequences


% Extension for viode file names
extVid = 'MOV';


%% Manage paths 
   
paths = givePaths;


%% Read data from catalog spreadsheet

% Throw error if spreadsheet not in place
if exist([paths.data filesep 'Weights_experiments.xlsx'],'file')==0
    error(['Missing spreadsheet: ' paths.data filesep ...
           'Weights_experiments.xlsx'])
end

% Read table
T = readtable([paths.data filesep 'Weights_experiments.xlsx']);

j = 1;seq = [];

% Step thru rows
for i = 1:length(T.date)
    
    % If using video at all and currently analyzing
    if T.use_video(i)==1 && T.ana_video(i)==1
       
        seq(j).dateNum        = datenum(T.date(i));
        seq(j).ext            = extVid;
        seq(j).fName_side     = [T.side_filename{i}];
        seq(j).fName_bot      = [T.bottom_filename{i}];
        seq(j).fName_calSide  = [T.side_cal_filename{i}];
        seq(j).fName_calBot   = [T.bot_cal_filename{i}];
        seq(j).expType        = T.exp_type{i};
        seq(j).indiv          = T.indiv_num(i);
        seq(j).addMass        = T.added_mass(i);
        seq(j).floatNum       = T.float_num(i);
        seq(j).bodyMass       = T.body_mass(i);
        
        % Check experiment type
        if strcmp(seq(j).expType,'c')
            seq(j).dirName = 'control';
        elseif strcmp(seq(j).expType,'w')
            seq(j).dirName = 'weights';
        elseif strcmp(seq(j).expType,'f')
            seq(j).dirName = 'floats';
        else
            error(['Do not recognize experiment type: ' seq(i).expType])
        end
        
        j = j + 1;
    end 
end

% Check for seq
if isempty(seq)
    error('No sequence data found from spreadsheet');
end

clear j T


%% Survey video and data files for match to sequences requested

% Check for presence of video directories
if ~isfolder([paths.vid filesep 'calibration'])
    error(['Missing folder: ' paths.vid filesep 'calibration'])
    
elseif ~isfolder([paths.vid filesep seq(i).dirName])
    error(['Missing folder: ' paths.vid filesep seq(i).dirName])
end

k = 1;

% Loop thru sequences to find matches among files
for i = 1:length(seq)
    
    % Video files for side view
    aSide = dir([paths.vid filesep seq(i).dirName filesep 'side' ...
                 filesep 'STUDIO*']);
    
    % Video files for bottom view
    aBot = dir([paths.vid filesep seq(i).dirName filesep 'bottom' ...
                filesep 'STUDIO*']);
    
     % Video files for bottom view
    aCalBot = dir([paths.vid filesep 'calibration' filesep 'bottom' ...
                   filesep 'STUDIO*']);
    
     % Video files for side view
    aCalSide = dir([paths.vid filesep 'calibration' filesep 'side' ...
                    filesep 'STUDIO*']);
    
    % Indicies 
    iMatch_bot      = nan;
    iMatch_side     = nan;
    iMatch_calSide  = nan;
    iMatch_calBot   = nan;
    
    % Loop thru bottom video files for match
    for j = 1:length(aBot)    
        if strcmp(aBot(j).name,[seq(i).fName_bot '.' seq(i).ext])
            iMatch_bot = j;
            break
        end
    end
    
    % Loop thru side video files for match
    for j = 1:length(aSide)    
        if strcmp(aSide(j).name,[seq(i).fName_side '.' seq(i).ext])
            iMatch_side = j;
            break
        end
    end
    
    % Loop thru bottom video files for match
    for j = 1:length(aCalBot)    
        if strcmp(aCalBot(j).name,[seq(i).fName_calBot '.' seq(i).ext])
            iMatch_calBot = j;
            break
        end
    end
    
    % Loop thru side video files for match
    for j = 1:length(aCalSide)    
        if strcmp(aCalSide(j).name,[seq(i).fName_calSide '.' seq(i).ext])
            iMatch_calSide = j;
            break
        end
    end
    
    % If there is a match of video to spreadsheet item
    if ~isnan(iMatch_bot)    &&  ~isnan(iMatch_side)    && ...
       ~isnan(iMatch_calBot) &&  ~isnan(iMatch_calSide)     
        
        % Store info
        vid_bot(k)      = aBot(iMatch_bot);
        vid_side(k)     = aSide(iMatch_side);
        vid_calBot(k)   = aCalBot(iMatch_calBot);
        vid_calSide(k)  = aCalSide(iMatch_calSide);
        
        % Advance index
        k = k + 1;
    end
    
    % If no match to bottom
    if isnan(iMatch_bot)
        warning(['No bottom video file matching ' seq(i).fName_bot ...
            ' from spreadsheet'])
    end
       
    % If no match to side
    if isnan(iMatch_side)
        warning(['No side video file matching ' seq(i).fName_side ...
            ' from spreadsheet'])
    end
    
    % If no match to bottom
    if isnan(iMatch_calBot)
        warning(['No bottom cal video file matching ' seq(i).fName_calBot ...
            ' from spreadsheet'])
    end
       
    % If no match to side
    if isnan(iMatch_calSide)
        warning(['No side cal video file matching ' seq(i).fName_calSide ...
            ' from spreadsheet'])
    end
    
    clear aSide aBot iMatch_bot iMatch_side j iMatch_calBot iMatch_calSide
end

% Check that there are matches
if ~exist('vid_side','var')
    error('No matches between spreadsheet and side videos . . .')
    
elseif ~exist('vid_bot','var')
    error('No matches between spreadsheet and bottom videos . . .')

elseif ~exist('vid_calBot','var')
    error('No matches between spreadsheet and bottom calibration videos . . .')
    
elseif ~exist('vid_calSide','var')
    error('No matches between spreadsheet and side calibration videos . . .')
end

clear i k


%% Make data directories, if necessary

% Make calibration folder, if not present
if ~isfolder([paths.data filesep 'calibration'])
    mkdir([paths.data filesep 'calibration'])
end

% Loop thur sequences
for i = 1:length(seq)
    
    % Make dirName, if not present
    if ~isfolder([paths.data filesep seq(i).dirName])
        mkdir([paths.data filesep seq(i).dirName])
        mkdir([paths.data filesep seq(i).dirName filesep 'side'])
        mkdir([paths.data filesep seq(i).dirName filesep 'bottom'])
    end
        
    % Current side directory
    currSide = [paths.data filesep seq(i).dirName filesep 'side' ...
                filesep seq(i).fName_side];
    
    % Make, if not there
    if ~isfolder(currSide)
        mkdir(currSide)
    end
    
    % Current side directory
    currBot = [paths.data filesep seq(i).dirName filesep 'bottom' ...
               filesep seq(i).fName_bot];
    
    % Make, if not there
    if ~isfolder(currBot)
        mkdir(currBot)
    end
    
    clear currSide currBot    
end


%% Make movies for DeepLabCut

% Loop thru sequences
for i = 1:length(seq)
    
    % Current directories
    dataPath = [seq(i).dirName filesep 'bottom' filesep seq(i).fName_bot];
    vidPath  = [seq(i).dirName filesep 'bottom' filesep seq(i).fName_bot  ...
                '.' seq(i).ext];
    
    % Run analysis
    movieMaster(dataPath,vidPath,'deep movie')
    
end







