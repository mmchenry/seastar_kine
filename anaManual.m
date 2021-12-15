function anaManual
% Analyzes manually collected data





%% Read data from catalog spreadsheet

if do.reImportData || do.anaAudioSync
    
    % Get list of all side csv files
%     dlc_path = [paths.data filesep 'data' filesep 'side_view' filesep '*.csv'];
    dlc_path = [paths.data filesep 'rawCSV' filesep '*.csv'];
    dlc_files = dir(dlc_path);
    
    % Check for dlc files
    if isempty(dlc_files)
        error(['No DLC csv files found in ' dlc_path]);
    end
    
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
        %     if T.use_video(i)==1 && T.ana_video(i)==1 && T.complete_video(i)==0
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
            seq(j).calConst       = T.cal_side(i);
            seq(j).fps            = T.frame_rate_side(i);
            seq(j).SW_percent     = T.percent_sw(i);
            
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
end
