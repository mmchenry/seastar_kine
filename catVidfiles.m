function cList = catVidfiles(vidPath,camName)
% Catalogs all video files in directory tree


a1 = dir(vidPath);


n = 1;
novid = 1;

% Loop thru oriention directories
for i = 1:length(a1)
    
    if a1(i).isdir && strcmp(a1(i).name,'Horizontal')
        
        currOrient = 'h';
        
    elseif a1(i).isdir && strcmp(a1(i).name,'Vertical')
        
        currOrient = 'v';
        
    elseif a1(i).isdir && ...
            (strcmp(a1(i).name,'Upside-down') || strcmp(a1(i).name,'UpsideDown'))
        
        currOrient = 'u';
        
    else
        currOrient = [];
    end
    
   % pause(0.1)
    
    if ~isempty(currOrient)
        
        a2 = dir([vidPath filesep a1(i).name filesep 'SS*']);

        % Loop thru juvenile directories
        for j = 1:length(a2)
            currAge = 'j';
            
            a3 = dir([vidPath filesep a1(i).name filesep a2(j).name filesep ...
                camName filesep '*.MOV']);
         
            indivNum = str2num(a2(j).name(3:4));
            
            cal.Type = [];
            
            m = 1;
            
            % Catalog, only if movies present
            if length(a3)>0
                
                movPath = [];
                movType = [];
                
                % Loop trhu individual directories
                for k = 1:length(a3)
                    
                    %                 if a3(k).isdir && length(a3(k).name)>2 &&...
                    %                         strcmp(a3(k).name(1:2),'SS')
                    %
                    %
                    %                     % Directory contents for individual
                    %                     a4 = dir([vidPath filesep a1(j).name filesep ...
                    %                         a2(j).name filesep a3(k).name]);
                    %
                    %                     m = 1;cal.Type=[];
                    
                    % Loop thru sequences for the indiviual
                    
                    % Extract file parts
                    %[pathstr,name,ext] = fileparts(movPath{f});
                    
                    % If video is a calibration
                    if length(a3(k).name) > 7 && ...
                            (strcmp(a3(k).name,'scale.MOV') || ...
                            strcmp(a3(k).name,'Scale.MOV') || ...
                            strcmp(a3(k).name,'Calibration.MOV') ||...
                            strcmp(a3(k).name,'calibration.MOV'))
                        
                        if strcmp(a3(k).name(end-2:end),'MOV')
                            cal.Type = 'mov';
                            
                        else
                            error('no match for calibration type');
                        end
                        
                        % Store calibration path
                        cal.Path = [a1(i).name filesep ...
                            a2(j).name filesep a3(k).name];
                        
                        % If video is an MOV file . . .
                    elseif strcmp(a3(k).name(1:2),'s0') && strcmp(a3(k).name(end-2:end),'MOV')
                        
                        movPath{m} = [a1(i).name filesep ...
                            a2(j).name filesep camName filesep a3(k).name ];
                        movType{m} = 'mov';
                        
                        m = m + 1;
                        %pause(0.1)
                    end
                end

                if ~isempty(movPath)
                    % Loop trhu sequences
                    for f = 1:length(movType)
                        
                        % Extract file parts
                        [pathstr,name,ext] = fileparts(movPath{f});
                        
                        % Store info on inidvidual
                        cList.age(n,1)    = currAge;
                        cList.indiv(n,1)  = indivNum;
                        cList.orient(n,1) = currOrient;
                        
                        cList.vidType{n} = movType{f};
                        cList.path{n}    = pathstr;
                        cList.fName{n}   = name;
                        cList.ext{n}     = ext;
                        
                        novid = 0;
                        
                        if isempty(cal.Type)
                            cList.calPath{n} = [];
                            %disp(' ')
                            warning(['No calibration for: ' movPath{f}]);
                        else
                            cList.calPath{n} = cal.Path;
                        end
                        
                        % Advance sequence index
                        n = n + 1;
                    end
                end
                
                
            else
                warning(['No videos in: ' vidPath filesep a1(i).name ...
                         filesep a2(j).name filesep camName ])
            end
        end
    end
end


if novid==1
    warning('No video files found')
    cList = [];
end
