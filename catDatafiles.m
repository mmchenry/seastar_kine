function cList = catDatafiles(dataPath,camName,datafName)
% Catalogs all video files in directory tree


a1 = dir(dataPath);


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
    
    
    if ~isempty(currOrient)
        
        a2 = dir([dataPath filesep a1(i).name filesep 'SS*']);

        % Loop thru juvenile directories
        for j = 1:length(a2)
            currAge = 'j';
            
            a3 = dir([dataPath filesep a1(i).name filesep a2(j).name filesep ...
                camName filesep 's*']);

            % Individual
            indivNum = str2num(a2(j).name(3:4));
            
            % Catalog, only if data present
            if length(a3)>0

                % Loop trhu individual directories
                for k = 1:length(a3)
                    
                    % Current path
                    currPath = [dataPath filesep a1(i).name filesep a2(j).name filesep ...
                                camName filesep a3(k).name filesep datafName '.mat'];
                    
                    if ~isempty(dir(currPath))
                        
                        % Sequence
                        seqNum = str2num(a3(k).name(end-1:end));
                        
                        % Extract file parts
                        [pathstr,name,ext] = fileparts(currPath);

                        % Store info on inidvidual
                        cList.age(n,1)    = currAge;
                        cList.indiv(n,1)  = indivNum;
                        cList.seq(n,1)    = seqNum;
                        cList.orient(n,1) = currOrient;
                        cList.path{n}     = pathstr;
                        cList.fName{n}    = name;
                        cList.ext{n}      = ext;

                        % Advance index
                        n = n + 1;

                        novid = 0;
                    end
                end
            end
        end
    end
end


if novid==1
    warning('No data files found')
    cList = [];
end
