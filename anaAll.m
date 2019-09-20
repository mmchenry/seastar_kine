function anaAll
% Analysis that takes all sequences into consideration



%% Manage paths (need to modify for new PC)
   
if ~isempty(dir(['/Users/mmchenry/Documents/Matlab code']))

    % Path to root dir of video (CSULB project, external drive)
    vidPath = '/Users/mmchenry/Documents/Video/Chip sea stars/prelim video';
    
    % Path to root of data
    dataPath = '/Users/mmchenry/Documents/Projects/Chip sea stars/prelim data';

elseif isfolder('C:\Users\tpo\Documents\seastar_kine')
    
    % Path to root dir of video (CSULB project, external drive)
    vidPath = 'C:\Users\tpo\Documents\Video\Chip sea stars\prelim video';

    % Path to root of data
    dataPath = 'C:\Users\tpo\Documents\Chip sea star data\prelim data';
    
else
    
    error('Do not recognize computer')
    
end


%% Load data

 % Get list of sequences
cList = dir([dataPath filesep 'S0*']);


% Loop trhu sequences
for i = 1:length(cList)
    
    % Get D structure
    load([dataPath filesep cList(i).name filesep 'summary stats'])
    
    % Store
    A(i).D = D;
    
end




%% Make boxplots

% Get listing of fields
flds = fields(A(1).D);

% Create containers
for i = 1:length(flds)  
    eval([flds{i} '.w = [];'])
    eval([flds{i} '.b = [];'])
end

figure;

% Loop trhu sequences
for i = 1:length(A)
    
    % Loop thru data fields
    for j = 1:length(flds)
        
        % Indicies for w and b fields
        iw = eval(['length([' flds{j} '.w]) + 1']);
        ib = eval(['length([' flds{j} '.b]) + 1']);
        
        % Store, if there is a value for w
        if ~eval(['isempty(A(i).D.' flds{j} '.w)'])
            eval([flds{j} '.w(iw,1) = A(i).D.' flds{j} '.w;'])
        end
        
        % Store, if there is a value for b
        if ~eval(['isempty(A(i).D.' flds{j} '.b)'])
            eval([flds{j} '.b(ib,1) = A(i).D.' flds{j} '.b;'])
        end
    end
end

for i = 1:length(flds)
    
    % Get data
    eval(['data = [' flds{i} '.w; ' flds{i} '.b ];'])
    
    % 
    g = repmat('w',eval(['length(' flds{i} '.w)']),1);
    g = [g; repmat('b',eval(['length(' flds{i} '.b)']),1)];
    
    subplot(length(flds),1,i)
    boxplot(data,g,'labels',{'walk','bounce'},'Colors','k')
    ylabel(flds{i})
    set(gca,'TickDir','out')
    
end

ttt=3;


