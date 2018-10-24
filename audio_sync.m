function delay = audio_sync(c_path,s_path)
% Save or load structure for syncing 3 video files from audio channels

%% Code execution

% Plot data to verify correct sync
plot_data = 1;

% Choose a portion of the audio
select_event = 1;


%% Parameters

% Interval for evaluating max intensity
dur_eval = 1;

% Duration of audio for visualizing the correction
%audio_dur = read_dur ;


%% Combine two signals

% Get audioinfo
aInfo_c = audioinfo(c_path);
aInfo_s = audioinfo(s_path);

% Extract audio
[tmp_c,Fs_c] = audioread(aInfo_c.Filename);
[tmp_s,Fs_s] = audioread(aInfo_s.Filename);

% Choose channel with greater signal (canon)
if mean(abs(tmp_c(:,1))) > mean(abs(tmp_c(:,2)))
    y_c = tmp_c(:,1);
else
    y_c = tmp_c(:,2);
end

% Choose channel with greater signal (sony)
aLevel = 0;
for i = 1:size(tmp_s,2)
    
    cVal = mean(abs(tmp_s(:,i)));
    
    if cVal > aLevel
        aLevel = cVal;
        iChannel = i;
    end
end

% Store best channel
y_s = tmp_s(:,iChannel);
 
% Normalize channels by maximum
y_s = y_s./max(abs(y_s));
y_c = y_c./max(abs(y_c));

% Time values
t_c = [0:(length(y_c)-1)]'./Fs_c;
t_s = [0:(length(y_s)-1)]'./Fs_s;

% Index of values
idx = 1:min([length(t_c) length(t_s)]);

% Audio matrix
y = [y_c(idx) y_s(idx)];

% Time matrix
if (Fs_c==Fs_s)
    t = t_c(idx);
    Fs = Fs_c;
else
    error('Unequal sample rates');
end
    
clear Fs_s Fs_c t_s t_c y_s y_c tmp_c tmp_s aInfo_c aInfo_s idx iChannel cVal
    

%% Find interval with max signal

% % Loop thru intervals to find mean values
% for j = 1:floor(max(t(:))/dur_eval)
%     
%     % Index of interval values
%     idx = (t>=(dur_eval*(j-1))) & (t<(dur_eval*j));
%     
%     % Mean product of  cameras over interval
%     yMean(j,1) = mean(abs(y(idx,1)).*abs(y(idx,2)));
%     
%     % Mean time value
%     tMean(j,1) = mean(t(idx));
% end
%     
% % Time index
% iTime = find(yMean==max(yMean),1,'first');
% 
% % Start time
% t_start = tMean(iTime) - 1.5*dur_eval;
% 
% % End time
% t_end = tMean(iTime) + 1.5*dur_eval;
% 
% % Index of values to interrogate
% idx = (t>=t_start) & (t<t_end);
% 
% % Trim to duration to be considered
% y = y(idx,:);
% t = t(idx,:);

% Find delays wrt first camera
delay = finddelay(y(:,1),y(:,2))./Fs;

% Store data
% aud.date_dir    = audio{i}.date_dir;
% aud.seq_dir     = audio{i}.seq_dir;
% aud.vid_file    = audio{i}.vid_file;
% aud.delay       = delay;

% Write data
%save([cDir filesep 'audio_data'],'aud');

if 0
    % Play audio (for troubleshooting)
    aud_player = audioplayer(y(:,1),Fs);
    play(aud_player)
end

% VISUALIZE RESULTS ---------------------------------------------------

figure;
%     subplot(3,1,1)
%     plot(tMean,yMean);
%     xlabel('t (s)')
%     ylabel('Mean audio (V)');

subplot(2,1,1)
plot(t,y);
xlabel('t (s)')
ylabel('Audio intensity');
title('Raw audio')
legend('Canon','Sony')


subplot(2,1,2)
plot(t,y(:,1),'-',...
    t-delay,y(:,2),'-');
xlabel('t (s)')
ylabel('Audio intensity');
title('Corrected for delay')

if delay>0
    disp(['Sony follows Canon by ' num2str(delay) ' s'])
else
    disp(['Canon follows Sony by ' num2str(-delay) ' s'])
end

pause(0.1)

% Capture graphs
%I = getframe(gcf);
%close

% Write frame
%imwrite(I.cdata,[cDir filesep 'Audio sync.jpeg'],'JPEG');
%close

% Update status
%disp(['Completed ' num2str(i) ' of ' num2str(length(audio))])







function prompt_play(aud,ch_num)

buttonName = questdlg('Play weak channel?',['Channel ' ch_num],...
                      'Yes','No','Cancel','Yes');

switch buttonName
    
    case 'Cancel'
        return
    case 'No'
        % Do nothing
    case 'Yes'
        disp(['Playing audio for the following:'])
        disp(aud.path)
        play(aud.player)
        
end



    