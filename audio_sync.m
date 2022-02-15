function [delay,info] = audio_sync(side_path,bot_path)
% Save or load structure for syncing 3 video files from audio channels

%% Code execution

% Plot data to verify correct sync
plot_data = 0;

% Audio channel(s) to monitor for syncing
aChan = [3 4];


%% Parameters

% Interval for evaluating max intensity
dur_eval = 1;

% Duration of audio for visualizing the correction
%audio_dur = read_dur ;


%% Combine two signals

% Get audioinfo
aInfo_side = audioinfo(side_path);
aInfo_bot = audioinfo(bot_path);

% Extract audio
[tmp_side,Fs_side] = audioread(aInfo_side.Filename);
[tmp_bot,Fs_bot]   = audioread(aInfo_bot.Filename);

% Choose channel with greater signal (canon)
% if mean(abs(tmp_side(:,1))) > mean(abs(tmp_side(:,2)))
%     y_side = tmp_side(:,1);
% else
%     y_side = tmp_side(:,2);
% end
y_side = mean(tmp_side(:,aChan),2);


% % Choose channel with greater signal (sony)
% aLevel = 0;
% for i = 1:size(tmp_bot,2)
%     
%     cVal = mean(abs(tmp_bot(:,i)));
%     
%     if cVal > aLevel
%         aLevel = cVal;
%         iChannel = i;
%     end
% end

% Store best channel
% y_bot = tmp_bot(:,iChannel);

y_bot = mean(tmp_bot(:,aChan),2);
 
% Normalize channels by maximum
y_bot = y_bot./max(abs(y_bot));
y_side = y_side./max(abs(y_side));

% Time values
t_side  = [0:(length(y_side)-1)]'./Fs_side;
t_bot   = [0:(length(y_bot)-1)]'./Fs_bot;

% Index of values
idx = 1:min([length(t_side) length(t_bot)]);

% Audio matrix
y = [y_side(idx) y_bot(idx)];

% Time matrix
if (Fs_side==Fs_bot)
    t = t_side(idx);
    Fs = Fs_side;
else
    error('Unequal sample rates');
end
    
clear Fs_bot Fs_side t_bot t_side y_bot y_side tmp_side tmp_bot aInfo_side aInfo_bot idx iChannel cVal
    

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
% t_bottart = tMean(iTime) - 1.5*dur_eval;
% 
% % End time
% t_end = tMean(iTime) + 1.5*dur_eval;
% 
% % Index of values to interrogate
% idx = (t>=t_bottart) & (t<t_end);
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

if delay>0
    info = ['Bottom follows Side by ' num2str(delay) ' s'];
else
    info = ['Side follows Bottom by ' num2str(delay) ' s'];
end

% VISUALIZE RESULTS ---------------------------------------------------
if plot_data
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
    plot(t+delay,y(:,1),'-',...
        t,y(:,2),'-');
    xlabel('t (s)')
    ylabel('Audio intensity');
    title('Corrected for delay')
    pause(0.1)
end

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



    