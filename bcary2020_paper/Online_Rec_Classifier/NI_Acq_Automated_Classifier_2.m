function NI_Acq_Automated_Classifier_2

set(0,'DefaultLineLineWidth',1,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',8,...
    'DefaultAxesBox','off',...
    'DefaultAxesFontWeight','Bold');

% edit to match path of commmon dependcies dir
addpath(genpath('G:\bcary2020_paper\common_dependencies'))

% edit to match paths of python on computer and path to 
hs.python_path = 'C:\\Users\\EEG\\Anaconda2\\python';
hs.py_vid_path = 'G:\bcary2020_paper\Online_Rec_Classifier\python_code\AutoVideoAcq.py';

daqreset
disp('Loading NIDAQ Session...')
hs.session = daq.createSession('ni');
a_in1 = addAnalogInputChannel(hs.session,'Dev1', 'ai0', 'Voltage');
a_in2 = addAnalogInputChannel(hs.session,'Dev1', 'ai1', 'Voltage');

a_in1.Range = [-10 10];
a_in2.Range = [-10 10];

hs.samp_rate = 500; % sample rate for EEG and EMG
hs.fps = 2; % frames per sec of video

%%%
% Default thresholds used for automated scoring
hs.nrem_thresh = 100;
hs.rem_thresh = 1;
hs.movt_thresh = 1;
hs.emg_thresh = 1;
%%%

hs.osc_wind_xrange = 0.01;
hs.window_sec = 10; %in sec, length of fft window
hs.gain = 1;
hs.cell_num = 1;
hs.trace_num = 1;
hs.first_run = 1;
hs.trace_name = ['Cell_',num2str(hs.cell_num), '_', num2str(hs.trace_num)];
hs.data_output = struct;
time_data = [];
hs.full_brain_state = [];
hs.wake_history_hr = [];
hs.wake_history_4hr = [];
hs.hist_time_4hr = [];

hs.sample_time = 60; % in sec !Should be divisible by hs.window_sec
hs.vid_buff_sec = 10; % stop video recording for these secs to allow processing
hs.epoch_size = round(hs.sample_time/hs.window_sec);

%%%filter for emg
hs.emg_filt = designfilt('bandpassiir',...
    'FilterOrder',4, ...
    'HalfPowerFrequency1', 5,...
    'HalfPowerFrequency2',150, ...
    'SampleRate',hs.samp_rate);

%%%messaging variables
hs.messaging_onoff = 'off';
hs.recipients = {1112223333,'your email here'};
hs.carrier    = 'sprint';
%%%

hs.save_round = 1;
hs.temp_data_save = [];

screen_size = get(0,'ScreenSize');
hs.s_wid = screen_size(3);
hs.s_ht = screen_size(4);

hs.nidaq_on = 0;

disp('Choose Data Save Folder...')
[hs.dir_path] = uigetdir('Choose Data Save Folder...');

% for debugging
% hs.dir_path = 'C:\Users\EEG\Desktop\Data\temp';
hs.vid_dir_path = [hs.dir_path, filesep, 'video'];

if ~isdir(hs.vid_dir_path)
    disp(['Creating dir: ', hs.vid_dir_path])
    mkdir(hs.dir_path,'video');
end

%%%%%%%%%%%%%%%


%%%%%
hs.bg_color1 = [0.8 0.8 0.8];
hs.bg_color2 = [0.8 0.95 0.8];
hs.bg_color3 = [0.8 0.8 0.95];

function make_fig
    
hs.osc = figure('Visible','on','NumberTitle','off',...
    'Color',hs.bg_color1,'Position',[.48*hs.s_wid .35*hs.s_ht .5*hs.s_wid .5*hs.s_ht],...
    'Visible', 'on');
hs.osc_ax_1 = axes('Parent',hs.osc,...
    'Position', [0.25 0.55 0.7 0.4],...
    'Units', 'Normalized');
hs.osc_ax_2 = axes('Parent',hs.osc,...
    'Position', [0.25 0.08 0.7 0.4],...
    'Units', 'Normalized');
%%%%%
hs.spect = figure('Visible','on','NumberTitle','off',...
    'Color',hs.bg_color2,'Position',[.05*hs.s_wid .2*hs.s_ht .35*hs.s_wid .7*hs.s_ht],...
    'Visible', 'on');
hs.spect_ax = axes('Parent',hs.spect,...
    'Position', [0.05 0.05 0.9 0.9],...
    'Units', 'Normalized');
%%%%%
hs.history = figure('Visible','on','NumberTitle','off',...
    'Color',hs.bg_color3,'Position',[.2*hs.s_wid .1*hs.s_ht .5*hs.s_wid .4*hs.s_ht],...
    'Visible', 'on');
hs.history_ax = axes('Parent',hs.history,...
    'Position', [0.05 0.05 0.9 0.9],...
    'Units', 'Normalized');

%%%%%%%%%%%

uicontrol(hs.osc,'style','text',...
    'string','Sample time (s)',...
    'Units','normalized',...
    'Position', [.025 .92 .15 .03],...
    'FontSize', 10,...
    'Visible', 'on');
    
hs.sample_time_box = uicontrol(hs.osc,'style','edit',...
    'Units','normalized',...
    'string', num2str(hs.sample_time),...
    'FontSize', 14,...
    'Position', [.025 .865 .1 .05],...
    'Visible', 'on');
%%%
uicontrol(hs.osc,'style','text',...
    'string','Ch.1 Gain',...
    'Units','normalized',...
    'Position', [.025 .82 .15 .03],...
    'FontSize', 10,...
    'Visible', 'on');
    
hs.ch1_gain_box = uicontrol(hs.osc,'style','edit',...
    'Units','normalized',...
    'string', '1000',...
    'FontSize', 14,...
    'Position', [.025 .765 .1 .05],...
    'Visible', 'on');
%%%
uicontrol(hs.osc,'style','text',...
    'string','Ch.2 Gain',...
    'Units','normalized',...
    'Position', [.025 .72 .15 .03],...
    'FontSize', 10,...
    'Visible', 'on');
    
hs.ch2_gain_box = uicontrol(hs.osc,'style','edit',...
    'Units','normalized',...
    'string', '5000',... 
    'FontSize', 14,...
    'Position', [.025 .665 .1 .05],...
    'Visible', 'on');
%%%
uicontrol(hs.osc,'style','text',...
    'string','NREM thresh',...
    'Units','normalized',...
    'Position', [.025 .620 .15 .03],...
    'FontSize', 10,...
    'Visible', 'on');
    
hs.dthresh_box = uicontrol(hs.osc,'style','edit',...
    'Units','normalized',...
    'string', num2str(hs.nrem_thresh),...
    'FontSize', 14,...
    'Position', [.025 .565 .1 .05],...
    'Visible', 'on');
%%%
uicontrol(hs.osc,'style','text',...
    'string','REM thresh',...
    'Units','normalized',...
    'Position', [.025 .52 .15 .03],...
    'FontSize', 10,...
    'Visible', 'on');
    
hs.tthresh_box = uicontrol(hs.osc,'style','edit',...
    'Units','normalized',...
    'string', num2str(hs.rem_thresh),...
    'FontSize', 14,...
    'Position', [.025 .465 .1 .05],...
    'Visible', 'on');
%%%
uicontrol(hs.osc,'style','text',...
    'string','EMG thresh',...
    'Units','normalized',...
    'Position', [.025 .42 .15 .03],...
    'FontSize', 10,...
    'Visible', 'on');
    
hs.emgthresh_box = uicontrol(hs.osc,'style','edit',...
    'Units','normalized',...
    'string', num2str(hs.emg_thresh),...
    'FontSize', 14,...
    'Position', [.025 .365 .1 .05],...
    'Visible', 'on');
%%%
uicontrol(hs.osc,'style','text',...
    'string','Movt thresh',...
    'Units','normalized',...
    'Position', [.025 .32 .15 .03],...
    'FontSize', 10,...
    'Visible', 'on');
    
hs.movtthresh_box = uicontrol(hs.osc,'style','edit',...
    'Units','normalized',...
    'string', num2str(hs.movt_thresh),...
    'FontSize', 14,...
    'Position', [.025 .265 .1 .05],...
    'Visible', 'on');
end



make_fig
% %%%%%%%%
%while loop
while ishandle(hs.osc)
    run_sample
end


function run_sample
    disp(clock)
    tic

%     try
    hs.session.Rate = hs.samp_rate;
    osc_wind_size = get(hs.sample_time_box,'String');
    osc_wind_size = str2double(osc_wind_size);

    ch1_gain = get(hs.ch1_gain_box,'String');
    ch1_gain = str2double(ch1_gain);

    ch2_gain = get(hs.ch2_gain_box,'String');
    ch2_gain = str2double(ch2_gain);

    hs.session.DurationInSeconds = osc_wind_size;

    vidsecs = hs.sample_time - hs.vid_buff_sec;
    %%%
    % Start acquisition
    % dos command to start video acquisition
    dos_command = [hs.python_path,' ',hs.py_vid_path,' ',...
        num2str(hs.save_round),' ',hs.dir_path,' ',num2str(vidsecs),' &'];
    hs.dos_status = dos(dos_command);
    
    [data,time] = hs.session.startForeground;
    eeg_data = data(:,2);
    emg_data = data(:,1);
    emg_data = filtfilt(hs.emg_filt, emg_data);

    emg_mult = 1e3; %in mV
    emg_toplot = (emg_data./ch1_gain).*emg_mult;
    plot(hs.osc_ax_1, time, emg_toplot,'k');
    xlabel(hs.osc_ax_1,'Time (secs)');
    ylabel(hs.osc_ax_1,'mV')
    ylim(hs.osc_ax_1, [prctile(emg_toplot,0.0005)-1, prctile(emg_toplot,99.9995)+1])
    box(hs.osc_ax_1,'off')
    title(hs.osc_ax_1,'Ch. 1','Units','Normalized',...
        'Position',[0.0500 1.0096 0.5000])

    eeg_mult = 1e6; %If you want in microVolts
    eeg_data_to_plot = (eeg_data./ch2_gain).*eeg_mult;
    plot(hs.osc_ax_2, time, eeg_data_to_plot,'k');
    xlabel(hs.osc_ax_2,'Time (secs)');
    ylabel(hs.osc_ax_2,'microV')
    ylim(hs.osc_ax_2, [prctile(eeg_data_to_plot,0.0005)-1, prctile(eeg_data_to_plot,99.9995)+1])
    box(hs.osc_ax_2,'off')
    
    %%%spectrogram stuff
    fs = hs.samp_rate;
    title(hs.osc_ax_2,'Ch. 2','Units','Normalized',...
        'Position',[0.0500 1.0096 0.5000])

    
    wind_size = hs.window_sec*fs;
    [s, f, t] = spectrogram(eeg_data, wind_size, 0, [], fs, 'yaxis');

    powers = real(s).^2 + imag(s).^2;
    powers(powers==0) = NaN; %sometimes getting 0 values in col's

    max_freq = 35;
    select_freqs = find(f<=max_freq & f>0.5);
    powers = (powers(select_freqs,:));

    freq_ticks = (f(select_freqs));

    power_to_plot = log10(powers); % take the log10 of powers for better visualization
%     blur_power = imgaussfilt(power_to_plot,0.5);
    hs.spect_ax = subplot(6,1,1:2,'Parent',hs.spect);
    imagesc(t,freq_ticks,power_to_plot,'Parent',hs.spect_ax) % plot the spectrogram results

    % change contrast of color axis
    low_perc = 62;
    high_perc = 98;
    caxis(hs.spect_ax,[prctile(power_to_plot(:),low_perc) prctile(power_to_plot(:),high_perc)])

    set(hs.spect_ax,'YDir','normal')
    spect_ymax = 35;
    ylim(hs.spect_ax,[0 spect_ymax])
    ylabel(hs.spect_ax,'EEG Spect. (Freq. Hz)')

    delta_freqs = find(freq_ticks<=4 & freq_ticks>1);
    theta_freqs = find(freq_ticks<=8 & freq_ticks>5);
    beta_freqs = find(freq_ticks<=34.5 & freq_ticks>20);

    d_power = mean(powers(delta_freqs,:)); 
    t_power = mean(powers(theta_freqs,:));
    b_power = mean(powers(beta_freqs,:));

    d_to_b_ratio = d_power./b_power;
    t_to_d_ratio = t_power./d_power;

    % close dos window
    if hs.dos_status == 0
        dos('taskkill /fi "WindowTitle eq C:\Windows\system32\cmd.exe"');
    else
        pause(3);
        disp('Extra Time for Video...')
        dos('taskkill /fi "WindowTitle eq C:\Windows\system32\cmd.exe"');
    end
    
    %read in video tracking data
    try
        this_round_mat = ['video_',num2str(hs.save_round),'_pymovement','.mat'];
        move_data = load([hs.vid_dir_path,filesep,this_round_mat]);
        raw_move_data = move_data.RAW.raw_movement;

        raw_move_data(1:2) = NaN;
        raw_move_data = [raw_move_data, nan(1,hs.fps*hs.vid_buff_sec)];
%         move_ax = subplot(5,1,5,'Parent',hs.spect);
%         move_time = 1/hs.fps:1/hs.fps:length(raw_move_data)/hs.fps;
%         plot(move_ax, move_time, raw_move_data,'k')
    catch ME
        disp(ME)
        disp('err loading movt data')
        raw_move_data = nan(1,size(powers,2)*hs.window_sec/hs.fps);
%         move_ax = subplot(5,1,5,'Parent',hs.spect);
%         move_time = 1/hs.fps:1/hs.fps:length(raw_move_data)/hs.fps;
%         plot(move_ax, move_time, raw_move_data,'k')
    end
    
    % prepare a movement "average" from weighted combination of avg and median
    move_average = zeros(1,hs.epoch_size);
    wind_avg = hs.fps*hs.window_sec;
    start_ind = 1;
    end_ind = start_ind + wind_avg;
    coef_bias_mean = 0.25;
    movt_mean = nan(1,hs.epoch_size);
    movt_mean(1:floor(length(raw_move_data)/wind_avg)) = slide_window_avg(raw_move_data,wind_avg);

    movt_med = nan(1,hs.epoch_size);
    movt_med(1:floor(length(raw_move_data)/wind_avg)) = ...
        nanmedian(reshape(raw_move_data(1:wind_avg * floor(numel(raw_move_data) / wind_avg)), wind_avg, []), 1);

    move_average = (coef_bias_mean*movt_mean + movt_med)/2;
    
    % Prepare EMG values
    emg_data(isnan(emg_data)) = nanmean(emg_data);
    
    avg_sec = 1;
    avg_int = round(hs.samp_rate*avg_sec);
    std_devs = 4;
    emg_std = nanstd(emg_data); 

    abs_emg = abs(emg_data);
    abs_emg(abs_emg<std_devs*emg_std) = NaN; % remove abs_emg vals below std thresh
    emg_avg = slide_window_avg(abs_emg,avg_int); % apply sliding window (non overlapping)
    emg_avg = emg_avg./emg_std - std_devs; % normalize by std
    emg_avg(isnan(emg_avg)) = 0;

    emg_average = slide_window_avg(emg_avg,hs.window_sec); % emg avg vector for thresh-based algorithm

    %%%%%%%%%%%%%%%
    % Sleep/wake classification

    %brain state:
    %1 = awake
    %2 = rem
    %3 = nrem

    %the standard thresholds for the state continuing after the previous state
    %was the same
    nrem_thresh = get(hs.dthresh_box,'String');
    hs.nrem_thresh = str2double(nrem_thresh);
    rem_thresh = get(hs.tthresh_box,'String');
    hs.rem_thresh = str2double(rem_thresh);
    emg_thresh = get(hs.emgthresh_box,'String');
    hs.emg_thresh = str2double(emg_thresh);
    movt_thresh = get(hs.movtthresh_box,'String');
    hs.movt_thresh = str2double(movt_thresh);
    
    feats = [d_to_b_ratio', t_to_d_ratio', emg_average', move_average'];
    labels = hs.full_brain_state;
    thresholds = [hs.nrem_thresh, hs.rem_thresh, hs.emg_thresh, hs.movt_thresh];
    [this_state_history] = Thresh_State_Algo(t,feats,labels,thresholds);

    plot_state_rects(this_state_history)
    
    %Clean up classification by making a loop to go through states and get rid
    %of REM islands for example

    %%%%%%%%%%%%%%%

%     this_state_history = brain_state(end-length(t)+1:end);
% 
%     nrem_times = t(this_state_history==3)./time_conversion;
%     nrem_ratio = d_to_b_ratio(this_state_history==3);
% 
%     rem_times = t(this_state_history==2)./time_conversion;
%     rem_ratio = t_to_d_ratio(this_state_history==2);

    min_d_toplot = hs.nrem_thresh/10; % arbitrary lower minimum for plotting given use of y log scale
    d_toplot = d_to_b_ratio; 
    d_toplot(d_toplot==0) = min_d_toplot;
    highd_times = t(d_toplot > hs.nrem_thresh); % d ratio pts above thresh
    highd_pts = d_toplot(d_toplot > hs.nrem_thresh);
    hs.delta_ax = subplot(6,1,3,'Parent',hs.spect); % define subplot axis
    cla(hs.delta_ax,'reset')
    hold(hs.delta_ax,'on')
    box(hs.delta_ax,'off')
    plot(hs.delta_ax,t, d_toplot,'Color',[0.2 0.7 1],'LineWidth',1.5)
    plot(hs.delta_ax,highd_times,highd_pts,'.','Color',[0.2 0.7 1],'MarkerSize',22)
    xlim(hs.delta_ax,[0 osc_wind_size])
    rd = refline(hs.delta_ax,0,hs.nrem_thresh);
    rd.LineStyle = '--';
    rd.Color = [0.2 0.7 1];
    try ylim(hs.delta_ax,[min_d_toplot max(d_toplot)]); catch; end
    set(hs.delta_ax,'yscale','log')
    ylabel(hs.delta_ax,'Delta/beta')

    % plot theta ratio same as delta
    min_t_toplot = hs.rem_thresh/10;
    t_toplot = t_to_d_ratio;
    t_toplot(t_toplot==0) = min_t_toplot;
    hight_times = t(t_toplot > hs.rem_thresh);
    hight_pts = t_toplot(t_toplot > hs.rem_thresh);
    hs.theta_ax = subplot(6,1,4,'Parent',hs.spect);
    cla(hs.theta_ax,'reset')
    hold(hs.theta_ax,'on')
    plot(hs.theta_ax,t, t_toplot,'Color',[1 0.25 1],'LineWidth',1.3)
    plot(hs.theta_ax,hight_times,hight_pts,'.','Color',[0.9 0.15 0.9],'MarkerSize',22)
    rt = refline(hs.theta_ax,0,hs.rem_thresh);
    rt.LineStyle = '--';
    rt.Color = [0.9 0.15 0.9];
    hold(hs.theta_ax,'off')
    xlim(hs.theta_ax,[0 osc_wind_size])
    try ylim(hs.theta_ax,[min_t_toplot, max(t_to_d_ratio)]); catch; end
    set(hs.theta_ax,'yscale','log')
    ylabel(hs.theta_ax,'Theta/delta')

    % plot movement/emg values
    max_movt_y = 15; % arbitrary max for visualization only
    highm_times = t(move_average > hs.movt_thresh);
    highm_pts = move_average(move_average > hs.movt_thresh);
    
    hs.move_ax = subplot(6,1,5,'Parent',hs.spect);
    cla(hs.move_ax,'reset')
    hold(hs.move_ax, 'on')
%     yyaxis(hs.move_ax,'left') % prime left y axis
    hs.move_ax.YColor = 'k';
    ylim(hs.move_ax,[0 max_movt_y]);
    ylabel(hs.move_ax,'Movt (pixels)')
    move_time = 1/hs.fps:1/hs.fps:length(raw_move_data)/hs.fps;
    plot(hs.move_ax, move_time, smooth(raw_move_data,4),'k')
    plot(hs.move_ax, t, move_average,'ko','MarkerSize',5)
    plot(hs.move_ax,highm_times,highm_pts,'go','MarkerSize',5)
    rm = refline(hs.move_ax,0,hs.movt_thresh);
    rm.LineStyle = '--';
    rm.Color = 'k';
    
    max_emg_y = 10; % arbitrary max for visualization only
    norm_emg_toplot = emg_average;
    highe_times = t(norm_emg_toplot > hs.movt_thresh);
    highe_pts = norm_emg_toplot(norm_emg_toplot > hs.movt_thresh);
    
    hs.emg_ax = subplot(6,1,6,'Parent',hs.spect);
    cla(hs.emg_ax,'reset')
    hs.emg_color = [0.9 0.07 0.07 0.7];
    ylabel(hs.emg_ax,'EMG (Norm.)')
    ylim(hs.emg_ax,[0 max_emg_y]);
    hold(hs.emg_ax, 'on')
%     yyaxis(hs.move_ax,'right') % prime right y axis
%     hs.emg_ax.YColor = hs.emg_color(1:3);
    plot(hs.emg_ax, t, norm_emg_toplot,'k-')
    plot(hs.emg_ax,highe_times,highe_pts,'o','Color',hs.emg_color,'MarkerSize',5)

    rm = refline(hs.emg_ax,0,hs.emg_thresh);
    rm.LineStyle = '--';
    rm.Color = hs.emg_color;
    hold(hs.emg_ax, 'off')
    box(hs.emg_ax,'off')
    xlim(hs.emg_ax, [0 osc_wind_size]);
    xlabel(hs.emg_ax,'Seconds')
    
    %%%
    % Brain State tracking
    hs.full_brain_state = [hs.full_brain_state this_state_history];
    run_time = length(hs.full_brain_state)*hs.window_sec; %in sec

    try
    history_hours = 1;
    hist_duration_sec = history_hours*3600;
    hist_indx_length = round(history_hours*3600/hs.window_sec);
    % try
    if run_time > hist_duration_sec
        disp('Hour...')
        wake_prop_1 = nanmean(hs.full_brain_state(end-hist_indx_length:end)==1);
        rem_prop_1 = nanmean(hs.full_brain_state(end-hist_indx_length:end)==2);
        nrem_prop_1 = nanmean(hs.full_brain_state(end-hist_indx_length:end)==3);
    else
        wake_prop_1 = nanmean(hs.full_brain_state==1);
        rem_prop_1 = nanmean(hs.full_brain_state==2);
        nrem_prop_1 = nanmean(hs.full_brain_state==3);
    end
    % end

    %%% second alert 
    history_hours = 2;
    hist_duration_sec = history_hours*3600;
    hist_indx_length = (history_hours*3600/hs.window_sec);
    % try
    if run_time > hist_duration_sec
        disp('Two Hour...')
        wake_prop_2 = nanmean(hs.full_brain_state(end-hist_indx_length:end)==1);
        rem_prop_2 = nanmean(hs.full_brain_state(end-hist_indx_length:end)==2);
        nrem_prop_2 = nanmean(hs.full_brain_state(end-hist_indx_length:end)==3);
    else
        wake_prop_2 = NaN;
        rem_prop_2 = NaN;
        nrem_prop_2 = NaN;
    end
    % end

    %%% third alert 
    history_hours = 4;
    hist_duration_sec = history_hours*3600;
    hist_indx_length = (history_hours*3600/hs.window_sec);
    % try
    if run_time > hist_duration_sec
        disp('Four Hour...')
        wake_prop_3 = nanmean(hs.full_brain_state(end-hist_indx_length:end)==1);
        rem_prop_3 = nanmean(hs.full_brain_state(end-hist_indx_length:end)==2);
        nrem_prop_3 = nanmean(hs.full_brain_state(end-hist_indx_length:end)==3);
    else
        wake_prop_3 = NaN;
        rem_prop_3 = NaN;
        nrem_prop_3 = NaN;
    end
    % end

    subject = 'Sleep/Wake History';
    email_message    = {'Now ~ ';...
        ['W:', num2str(nanmean(this_state_history==1),2),...
        ' R:', num2str(nanmean(this_state_history==2),2),...
        ' N:', num2str(nanmean(this_state_history==3),2)];...
        'One Hour (or less)';...
        ['Wake: ', num2str(wake_prop_1)];...
        ['REM: ', num2str(rem_prop_1)];...
        ['NREM: ', num2str(nrem_prop_1)];...
        'Two Hour';...
        ['Wake: ', num2str(wake_prop_2)];...
        ['REM: ', num2str(rem_prop_2)];...
        ['NREM: ', num2str(nrem_prop_2)];...
        'Four Hour';...
        ['Wake: ', num2str(wake_prop_3)];...
        ['REM: ', num2str(rem_prop_3)];...
        ['NREM: ', num2str(nrem_prop_3)]};
    text_message    = {'Now ~ ';...
        ['W:', num2str(nanmean(this_state_history==1),2),...
        ' R:', num2str(nanmean(this_state_history==2),2),...
        ' N:', num2str(nanmean(this_state_history==3),2)];...
        'One Hour (or less) ~ ';...
        ['W:', num2str(wake_prop_1,2),...
        ' R:', num2str(rem_prop_1,2),...
        ' N:', num2str(nrem_prop_1,2)];...
        '2 Hr. ~ ';...
        ['W:', num2str(wake_prop_2,2),...
        ' R:', num2str(rem_prop_2,2),...
        ' N:', num2str(nrem_prop_2,2)];...
        '4 Hr. ~ ';...
        ['W:', num2str(wake_prop_3,2),...
        ' R:', num2str(rem_prop_3,2),...
        ' N:', num2str(nrem_prop_3,2)]};

    if strcmp(hs.messaging_onoff,'on')
        send_msg(hs.recipients(2), subject, email_message, hs.carrier)
        send_msg(hs.recipients(1), {}, text_message, hs.carrier)
    end
    
    disp(email_message)
    
    catch ME
        disp(ME)
        disp('text error')
    end

    %%%%%%
    
    %Track history of wake/sleep ratio
    try
    hs.hist_time = (hs.sample_time/3600):(hs.sample_time/3600):run_time(end)/3600;
    
    hs.wake_history_hr(end+1) = wake_prop_1;
    plot(hs.history_ax, hs.hist_time, hs.wake_history_hr, 'k');
    if exist('wake_prop_3','var')
        hs.wake_history_4hr(end+1) = wake_prop_3;
        hs.hist_time_4hr = [hs.hist_time_4hr hs.hist_time(end)];
        hold(hs.history_ax,'on');
        plot(hs.history_ax, hs.hist_time_4hr, hs.wake_history_4hr, 'r');
        hold(hs.history_ax,'off');        
    end
    ylim(hs.history_ax,[-.1 1.1])
%     
%     plot(hs.history_ax, hs.hist_time, wake_history, 'k');
%     
%     
%     avg_window_mins = 240;
%     smooth_ind = avg_window_mins*60/hs.window_sec;
%     if run_time > hist_duration_sec
%         fr_hr_wake_history = smooth(double(hs.full_brain_state==1),smooth_ind);
%         hold(hs.history_ax,'on');
%         plot(hs.history_ax, hs.hist_time, fr_hr_wake_history, 'r');
%         hold(hs.history_ax,'off');
%     end
    
    catch ME 
        disp(ME)
        disp('history window error')
    end
    


    %%%%%%%%%%%%%%%

    %%% save stuff
    savename = ['data_',num2str(hs.save_round)];
    assignin('base',savename, data)
    save([hs.dir_path filesep savename], 'data')

    try
    time_data{end+1,1} = toc - osc_wind_size;
    time_data{end,2} = clock;

    if hs.save_round ~= 1
        time1 = time_data{end-1,2}(4)*3600 +...
            time_data{end-1,2}(5)*60 +...
            time_data{end-1,2}(6);
        time2 = time_data{end,2}(4)*3600 +...
            time_data{end,2}(5)*60 +...
            time_data{end,2}(6);
        time_data{end,3} = time2 - time1;
    else
        time_data{end,3} = NaN;
    end

    catch ME
        disp(ME)
        disp('Time save error')
    end

    savename = 'time_data';
    assignin('base',savename, time_data)
    save([hs.dir_path filesep savename], 'time_data')
    evalin('base','clear')

    brain_states_save = hs.full_brain_state;
    savename = 'sleep_state_data';
    assignin('base',savename, hs.full_brain_state)
    save([hs.dir_path filesep savename], 'brain_states_save')
    evalin('base','clear')

    hs.save_round = hs.save_round + 1;
    hs.first_run = 0;

    %slow
%     drawnow

    toc
% 
%     catch ME
%         disp(ME)
%         disp('General Error within loop')
%     end

end


%%%
function plot_state_rects(state_data)
    
    y_coord = 31;
    rect_height = 3;
    st_width = hs.window_sec; %in sec
    st_colors = [0 1 0; 0.9 0.15 0.9; 0.2 0.7 1];
    counter = 0;
    for st_num = 1:length(state_data)

        st_width = hs.window_sec;
        if st_num == length(state_data)
            x_coord = (st_num-counter)*st_width - hs.window_sec;
            st_width = st_width*(counter+1);
            rectangle('Position',[x_coord, y_coord,...
                st_width, rect_height],'FaceColor',st_colors(state_data(st_num-1),:),...
                'linestyle','none','Parent',hs.spect_ax);
        elseif (st_num ~= 1) && (state_data(st_num) == state_data(st_num-1))
            counter = counter + 1;
        elseif (st_num == 1) || (state_data(st_num) ~= state_data(st_num-1))
            st_width = hs.window_sec; %in sec
            x_coord = (st_num-counter)*st_width - hs.window_sec;
            st_width = st_width*(counter);
            
            if (st_num ~= 1)
                rectangle('Position',[x_coord, y_coord,...
                    st_width, rect_height],'FaceColor',st_colors(state_data(st_num-1),:),...
                    'linestyle','none','Parent',hs.spect_ax);
            else
                rectangle('Position',[x_coord, y_coord,...
                    st_width, rect_height],'FaceColor',st_colors(state_data(st_num),:),...
                    'linestyle','none','Parent',hs.spect_ax);
            end
            
            counter = 1;
        end
        
    end

end




end






