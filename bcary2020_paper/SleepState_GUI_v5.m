function Sleep_state_coder_v4_SLEEPRIG

screen_size = get(0,'ScreenSize');
s_wid = screen_size(3);
s_ht = screen_size(4);

set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',12,...
    'DefaultAxesBox','off',...
    'DefaultAxesFontWeight','Bold');

% path for folder containing raw data
start_path = 'C:\Users\SleepRig-2\Desktop\Data\';

% folder containing python functions used for random forest (RF) algorithm
pywrap_dir = 'C:\Users\SleepRig-2\Documents\MATLAB\SleepState_Code\python_callers';
addpath(genpath(pywrap_dir));

% folder containing python program
python_path = 'C:\\Users\\SleepRig-2\\Anaconda2\\python';
addpath(genpath(python_path));

% path used for temporary python RF files
hs.ml_model_path='C:\Users\SleepRig-2\Documents\MATLAB\tempPyVidcodeFiles\';

%%%%%%%%
%Params%
%%%%%%%%

% global variables used for analysis. These need to be set correctly
% according to the format of raw input files
hs.sample_rate  = 500;
hs.viewer_mins = 16; %length in min of GUI viewing window
hs.filesave_mins = 4; % length in min of eeg/emg files
hs.num_files_to_load = round(hs.viewer_mins/hs.filesave_mins);
hs.trace_length = hs.filesave_mins*hs.sample_rate*60;
hs.data_end = hs.viewer_mins*60*hs.sample_rate;

% video fps
hs.fps = 2;
movt_data_length_min = 4; % length in min of video files
hs.movt_data_length = hs.fps*movt_data_length_min*60;
hs.num_vidfiles_to_load = round(hs.viewer_mins/movt_data_length_min);
hs.length_movt_toload = hs.viewer_mins*60*hs.fps;
hs.video_beg_secs = [];

hs.sleep_state_data = [];
hs.window_sec = 10; % sec for each epoch used for classification
hs.epoch_size = hs.viewer_mins*(60/hs.window_sec);
hs.coords_chosen = 0;
hs.ind1 = 0;
hs.ind2 = 0;
hs.time_conv = 60; % 3600 for hours, 60 for mins
hs.have_analyzed = [];
hs.num_inds_changd = 0;
hs.running_ep_num = 0;

% plot emg or gamma power in GUI --> 1=EMG, 0=gamma
hs.EMG_toplot = 1;
hs.emg_color = [0.9 0.07 0.07 0.7];

hs.first_run = 1;

% filter for emg
hs.emg_filt = designfilt('bandpassiir',...
    'FilterOrder',4, ...
    'HalfPowerFrequency1', 20,...
    'HalfPowerFrequency2',200, ...
    'SampleRate',hs.sample_rate);

% variables concerning RF machine learning (ML). 
hs.num_prev_states = 3;
hs.num_2nd_for_states = 2;
hs.ML_onoff = 0;
hs.auto_onoff = 0;

hs.mlmode = 'RF';
hs.nTrees = 200; % number of random forest trees used in the ML model
hs.min_state_size = 10; % minimum seconds of state epochs to for final output

%ML Vars
hs.ML_SSD = [];
hs.ML_SSD_prev1 = [];
hs.ML_SSD_prev2 = [];
hs.ML_delta = [];
hs.ML_theta = [];
hs.g_to_d_ratio = [];
hs.ML_movt = [];
hs.ML_EMG = [];

hs.full_confidence = [];

hs.state_sec = 1;
hs.upsamp_by = hs.window_sec/hs.state_sec;

%%%
% Can set starting hour. Default should be 0.
hs.hour_start = 0;
hs.epoch_num_start = floor(hs.hour_start*(60/hs.viewer_mins));

hs.darkhour = 19; % 0-24 hr at which the lights turn off in the video

disp(['Starting close to hour: ',num2str(hs.hour_start)])

%%%
% Default thresholds used for automated scoring
hs.nrem_thresh = 100;
hs.rem_thresh = 1;
hs.movt_thresh = 1;
hs.emg_thresh = 4;
hs.volt_thresh = 2000; %microV
%%%

% function that creates the GUI. Creates button variables with callbacks
function create_window
hs.fig = figure('Visible','on','NumberTitle','off',...
    'Color',[.8 .8 .95],'Position',[.15*s_wid .25*s_ht .65*s_wid .6*s_ht],...
    'KeyPressFcn', @key_catcher);

% switch ML/thresh algo
hs.switch_mode_but = uicontrol(hs.fig,'Style','togglebutton','String','ML Off',...
        'Units','normalized', 'Callback', @switch_mode_but,...
        'Position',[.92 0.71 0.065 0.1]);
if hs.ML_onoff
    set(hs.switch_mode_but,'String','ML On')
else
    set(hs.switch_mode_but,'String','ML Off')
end

% train ML
hs.trainML_but = uicontrol(hs.fig,'Style','pushbutton','String','Train ML',...
        'Units','normalized', 'Callback', @trainML_but,...
        'Position',[.92 0.59 0.065 0.1]);
% next
hs.next_but = uicontrol(hs.fig,'Style','pushbutton','String','Next',...
        'Units','normalized', 'Callback', @next_but,...
        'Position',[.92 0.47 0.065 0.1]);
% back
hs.back_but = uicontrol(hs.fig,'Style','pushbutton','String','Back',...
        'Units','normalized', 'Callback', @back_but,...
        'Position',[.92 0.35 0.065 0.1]);

hs.reload_but = uicontrol(hs.fig,'Style','pushbutton','String',{'Reload'},...
        'Units','normalized', 'Callback', @reload_next_but,...
        'Position',[.92 0.05 0.065 0.1]);

hs.auto_but = uicontrol(hs.fig,'Style','togglebutton','String',{'Auto Off'},...
        'Units','normalized', 'Callback', @auto_next_but,...
        'Position',[.92 0.17 0.065 0.075]);    
if hs.auto_onoff
    set(hs.auto_but,'String','Auto On')
else
    set(hs.auto_but,'String','Auto Off')
end    

% text instructions
hs.text_instructions = uicontrol(hs.fig,'Style','text','String','Wake (w/p), NREM (e/[), REM (r/]), change (space/o)',...
        'Units','normalized',...
        'Position',[0.004 0.81 0.085 0.1]);

%%%
uicontrol(hs.fig,'style','text',...
    'string','Delta Thresh',...
    'Units','normalized',...
    'Position', [0.004 0.56 0.05 0.061],...
    'FontSize', 10,...
    'Visible', 'on');
    
hs.dthresh_box = uicontrol(hs.fig,'style','edit',...
    'Units','normalized',...
    'string', num2str(hs.nrem_thresh),...
    'FontSize', 12,...
    'Position', [.058 .56 .04 .061],...
    'Visible', 'on');
%%%
uicontrol(hs.fig,'style','text',...
    'string','Theta Thresh',...
    'Units','normalized',...
    'Position', [0.004 0.49 0.05 0.061],...
    'FontSize', 10,...
    'Visible', 'on');
    
hs.tthresh_box = uicontrol(hs.fig,'style','edit',...
    'Units','normalized',...
    'string', num2str(hs.rem_thresh),...
    'FontSize', 12,...
    'Position', [.058 .49 .04 .061],...
    'Visible', 'on');
%%%
uicontrol(hs.fig,'style','text',...
    'string','Movt. Thresh',...
    'Units','normalized',...
    'Position', [0.004 0.42 0.05 0.061],...
    'FontSize', 10,...
    'Visible', 'on');
    
hs.movtthresh_box = uicontrol(hs.fig,'style','edit',...
    'Units','normalized',...
    'string', num2str(hs.movt_thresh),...
    'FontSize', 12,...
    'Position', [.058 .42 .04 .061],...
    'Visible', 'on');
%%%
uicontrol(hs.fig,'style','text',...
    'string','EMG Thresh',...
    'Units','normalized',...
    'Position', [0.004 0.35 0.05 0.061],...
    'FontSize', 10,...
    'Visible', 'on');
    
hs.emgthresh_box = uicontrol(hs.fig,'style','edit',...
    'Units','normalized',...
    'string', num2str(hs.emg_thresh),...
    'FontSize', 12,...
    'Position', [.058 .35 .04 .061],...
    'Visible', 'on');
%%%
uicontrol(hs.fig,'style','text',...
    'string','microV Thresh',...
    'Units','normalized',...
    'Position', [0.004 0.28 0.05 0.061],...
    'FontSize', 10,...
    'Visible', 'on');
    
hs.voltthresh_box = uicontrol(hs.fig,'style','edit',...
    'Units','normalized',...
    'string', num2str(hs.volt_thresh),...
    'FontSize', 12,...
    'Position', [.058 .28 .04 .061],...
    'Visible', 'on');
end

%%%
% create the first GUI window
create_window
drawnow
%%%

% load a finder window to pick a dir containing the raw EEG/EMG files with
% video files in a subdir
[hs.dir_path] = uigetdir('Find Raw EEG Files...');
hs.directory = dir(hs.dir_path);

[~, dir_name] = fileparts(hs.dir_path);
hs.anim = dir_name; % use dir name as the save name for the experiment
hs.anim_ssd = [hs.anim,'_ssd'];
disp(['Exp. save name is: ',hs.anim])

hs.mat_files = dir([hs.dir_path, filesep,'data_*']); %collect mat files in dir
mat_datenums = [];
for ii = 1:length(hs.mat_files)
    mat_datenums(ii) = round(hs.mat_files(ii).datenum,10);
end
mat_datenums = sort(mat_datenums);
exp_length = (mat_datenums(end)-mat_datenums(1))*24*3600 + hs.filesave_mins*60;

hs.num_data_files = length(mat_datenums);
hs.have_analyzed = zeros(1,round(hs.num_data_files/hs.num_files_to_load + 1));

% loop to load all of the recording's data at once and then shifts the
% files accordings to datenum's of the mat files in the dir
hs.full_eeg_data = nan(1,round(exp_length*hs.sample_rate) + hs.trace_length);
hs.full_emg_data = nan(1,round(exp_length*hs.sample_rate) + hs.trace_length);
textprogressbar(char('Loading Matlab Files: '));
for mat_num = 1:length(hs.mat_files)

    try
    textprogressbar(mat_num*100/(length(hs.mat_files)))
    if mat_num == 1
        mat_end_time = hs.filesave_mins*60;
        first_abs_end_time = mat_datenums(mat_num)*24*3600;
    else
        mat_end_time = mat_datenums(mat_num)*(24*3600)...
            - first_abs_end_time + hs.filesave_mins*60;
    end
     
    this_eeg_name = ['data_',num2str(mat_num),'.mat'];
    this_data = load([hs.dir_path,filesep,this_eeg_name]);
    this_data = this_data.data;
    
    start_ind = round((mat_end_time*hs.sample_rate-length(this_data)+1));
    end_ind = round(mat_end_time*hs.sample_rate);
    shift_ind = start_ind:end_ind;
    
    hs.full_eeg_data(shift_ind) = this_data(:,2);
    hs.full_emg_data(shift_ind) = this_data(:,1);
    catch
        disp('loading err')
        textprogressbar('done');
    end
end
textprogressbar('done');

% convert to mV from V
hs.full_eeg_data = hs.full_eeg_data*1000;
hs.full_emg_data = hs.full_emg_data*1000;

% find video dir which should be labeled simply 'video'
hs.vid_dir_path = [hs.dir_path, filesep, 'video'];
vid_directory = dir(hs.vid_dir_path);

hs.csv_files = dir([hs.vid_dir_path, filesep,'*.csv']);
hs.avi_datenums = [];
hs.full_vid_fr_times = [];
hs.avi_files = dir([hs.vid_dir_path, filesep,'*.avi']);
hs.avi_datenums = [];
hs.avi_filenames = {};
no_vid_secs = 10; % seconds of video loss each data file
fr_time_length = hs.movt_data_length-no_vid_secs*hs.fps;

% Load in the video file frame times from csv
textprogressbar(char('Loading Movement Times: '));
try
for ii = 1:length(hs.csv_files)
    textprogressbar(ii*100/(length(hs.csv_files)))
    csv_path = [hs.csv_files(ii).folder filesep hs.csv_files(ii).name];
    fr_times = csvread(csv_path);
    
    if length(fr_times) > fr_time_length
        fr_times = linspace(fr_times(1),fr_times(end),fr_time_length)';
    elseif length(fr_times) < fr_time_length
        fr_times(end:fr_time_length) = NaN;
    end
    hs.full_vid_fr_times = [hs.full_vid_fr_times fr_times'];
    
    hs.avi_datenums(ii) = round(hs.avi_files(ii).datenum,10);
    hs.avi_filenames{ii} = hs.avi_files(ii).name;
end
catch
    textprogressbar('crashed');
end
textprogressbar('done');

hs.full_vid_fr_times = sort(hs.full_vid_fr_times);
hs.full_vid_fr_times = hs.full_vid_fr_times - hs.full_vid_fr_times(1);
[hs.avi_datenums sort_ind] = sort(hs.avi_datenums);
hs.avi_filenames = hs.avi_filenames(sort_ind);

hs.full_movt_data = [];

% load in video file movement mat files from python motion detection
try
textprogressbar(char('Loading Movement Files: '));
for avi_num = 1:length(hs.avi_datenums)
        
    textprogressbar(avi_num*100/(length(hs.avi_datenums)))
    try
        if avi_num == 1
            avi_end_time = movt_data_length_min*60;
            first_abs_end_time = hs.avi_datenums(avi_num)*(24*3600);
        else
            avi_end_time = hs.avi_datenums(avi_num)*(24*3600)...
                - first_abs_end_time + movt_data_length_min*60;
        end

        this_movt_name = ['video_',num2str(avi_num),'_pymovement','.mat'];
        this_move_data = load([hs.vid_dir_path,filesep,this_movt_name]);
        raw_move_data = this_move_data.RAW.raw_movement;
        first_movt = find(raw_move_data>0,1);
        raw_move_data(first_movt) = NaN;
        
        hs.full_movt_data = [hs.full_movt_data raw_move_data];
        hs.video_beg_secs(avi_num) = avi_end_time-hs.movt_data_length/hs.fps;
    catch
        disp(['Unable to load file ',num2str(avi_num)])
    end
end
catch
    textprogressbar('crashed');
end
textprogressbar('done');

if length(hs.full_movt_data) ~= length(hs.full_vid_fr_times)
    disp('num frames and times do not match')
    keyboard
end

hs.abs_file_starttime = datenum(hs.mat_files(1).date)*(24*3600)- hs.filesave_mins*60; %File datenum start time in seconds

%%%
hs.epoch_num = hs.epoch_num_start - 1;
next_but
%%%

% Load the next window of data (epoch)
function next_but(~,~)
disp('next..')

% %save
% epoch_indx = (hs.epoch_num*hs.epoch_size+1):hs.epoch_num*hs.epoch_size+hs.epoch_size;
% hs.sleep_state_data(epoch_indx) = hs.brain_state;
% clock_time = clock;
% savename = ['SSData_',num2str(clock_time(2)),'_',num2str(clock_time(3))];
% assignin('base',savename,hs.sleep_state_data);

% Report the running accuracy of classification based on changed segments
if hs.first_run ~= 1
    hs.running_ep_num = hs.running_ep_num + 1;
    hs.num_inds_changd = hs.num_inds_changd + sum(hs.temp_BS~=hs.brain_state);
    run_acc = (1 - hs.num_inds_changd/(hs.running_ep_num*hs.viewer_mins*60*hs.upsamp_by/hs.window_sec));
    disp(['Running Acc. : ',num2str(run_acc*100),'%'])
end

hs.epoch_num = hs.epoch_num + 1;

hs.eeg_data = [];
hs.eeg_data = nan([1,round(hs.num_files_to_load*hs.trace_length)]);
hs.emg_data = [];
hs.emg_data = nan([1,round(hs.num_files_to_load*hs.trace_length)]);

first_file_num = (hs.epoch_num*hs.num_files_to_load+1);
last_file_num = (first_file_num + hs.num_files_to_load);

ind_st = (first_file_num-1)*(hs.trace_length)+1;
ind_end = ind_st + hs.num_files_to_load*(hs.trace_length) - 1;
hs.eeg_data = hs.full_eeg_data(ind_st:ind_end)-nanmean(hs.full_eeg_data(ind_st:ind_end)); % center around 0
hs.eeg_data(isnan(hs.eeg_data)) = nanmean(hs.eeg_data);
hs.emg_data = hs.full_emg_data(ind_st:ind_end)-nanmean(hs.full_emg_data(ind_st:ind_end)); % center around 0

% Collect the video fr times within this viewer window
start_sec = hs.epoch_num*hs.viewer_mins*60;
end_sec = start_sec + hs.viewer_mins*60;
epoch_fr_ind = (hs.full_vid_fr_times>=start_sec & hs.full_vid_fr_times<=end_sec);
hs.vid_fr_times = hs.full_vid_fr_times(epoch_fr_ind) - start_sec;

% Use interp function to shift movt data according to fr times
hs.temp_movt_data = hs.full_movt_data(epoch_fr_ind);
hs.movt_data = zeros(1,hs.viewer_mins*60*hs.fps);
sec_int = 1/hs.fps;
st_sec = 0;
end_sec = sec_int;
hs.movt_data = interp1(hs.vid_fr_times,hs.temp_movt_data,sec_int:sec_int:hs.viewer_mins*60);

% time variable used for plotting
hs.time = (hs.epoch_num*hs.viewer_mins*60 + ((1/hs.sample_rate):(1/hs.sample_rate ):hs.viewer_mins*60));
hs.time = hs.time./60;

analyze_and_plot % run analysis and plotting now that this view window's data is loaded
end

% Load the previous window of data (epoch)
function back_but(~,~)
disp('back..')

if hs.epoch_num > hs.epoch_num_start
    hs.epoch_num = hs.epoch_num - 1;
    hs.running_ep_num = hs.running_ep_num - 1;
else
    hs.epoch_num = hs.epoch_num_start;
end

%%%
hs.eeg_data = [];
hs.eeg_data = nan([1,round(hs.num_files_to_load*hs.trace_length)]);
hs.emg_data = [];
hs.emg_data = nan([1,round(hs.num_files_to_load*hs.trace_length)]);

first_file_num = (hs.epoch_num*hs.num_files_to_load+1);
last_file_num = (first_file_num + hs.num_files_to_load);

ind_st = (first_file_num-1)*(hs.trace_length)+1;
ind_end = ind_st + hs.num_files_to_load*(hs.trace_length) - 1;
hs.eeg_data = hs.full_eeg_data(ind_st:ind_end)-nanmean(hs.full_eeg_data(ind_st:ind_end)); % center around 0
hs.eeg_data(isnan(hs.eeg_data)) = nanmean(hs.eeg_data);
hs.emg_data = hs.full_emg_data(ind_st:ind_end)-nanmean(hs.full_emg_data(ind_st:ind_end)); % center around 0

% Collect the video fr times within this viewer window
start_sec = hs.epoch_num*hs.viewer_mins*60;
end_sec = start_sec + hs.viewer_mins*60;
epoch_fr_ind = (hs.full_vid_fr_times>=start_sec & hs.full_vid_fr_times<=end_sec);
hs.vid_fr_times = hs.full_vid_fr_times(epoch_fr_ind) - start_sec;

% Use interp function to shift movt data according to fr times
hs.temp_movt_data = hs.full_movt_data(epoch_fr_ind);
hs.movt_data = zeros(1,hs.viewer_mins*60*hs.fps);
sec_int = 1/hs.fps;
st_sec = 0;
end_sec = sec_int;
hs.movt_data = interp1(hs.vid_fr_times,hs.temp_movt_data,sec_int:sec_int:hs.viewer_mins*60);

% time variable used for plotting
hs.time = (hs.epoch_num*hs.viewer_mins*60 + ((1/hs.sample_rate):(1/hs.sample_rate ):hs.viewer_mins*60));
hs.time = hs.time./60;

if hs.time(1) < 0 % prevent going 'back' if at first view window of data
    hs.time = ((1/hs.sample_rate ):(1/hs.sample_rate ):length(hs.eeg_data)/hs.sample_rate)/60;
end

analyze_and_plot % run analysis and plotting now that this view window's data is loaded
end

% core function that processes the loaded data and plots figures to GUI
function analyze_and_plot(~,~)
% try

% get voltage threshold from 'box' in GUI
voltthresh = get(hs.voltthresh_box,'String');
hs.volt_thresh = str2double(voltthresh);

bin_sec = hs.window_sec;
bin_ind_length = bin_sec*hs.sample_rate;

bin_stds = [];
bin_st = 1;
bin_end = bin_ind_length;

% remove eeg data beyond the voltage thresh
hs.eeg_data(hs.eeg_data > hs.volt_thresh |...
    hs.eeg_data < -hs.volt_thresh) = 0;

% collect eeg_max_min for each window_sec
hs.eeg_max_min = [];
for bin = 1:(size(hs.eeg_data,2)/bin_ind_length)
    this_data = hs.eeg_data(bin_st:bin_end);
    
%     hs.eeg_max_min(bin) = prctile(this_data,99)-prctile(this_data,1);
    hs.eeg_max_min(bin) = max(this_data)-min(this_data);
    bin_st = bin_st + bin_ind_length;
    bin_end = bin_end + bin_ind_length;
end
hs.eeg_data(isnan(hs.eeg_data)) = 0;

% perform fft on eeg data with window_sec time bin resolution,
% non-overlapping
wind_size = hs.window_sec*hs.sample_rate;
[s, f, t] = spectrogram(hs.eeg_data, wind_size, 0, [], hs.sample_rate , 'yaxis');

hs.spec_time = t;

powers = real(s).^2 + imag(s).^2; %convert power into positive real values
powers(powers==0) = NaN;
 
hs.max_freq = 81; % set a max and min frequency for power analysis
min_freq = 0.5;
select_freqs = find(f<=hs.max_freq & f>min_freq);
powers = (powers(select_freqs,:));

freq_ticks = f(select_freqs);
hs.freq_ticks = freq_ticks;

%%%
% Plot the raw EEG trace at the top of viewer window, with optional
% 'skip_ind' to lower plotting load
skip_ind = 1;
hs.raw_trace_ax = subplot(11,1,1:2,'Parent',hs.fig);
plot(hs.raw_trace_ax, hs.time(1:skip_ind:end)./60,hs.eeg_data(1:skip_ind:end),'k')
box(hs.raw_trace_ax,'off')
ylabel(hs.raw_trace_ax,'EEG (V)')
try ylim(hs.raw_trace_ax, [prctile(hs.eeg_data,0.0005), prctile(hs.eeg_data,99.9995)]); catch; end
xlim(hs.raw_trace_ax,[hs.time(1)./60 hs.time(end)./60])

power_to_plot = log10(powers); % take the log10 of powers for better visualization
blur_power = imgaussfilt(power_to_plot,0.5);
hs.spect_ax = subplot(11,1,3:5,'Parent',hs.fig);
imagesc(t./hs.time_conv,freq_ticks,blur_power,'Parent',hs.spect_ax) % plot the spectrogram results

% change contrast of color axis
low_perc = 62;
high_perc = 98;
caxis(hs.spect_ax,[prctile(power_to_plot(:),low_perc) prctile(power_to_plot(:),high_perc)])

set(gca,'YDir','normal')
hs.spect_ymax = 35;
ylim(hs.spect_ax,[0 hs.spect_ymax])
ylabel(hs.spect_ax,'EEG Spect. (Freq. Hz)')

% define different power bands
d_powers = powers(freq_ticks<=4 & freq_ticks>0.5,:);
t_powers = powers(freq_ticks<=8 & freq_ticks>5,:);
b_powers = powers(freq_ticks<=34.5 & freq_ticks>20,:);
g_powers = powers(freq_ticks<=80 & freq_ticks>30,:);

total_power = nanmean(powers);

% create minorly smoothed power means/medians from power bands
d_power = nanmean(d_powers);
d_power = smooth(d_power,3);
d_power = d_power';

b_power = nanmedian(b_powers);
b_power = smooth(b_power,4);
b_power = b_power';

t_power = nanmean(t_powers);
t_power = smooth(t_power,3);
t_power = t_power';

g_power = nanmedian(g_powers);
g_power = smooth(g_power,4);
g_power = g_power';

% create power ratios
d_to_b_ratio = d_power./b_power;
t_to_d_ratio = t_power./d_power;
g_to_d_ratio = g_power./total_power;

% perform minor smoothing on ratios and remove negative values
d_to_b_ratio = smooth(d_to_b_ratio,2);
t_to_d_ratio = smooth(t_to_d_ratio,3);
g_to_d_ratio = smooth(g_to_d_ratio,3);
d_to_b_ratio(d_to_b_ratio<0) = 0;
t_to_d_ratio(t_to_d_ratio<0) = 0;
g_to_d_ratio(g_to_d_ratio<0) = 0;

% prepare a movement "average" from weighted combination of avg and median
hs.move_average = zeros(1,hs.epoch_size);
start_ind = 1;
end_ind = start_ind + hs.fps*hs.window_sec;
coef_bias_mean = 0.25;
wind_avg = hs.fps*hs.window_sec;
movt_mean = nan(1,hs.epoch_size);
movt_mean(1:floor(length(hs.movt_data)/wind_avg)) = slide_window_avg(hs.movt_data,wind_avg);

movt_med = nan(1,hs.epoch_size);
movt_med(1:floor(length(hs.movt_data)/wind_avg)) = ...
    nanmedian(reshape(hs.movt_data(1:wind_avg * floor(numel(hs.movt_data) / wind_avg)), wind_avg, []), 1);

hs.move_average = (coef_bias_mean*movt_mean + movt_med)/2;
hs.sm_mov_avg = smooth(hs.move_average,3);


hs.spect_ax = subplot(11,1,3:5,'Parent',hs.fig);
hold(hs.spect_ax,'on')

% plot delta ratio values
min_d_toplot = 10; % arbitrary lower minimum for plotting given use of y log scale
d_toplot = d_to_b_ratio; 
d_toplot(d_toplot==0) = min_d_toplot;
highd_times = t(d_toplot > hs.nrem_thresh)./hs.time_conv; % d ratio pts above thresh
highd_pts = d_toplot(d_toplot > hs.nrem_thresh);
hs.delta_ax = subplot(11,1,6:7,'Parent',hs.fig); % define subplot axis
cla(hs.delta_ax,'reset')
hold(hs.delta_ax,'on')
box(hs.delta_ax,'off')
plot(hs.delta_ax,t./hs.time_conv, d_toplot,'Color',[0.2 0.7 1],'LineWidth',1.5)
plot(hs.delta_ax,highd_times,highd_pts,'.','Color',[0.2 0.7 1],'MarkerSize',22)
xlim(hs.delta_ax,[0 hs.viewer_mins])
rd = refline(0,hs.nrem_thresh);
rd.LineStyle = '--';
rd.Color = [0.2 0.7 1];
try ylim(hs.delta_ax,[min_d_toplot max(d_toplot)]); catch; end
set(hs.delta_ax,'yscale','log')
ylabel(hs.delta_ax,'Delta/beta')

% plot theta ratio same as delta
min_t_toplot = 0.15;
t_toplot = t_to_d_ratio;
t_toplot(t_toplot==0) = min_t_toplot;
hight_times = t(t_toplot > hs.rem_thresh)./hs.time_conv;
hight_pts = t_toplot(t_toplot > hs.rem_thresh);
hs.theta_ax = subplot(11,1,8:9,'Parent',hs.fig);
cla(hs.theta_ax,'reset')
hold(hs.theta_ax,'on')
plot(hs.theta_ax,t./hs.time_conv, t_toplot,'Color',[1 0.25 1],'LineWidth',1.3)
plot(hs.theta_ax,hight_times,hight_pts,'.','Color',[0.9 0.15 0.9],'MarkerSize',22)
rt = refline(0,hs.rem_thresh);
rt.LineStyle = '--';
rt.Color = [0.9 0.15 0.9];
hold(hs.theta_ax,'off')
xlim(hs.theta_ax,[0 hs.viewer_mins])
try ylim(hs.theta_ax,[min_t_toplot, max(t_to_d_ratio)]); catch; end
set(hs.theta_ax,'yscale','log')
ylabel(hs.theta_ax,'Theta/delta')


% Prepare EMG values
emg_data(isnan(emg_data)) = nanmean(emg_data);
emg_data = filtfilt(hs.emg_filt,emg_data);

avg_int = round(hs.sample_rate*hs.state_sec);
std_devs = 4;
emg_std = nanstd(emg_data); 

abs_emg = abs(emg_data);
abs_emg(abs_emg<std_devs*emg_std) = NaN; % remove abs_emg vals below std thresh
emg_avg = slide_window_avg(abs_emg,avg_int); % apply sliding window (non overlapping)
emg_avg = emg_avg./emg_std - std_devs; % normalize by std
emg_avg(isnan(emg_avg)) = 0;

hs.emg_average = slide_window_avg(emg_avg,round(hs.window_sec/hs.state_sec)); % emg avg vector for thresh-based algorithm

% prepare movt data
avg_int = round(hs.fps/hs.state_sec);
movt_data_forML = slide_window_avg(hs.movt_data,avg_int); % average movt data to stanard state_sec size
movt_data_forML(isnan(movt_data_forML)) = 0;

% plot movement/emg values
max_movt_y = 15; % arbitrary max for visualization only
max_emg_y = 10; % arbitrary max for visualization only
gamma_toplot = g_to_d_ratio;
emg_toplot = emg_avg;

hs.move_ax = subplot(11,1,10:11,'Parent',hs.fig);
cla(hs.move_ax,'reset')
hold(hs.move_ax, 'on')
yyaxis(hs.move_ax,'left') % prime left y axis
hs.move_ax.YColor = 'k';
ylim(hs.move_ax,[0 max_movt_y]);
move_time = 1/hs.state_sec:1/hs.state_sec:length(movt_data_forML)/hs.state_sec;
plot(hs.move_ax, (move_time/60+hs.time(1))/hs.time_conv, smooth(movt_data_forML,4),'k')
plot(hs.move_ax, (t/60+hs.time(1))/hs.time_conv, hs.move_average,'go','MarkerSize',5)
rm = refline(hs.move_ax,0,hs.movt_thresh);
rm.LineStyle = '--';
rm.Color = 'k';
if hs.EMG_toplot % choose iether emg or gamma for label
    ylabel(hs.move_ax,'Movt. (black)/ EMG (red)')
    ylim(hs.move_ax,[0 max_emg_y]);
else
    ylabel(hs.move_ax,'Movt. (black)/ gamma sc. (red)')
end

yyaxis(hs.move_ax,'right') % prime right y axis
hs.move_ax.YColor = hs.emg_color(1:3);
if hs.EMG_toplot % plot either emg or gamma on right y axis
    plot(hs.move_ax, (move_time/60+hs.time(1))/hs.time_conv, smooth(emg_toplot,2),'Color',hs.emg_color)
else
    plot(hs.move_ax, (t/60+hs.time(1))/hs.time_conv, smooth(gamma_toplot,2),'Color',hs.emg_color)
end
rm = refline(hs.move_ax,0,hs.emg_thresh);
rm.LineStyle = '--';
rm.Color = hs.emg_color;
hold(hs.move_ax, 'off')
box(hs.move_ax,'off')
xlim(hs.move_ax,[hs.time(1)/hs.time_conv (hs.viewer_mins+hs.time(1))/hs.time_conv]);
xlabel('Hours since start')

%%%%%%%%%%%%%%%%%%%%
%%%prepare_ML_vars%%
%%%%%%%%%%%%%%%%%%%%

% define index values for viewer window (or 'epoch') to pull from
% previously logged values or add new ones
epoch_indx_up = (hs.epoch_num*hs.epoch_size*hs.upsamp_by+1):...
        hs.epoch_num*hs.epoch_size*hs.upsamp_by+hs.epoch_size*hs.upsamp_by;
    
% prepare some machine learning features with consistent number of bins
hs.movt_data_forML(epoch_indx_up) = movt_data_forML;
hs.emg_avg(epoch_indx_up) = emg_avg;
hs.d_to_b_ratio(epoch_indx_up) = repelem(d_to_b_ratio,hs.upsamp_by);
hs.t_to_d_ratio(epoch_indx_up) = repelem(t_to_d_ratio,hs.upsamp_by);
hs.g_to_d_ratio(epoch_indx_up) = repelem(g_to_d_ratio,hs.upsamp_by);

% if viewer window ('epoch') hasn't been previously analyzed, log all of
% the ML features in the matrix containing all ML variables
sh_ind = hs.num_prev_states;
if hs.have_analyzed(hs.epoch_num + 1) == 0
if hs.epoch_num ~= 0
    try
        hs.ML_indx = (hs.epoch_num*hs.epoch_size*hs.upsamp_by+1-sh_ind*hs.upsamp_by):...
            hs.epoch_num*hs.epoch_size*hs.upsamp_by+hs.epoch_size*hs.upsamp_by-sh_ind*hs.upsamp_by;

        hs.ML_delta(hs.ML_indx) = repelem(d_to_b_ratio,hs.upsamp_by);
        hs.ML_delta_b1(hs.ML_indx) = hs.d_to_b_ratio(...
            end-hs.epoch_size*hs.upsamp_by-hs.upsamp_by+1:end-hs.upsamp_by);
        hs.ML_delta_b2(hs.ML_indx) = hs.d_to_b_ratio(...
            end-hs.epoch_size*hs.upsamp_by-(sh_ind-2)*hs.upsamp_by+1:end-(sh_ind-2)*hs.upsamp_by);
        hs.ML_delta_b3(hs.ML_indx) = hs.d_to_b_ratio(...
            end-hs.epoch_size*hs.upsamp_by-(sh_ind-1)*hs.upsamp_by+1:end-(sh_ind-1)*hs.upsamp_by);
        
        hs.ML_theta(hs.ML_indx) = repelem(t_to_d_ratio,hs.upsamp_by);
        hs.ML_theta_b1(hs.ML_indx) = hs.t_to_d_ratio(...
            end-hs.epoch_size*hs.upsamp_by-hs.upsamp_by+1:end-hs.upsamp_by);
        hs.ML_theta_b2(hs.ML_indx) = hs.t_to_d_ratio(...
            end-hs.epoch_size*hs.upsamp_by-(sh_ind-2)*hs.upsamp_by+1:end-(sh_ind-2)*hs.upsamp_by);
         
        hs.ML_gamma(hs.ML_indx) = repelem(g_to_d_ratio,hs.upsamp_by);
        hs.ML_gamma_b1(hs.ML_indx) = hs.g_to_d_ratio(...
            end-hs.epoch_size*hs.upsamp_by-hs.upsamp_by+1:end-hs.upsamp_by);
        hs.ML_gamma_b2(hs.ML_indx) = hs.g_to_d_ratio(...
            end-hs.epoch_size*hs.upsamp_by-(sh_ind-2)*hs.upsamp_by+1:end-(sh_ind-2)*hs.upsamp_by);

        hs.ML_movt(hs.ML_indx) = movt_data_forML;
        hs.ML_movt_b1(hs.ML_indx) = hs.movt_data_forML(...
            end-hs.epoch_size*hs.upsamp_by-hs.upsamp_by+1:end-hs.upsamp_by);
        hs.ML_movt_b2(hs.ML_indx) = hs.movt_data_forML(...
            end-hs.epoch_size*hs.upsamp_by-(sh_ind-2)*hs.upsamp_by+1:end-(sh_ind-2)*hs.upsamp_by);

        hs.d_power(hs.ML_indx) = repelem(d_power,hs.upsamp_by);
        hs.t_power(hs.ML_indx) = repelem(t_power,hs.upsamp_by);
        hs.sm_mov_avg_forML(hs.ML_indx) = repelem(hs.sm_mov_avg,hs.upsamp_by);
        hs.move_average_forML(hs.ML_indx) = repelem(hs.move_average,hs.upsamp_by);
        hs.eeg_max_min_forML(hs.ML_indx) = repelem(hs.eeg_max_min,hs.upsamp_by);
        
        hs.ML_EMG(hs.ML_indx) = emg_avg;
        hs.ML_EMG_b1(hs.ML_indx) = hs.emg_avg(...
            end-hs.epoch_size*hs.upsamp_by-hs.upsamp_by+1:end-hs.upsamp_by);
        hs.ML_EMG_b2(hs.ML_indx) = hs.emg_avg(...
            end-hs.epoch_size*hs.upsamp_by-(sh_ind-2)*hs.upsamp_by+1:end-(sh_ind-2)*hs.upsamp_by);
        
    catch
       disp('ML var loading err')
       keyboard
    end
else
    hs.ML_indx = 1:hs.epoch_num*hs.epoch_size*hs.upsamp_by...
        +hs.epoch_size*hs.upsamp_by-sh_ind*hs.upsamp_by;

    hs.ML_delta(hs.ML_indx) = repelem(d_to_b_ratio(sh_ind+1:end),hs.upsamp_by);
    hs.ML_delta_b1(hs.ML_indx) = repelem(d_to_b_ratio(...
         end-length(d_to_b_ratio)+sh_ind:end-sh_ind+2),hs.upsamp_by);
    hs.ML_delta_b2(hs.ML_indx) = repelem(d_to_b_ratio(...
         end-length(d_to_b_ratio)+sh_ind-1:end-sh_ind+1),hs.upsamp_by);
    hs.ML_delta_b3(hs.ML_indx) = repelem(d_to_b_ratio(...
         end-length(d_to_b_ratio)+sh_ind-2:end-sh_ind),hs.upsamp_by);
     
    hs.ML_theta(hs.ML_indx) = repelem(t_to_d_ratio(sh_ind+1:end),hs.upsamp_by);
    hs.ML_theta_b1(hs.ML_indx) = repelem(t_to_d_ratio(...
         end-length(d_to_b_ratio)+sh_ind:end-1),hs.upsamp_by);
    hs.ML_theta_b2(hs.ML_indx) = repelem(t_to_d_ratio(...
         end-length(d_to_b_ratio)+sh_ind-1:end-2),hs.upsamp_by);

    hs.ML_gamma(hs.ML_indx) = repelem(g_to_d_ratio(sh_ind+1:end),hs.upsamp_by);
    hs.ML_gamma_b1(hs.ML_indx) = repelem(g_to_d_ratio(...
         end-length(d_to_b_ratio)+sh_ind:end-1),hs.upsamp_by);
    hs.ML_gamma_b2(hs.ML_indx) = repelem(g_to_d_ratio(...
         end-length(d_to_b_ratio)+sh_ind-1:end-2),hs.upsamp_by);

    hs.ML_movt(hs.ML_indx) = movt_data_forML(sh_ind*hs.upsamp_by+1:end);
    hs.ML_movt_b1(hs.ML_indx) = hs.movt_data_forML(...
         end-length(movt_data_forML)+(sh_ind-1)*hs.upsamp_by+1:end-(sh_ind-2)*hs.upsamp_by);
    hs.ML_movt_b2(hs.ML_indx) = hs.movt_data_forML(...
         end-length(movt_data_forML)+(sh_ind-2)*hs.upsamp_by+1:end-(sh_ind-1)*hs.upsamp_by);

    hs.d_power(hs.ML_indx) = repelem(d_power(sh_ind+1:end),hs.upsamp_by);
    hs.t_power(hs.ML_indx) = repelem(t_power(sh_ind+1:end),hs.upsamp_by);
    hs.sm_mov_avg_forML(hs.ML_indx) = repelem(hs.sm_mov_avg(sh_ind+1:end),hs.upsamp_by);
    hs.move_average_forML(hs.ML_indx) = repelem(hs.move_average(sh_ind+1:end),hs.upsamp_by);
    hs.eeg_max_min_forML(hs.ML_indx) = repelem(hs.eeg_max_min(sh_ind+1:end),hs.upsamp_by);
     
    hs.ML_EMG(hs.ML_indx) = emg_avg(sh_ind*hs.upsamp_by+1:end);
    hs.ML_EMG_b1(hs.ML_indx) = hs.emg_avg(...
         end-length(emg_avg)+(sh_ind-1)*hs.upsamp_by+1:end-(sh_ind-2)*hs.upsamp_by);
    hs.ML_EMG_b2(hs.ML_indx) = hs.emg_avg(...
         end-length(emg_avg)+(sh_ind-2)*hs.upsamp_by+1:end-(sh_ind-1)*hs.upsamp_by);    
       
end
end

hs.ML_input_vars = [hs.ML_delta;hs.ML_delta_b1;hs.ML_delta_b2;hs.ML_delta_b3;
                    hs.ML_theta;hs.ML_theta_b1;hs.ML_theta_b2;
                    hs.ML_gamma;hs.ML_gamma_b1;hs.ML_gamma_b2;
                    hs.ML_movt;hs.ML_movt_b1;hs.ML_movt_b2;
                    hs.d_power;hs.t_power;hs.sm_mov_avg_forML;hs.move_average_forML;hs.eeg_max_min_forML;
                    hs.ML_EMG;hs.ML_EMG_b1;hs.ML_EMG_b2];
                                
%%%%%%%%%%%%%%%%%%%%%% 
%%%%Classify states%%%
%%%%%%%%%%%%%%%%%%%%%%

% run classification (RF ML or treshold-based)
if hs.have_analyzed(hs.epoch_num + 1) == 0
    if hs.ML_onoff == 1 && hs.epoch_num ~= 0
        ML_indx = (hs.epoch_num*hs.epoch_size*hs.upsamp_by+1-hs.num_prev_states*hs.upsamp_by):...
            hs.epoch_num*hs.epoch_size*hs.upsamp_by+hs.epoch_size*hs.upsamp_by-hs.num_prev_states*hs.upsamp_by;
        feature_array = hs.ML_input_vars(:,ML_indx)';
        ML_classifier(feature_array,hs.anim) % run ML classifier
        hs.brain_state = remove_small_states(hs.brain_state);
        temp_conf = hs.confidence;
        
        %run classifier again
        f_sh_ind = hs.num_2nd_for_states;
        b_sh_ind = hs.num_prev_states;
        hs.SSD_b1 = hs.brain_state((b_sh_ind-1)*hs.upsamp_by+1:end-(b_sh_ind-2)*hs.upsamp_by-f_sh_ind*hs.upsamp_by);
        hs.SSD_b2 = hs.brain_state((b_sh_ind-2)*hs.upsamp_by+1:end-(b_sh_ind-1)*hs.upsamp_by-f_sh_ind*hs.upsamp_by);
        hs.SSD_b3 = hs.brain_state((b_sh_ind-3)*hs.upsamp_by+1:end-(b_sh_ind-0)*hs.upsamp_by-f_sh_ind*hs.upsamp_by);

        hs.SSD_f1 = hs.brain_state((b_sh_ind+f_sh_ind-1)*hs.upsamp_by+1:end-(f_sh_ind-1)*hs.upsamp_by);
        hs.SSD_f2 = hs.brain_state((b_sh_ind+f_sh_ind)*hs.upsamp_by+1:end-(f_sh_ind-2)*hs.upsamp_by);

        ML_indx = (hs.epoch_num*hs.epoch_size*hs.upsamp_by+1-hs.num_prev_states*hs.upsamp_by):...
                    hs.epoch_num*hs.epoch_size*hs.upsamp_by+hs.epoch_size*hs.upsamp_by-hs.num_prev_states*hs.upsamp_by;
        ML_indx = ML_indx(b_sh_ind*hs.upsamp_by+1:end-f_sh_ind*hs.upsamp_by);
        ML_input_vars_short = hs.ML_input_vars(:,ML_indx);
        feature_array = [ML_input_vars_short;
            hs.SSD_b1; hs.SSD_b2; hs.SSD_b3;
            hs.SSD_f1; hs.SSD_f2]';        
        
        temp_bs = hs.brain_state;
        ML_classifier(feature_array,hs.anim_ssd)

        temp_bs(b_sh_ind*hs.upsamp_by+1:end-f_sh_ind*hs.upsamp_by) = hs.brain_state;
        hs.brain_state = remove_small_states(temp_bs); % remove states that are isolated or too small
    else
        % if ML mode isn't on apply tresh-based algorithm
        Thresh_State_Algo(t,d_to_b_ratio,t_to_d_ratio,hs.emg_average,hs.move_average)
    end
else
    hs.brain_state = hs.sleep_state_data((hs.epoch_num*hs.epoch_size*hs.upsamp_by+1):...
        hs.epoch_num*hs.epoch_size*hs.upsamp_by+hs.epoch_size*hs.upsamp_by);
end

% Add SSD classification to ML_vars
epoch_indx_up = (hs.epoch_num*hs.epoch_size*hs.upsamp_by+1):...
        hs.epoch_num*hs.epoch_size*hs.upsamp_by+hs.epoch_size*hs.upsamp_by;
hs.sleep_state_data(epoch_indx_up) = hs.brain_state;

% if ML was used, plot the confidences from the ML output
if hs.ML_onoff
    if hs.have_analyzed(hs.epoch_num + 1) == 0
        f_sh_ind = hs.num_2nd_for_states;
        b_sh_ind = hs.num_prev_states;
        temp_conf(b_sh_ind*hs.upsamp_by+1:end-f_sh_ind*hs.upsamp_by,:) = hs.confidence;
        hs.confidence = temp_conf;
    else
        hs.confidence = hs.full_confidence(epoch_indx_up,:);
    end

    %Plot confidences
    y_lim = [22.5 27.5];
    conf_toplot = hs.confidence*(y_lim(2)-y_lim(1))+y_lim(1);
    wake_conf = conf_toplot(:,1);
    nrem_conf = conf_toplot(:,3);
    rem_conf = conf_toplot(:,2);

    rectangle(hs.spect_ax,'Position',[0, y_lim(1)-1,...
        16, y_lim(2)-y_lim(1)+2],'FaceColor',[1 1 1 0.6],...
        'linestyle','none');
    plot(hs.spect_ax,move_time./hs.time_conv,wake_conf,'g','LineWidth',2);
    plot(hs.spect_ax,move_time./hs.time_conv,nrem_conf,'Color',[0.2 0.7 1],'LineWidth',2);
    plot(hs.spect_ax,move_time./hs.time_conv,rem_conf,'m','LineWidth',2);

    hs.full_confidence(epoch_indx_up,:) = hs.confidence;
end

% log the classified states in the ML labels vector
sh_ind = hs.num_prev_states;
if hs.epoch_num ~= 0
        hs.ML_indx = (hs.epoch_num*hs.epoch_size*hs.upsamp_by+1-sh_ind*hs.upsamp_by):...
            hs.epoch_num*hs.epoch_size*hs.upsamp_by+hs.epoch_size*hs.upsamp_by-sh_ind*hs.upsamp_by;

        hs.ML_SSD(hs.ML_indx) = hs.sleep_state_data(hs.ML_indx);
else
    hs.ML_indx = 1:hs.epoch_num*hs.epoch_size*hs.upsamp_by...
        +hs.epoch_size*hs.upsamp_by-sh_ind*hs.upsamp_by;

    hs.ML_SSD(hs.ML_indx) = hs.brain_state(sh_ind*hs.upsamp_by+1:end);
end

%%%
% Plot the classified states in the GUI
hs.temp_BS = hs.brain_state;
states_toplot = hs.brain_state;
plot_state_rects(states_toplot)
%%%

% save variables to workspace
clock_time = clock;
savename = ['SSData_',num2str(clock_time(2)),'_',num2str(clock_time(3))];
assignin('base',savename,hs.sleep_state_data);

assignin('base','ML_input_vars',hs.ML_input_vars)
assignin('base','ML_confidences',hs.full_confidence)

metadata.starttime = hs.abs_file_starttime + hs.epoch_num_start*hs.viewer_mins*60;
metadata.darkhour = hs.darkhour;
assignin('base','metadata',metadata)

hs.first_run = 0;
hs.have_analyzed(hs.epoch_num + 1) = 1;

drawnow

end

% Function used apply hotkeys within GUI window
function key_catcher(~,eventdata)

    try
    button = eventdata.Character; % retrieve pressed button identity
    hold(hs.spect_ax,'on')
    state_inds = 1:length(hs.brain_state);
    
    % select state decisions for revision/review
    if strcmp(button, ' ') || strcmp(button, 'o') % if button pressed matches these options
        
        [x,y] = ginput;
        
        if length(x) == 2
            if x(2) > x(1)
                hs.ind1 = (x(1)*(60/(hs.window_sec/hs.upsamp_by)));
                hs.ind2 = (x(2)*(60/(hs.window_sec/hs.upsamp_by)));
            else
                hs.ind2 = (x(1)*(60/(hs.window_sec/hs.upsamp_by)));
                hs.ind1 = (x(2)*(60/(hs.window_sec/hs.upsamp_by)));
            end
            
            if hs.ind1 < 0
                hs.ind1 = 0;
            end
            if hs.ind2 > hs.epoch_size*hs.upsamp_by
                hs.ind2 = hs.epoch_size*hs.upsamp_by;
            end        
            
            hs.inds_to_change = state_inds-0.5 > hs.ind1 & state_inds-0.5 < hs.ind2;
            
            x(x<0) = (hs.window_sec/hs.upsamp_by)/(2*60);
            hs.plotted_coords = scatter(hs.spect_ax,x,y,'rx','LineWidth',3);

            hs.coords_chosen = 1;
        end
    end
    
    if strcmp(button, '=') % run next_but function to continue
        next_but
    end

    if strcmp(button, '-')
        back_but
    end
    
    % change selected states to wake
    if (strcmp(button, 'w') || strcmp(button, 'p')) && hs.coords_chosen == 1
        hs.brain_state(hs.inds_to_change) = 1;
        hs.coords_chosen = 0;
    end
    
    % change selected states to nrem
    if (strcmp(button, 'e') || strcmp(button, '['))  && hs.coords_chosen == 1
        hs.brain_state(hs.inds_to_change) = 3;
        hs.coords_chosen = 0;
    end
    
    % change selected states to rem
    if (strcmp(button, 'r') || strcmp(button, ']'))  && hs.coords_chosen == 1
        hs.brain_state(hs.inds_to_change) = 2;
        hs.coords_chosen = 0;        
    end
    
    % retrieve video for selected time
    if (strcmp(button, 'k') || strcmp(button, 'v'))  && hs.coords_chosen == 1
        start_sec = hs.ind1*(hs.window_sec/hs.upsamp_by)+hs.epoch_num*hs.viewer_mins*60;
        end_sec = hs.ind2*(hs.window_sec/hs.upsamp_by)+hs.epoch_num*hs.viewer_mins*60;
        display_video(start_sec, end_sec)
    end
    
    if hs.coords_chosen == 0 && isfield(hs,'plotted_coords')
        delete(hs.plotted_coords)
    end
    
    epoch_indx_up = (hs.epoch_num*hs.epoch_size*hs.upsamp_by+1):...
            hs.epoch_num*hs.epoch_size*hs.upsamp_by+hs.epoch_size*hs.upsamp_by;
    hs.sleep_state_data(epoch_indx_up) = hs.brain_state;

    plot_state_rects(hs.brain_state)
    
    catch
        disp('Button err...')
    end
end

%%%

% Reload GUI window and reanalyze viewer window
function reload_next_but(~,~)

    voltthresh = get(hs.voltthresh_box,'String');
    hs.volt_thresh = str2double(voltthresh);    
    nrem_thresh = get(hs.dthresh_box,'String');
    hs.nrem_thresh = str2double(nrem_thresh);
    rem_thresh = get(hs.tthresh_box,'String');
    hs.rem_thresh = str2double(rem_thresh);
    movt_thresh = get(hs.movtthresh_box,'String');
    hs.movt_thresh = str2double(movt_thresh);
    emg_thresh = get(hs.emgthresh_box,'String');
    hs.emg_thresh = str2double(emg_thresh);
    
    close(hs.fig)
    create_window
    hs.have_analyzed(hs.epoch_num + 1) = 0;
    hs.epoch_num = hs.epoch_num - 1;
    next_but

end

%%%

% auto continue and analyze data windows
function auto_next_but(~,~)
    hs.auto_onoff = hs.auto_but.Value;

    if hs.auto_onoff
        set(hs.auto_but,'String','Auto On')
    else
        set(hs.auto_but,'String','Auto Off')
    end
    
    while hs.auto_onoff
        pause(0.5)
        next_but
    end
    
end

%%%

% toggle between ML On and Off
function switch_mode_but(~,~)

    hs.running_ep_num = 0;
    hs.num_inds_changd = 0;
    
    hs.ML_onoff = hs.switch_mode_but.Value;

    if hs.ML_onoff
        set(hs.switch_mode_but,'String','ML On')
    else
        set(hs.switch_mode_but,'String','ML Off')
    end
end

%%%

% Train the ML model using the ML_vars features and SSD labels
function trainML_but(~,~)
    
    % remove previously saved temp model variables 
    disp('Clearing temp dir...')
    rmdir(hs.ml_model_path,'s')
    mkdir(hs.ml_model_path)
    
    % run the ML training by calling the matlab python wrapper
    train_data = [hs.ML_input_vars; hs.sleep_state_data(hs.num_prev_states*hs.upsamp_by+1:end)]';
    [hs.pyMdl_path,hs.OOBerr,hs.saveMdl_path] = ...
        AutoVidCode_trainModel_PyWrap(train_data,hs.anim,hs.mlmode,hs.nTrees);
    
    disp(['Train Acc. is: ',num2str(100*hs.OOBerr)])
    
    % collect classified states and create new features out of them for
    % second training round
    f_sh_ind = hs.num_2nd_for_states;
    b_sh_ind = hs.num_prev_states;
    hs.SSD = hs.sleep_state_data(b_sh_ind*hs.upsamp_by+1:end-f_sh_ind*hs.upsamp_by);
    hs.SSD_b1 = hs.sleep_state_data((b_sh_ind-1)*hs.upsamp_by+1:end-(b_sh_ind-2)*hs.upsamp_by-f_sh_ind*hs.upsamp_by);
    hs.SSD_b2 = hs.sleep_state_data((b_sh_ind-2)*hs.upsamp_by+1:end-(b_sh_ind-1)*hs.upsamp_by-f_sh_ind*hs.upsamp_by);
    hs.SSD_b3 = hs.sleep_state_data((b_sh_ind-3)*hs.upsamp_by+1:end-(b_sh_ind-0)*hs.upsamp_by-f_sh_ind*hs.upsamp_by);
    
    hs.SSD_f1 = hs.sleep_state_data((b_sh_ind+f_sh_ind-1)*hs.upsamp_by+1:end-(f_sh_ind-1)*hs.upsamp_by);
    hs.SSD_f2 = hs.sleep_state_data((b_sh_ind+f_sh_ind)*hs.upsamp_by+1:end-(f_sh_ind-2)*hs.upsamp_by);
    
    ML_input_vars_short = hs.ML_input_vars(:,1:end-f_sh_ind*hs.upsamp_by);
    train_data = [ML_input_vars_short;
        hs.SSD_b1; hs.SSD_b2; hs.SSD_b3;
        hs.SSD_f1; hs.SSD_f2;
        hs.SSD]';
    
    % train once more with SSD features
    [hs.pyMdl_path,hs.OOBerr,hs.saveMdl_path] = ...
            AutoVidCode_trainModel_PyWrap(train_data,hs.anim_ssd,hs.mlmode,hs.nTrees);

    disp(['Train Acc. is: ',num2str(100*hs.OOBerr)])
end

%%%

%
function ML_classifier(feature_array,anim)
    
    % save to temp file
    homefolder = getenv('HOME');
    tempdir = [homefolder filesep 'Documents' filesep 'MATLAB' filesep 'tempPyVidcodeFiles'];
    if ~exist(tempdir, 'dir'), mkdir(tempdir); end
    tempfeat_code = 'temp_feature_input.mat';
    tempfeat = fullfile(tempdir,tempfeat_code);
    if exist(tempfeat,'file'), delete(tempfeat); end
    save(tempfeat,'feature_array');
    
    % call to python script for scoring
    ml_model_path = [hs.ml_model_path anim '_pyMdl.pkl'];
    py_score_call = [python_path ' ' pywrap_dir filesep ...
                'auto_video_score_Python.py -f ' tempfeat ' -sd ' tempdir ' -mdl ' ml_model_path];
    py_score_status = system(py_score_call);
    output_code = 'pyOut.mat';

    % now get labels and confidence scores from pyOut.mat in TEMP folder
    output_load = load([tempdir filesep output_code]);
    labels = output_load.labels;
    hs.confidence = output_load.confidence;
    hs.brain_state = labels;
        
end

%%%
function plot_state_rects(state_data)
    
    y_coord = 31;
    rect_height = 3;
    st_width = hs.window_sec/hs.upsamp_by; %in sec
    st_colors = [0 1 0; 0.9 0.15 0.9; 0.2 0.7 1];
    counter = 0;
    for st_num = 1:length(state_data)
        
        if st_num == length(state_data)
            st_width = hs.window_sec/hs.upsamp_by; %in sec
            x_coord = (st_num-counter)*st_width/60 - 1/60;
            st_width = st_width*(counter+1);
            rectangle(hs.spect_ax,'Position',[x_coord, y_coord,...
                st_width, rect_height],'FaceColor',st_colors(state_data(st_num-1),:),...
                'linestyle','none');
        elseif (st_num ~= 1) && (state_data(st_num) == state_data(st_num-1))
            counter = counter + 1;
        elseif (st_num == 1) || (state_data(st_num) ~= state_data(st_num-1))
            st_width = hs.window_sec/hs.upsamp_by; %in sec
            x_coord = (st_num-counter)*st_width/60 - 1/60;
            st_width = st_width*(counter+1);
            
            if (st_num ~= 1)
                rectangle(hs.spect_ax,'Position',[x_coord, y_coord,...
                    st_width, rect_height],'FaceColor',st_colors(state_data(st_num-1),:),...
                    'linestyle','none');
            else
                rectangle(hs.spect_ax,'Position',[x_coord, y_coord,...
                    st_width, rect_height],'FaceColor',st_colors(state_data(st_num),:),...
                    'linestyle','none');
            end
            
            counter = 1;
        end
        
    end

end

%%%
function [score_out] = remove_small_states(label_input)
        
    %clear small states
    autoscore_out = label_input';
    ticker = 1;
    while ticker == 1
        % eliminate states that only last for one 10-sec bin
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % copy the state assignment
        temp = autoscore_out(:,1);
        %find the differential. zero when a state is constant
        temp2 = diff(temp);
        % make it binary. all changes are 1s and constants are 0s
        temp2(temp2~=0) = 1;
        % index the changes
        temp3 = find(temp2 == 1);
        % find the amount of bins between changes
        temp4 = diff(temp3);
        % look for bins changes that are only seperated by 1
        temp5 = find(temp4<hs.min_state_size);
        % if you have staggered 1-bin states flanking one another, this will
        % need to run a couple of times to iteratively clean them out. move on
        % when finished:
        if isempty(temp5)
            ticker = 0;
        end
        % get the indices of the bins that have a 1 bin state
        kills = temp3(temp5);
        % overwrite the 1 bin state with the preceding state
        temp(kills+1) = temp(kills);
        % reassign temp into the autocat variable
        autoscore_out(:,1) = temp;
    end
    
    state_count = 0;
    new_score_out = autoscore_out;
    min_isolated_size = 20;
    
    for ii = 2:length(autoscore_out)
        if autoscore_out(ii-1) == autoscore_out(ii)
            state_count = state_count + 1;
        else
            small_yes = state_count <= min_isolated_size;
            not_wake = autoscore_out(ii-1) ~= 1;
            wake_flank_yes = (autoscore_out(ii-state_count-1) == 1)...
                && (autoscore_out(ii) == 1);
            if small_yes && not_wake && wake_flank_yes
                new_score_out(ii-state_count:ii-1) = 1;
            end
            
            state_count = 1;
        end
    end

    score_out = new_score_out';
end

%%%
function display_video(start_sec, end_sec)
    disp('Loading Video...')
    
    video1_ind = find(start_sec < hs.video_beg_secs,1)-1;
    video2_ind = find(end_sec < hs.video_beg_secs,1)-1;
    
    videos_toload = video1_ind:video2_ind;
    
    vidpaths = {};
    for ii = 1:length(videos_toload)
        vidpaths{ii} = [hs.vid_dir_path filesep hs.avi_filenames{videos_toload(ii)}];
    end
        
    first_frame = round(start_sec-hs.video_beg_secs(video1_ind))*hs.fps;
    last_frame = round(end_sec-hs.video_beg_secs(video2_ind))*hs.fps;
    
    full_frames = [];
    for ii = 1:length(videos_toload)
        c = VideoReader(vidpaths{ii});
        frame_total = c.Duration*c.FrameRate;
        
        if first_frame > frame_total
            first_frame = frame_total-8;
            disp('first frame out of range, fixing...')
        end
        if last_frame > frame_total
            disp('last frame out of range, fixing...')
            last_frame = frame_total;
        end
        
        if (ii == length(videos_toload)) && (ii == 1)
            num_frames = last_frame-first_frame+1;
            full_frames(:,:,:,end+1:end+num_frames) = read(c,[first_frame last_frame]);
        elseif ii == length(videos_toload)
            num_frames = last_frame;
            full_frames(:,:,:,end+1:end+num_frames) = read(c,[1 num_frames]);
        elseif ii == 1
            num_frames = frame_total - first_frame;
            this_frames = read(c,[first_frame frame_total]);
            full_frames(:,:,:,1:num_frames+1) = this_frames;
        else
            num_frames = frame_total;
            full_frames(:,:,:,end+1:num_frames) = read(c,[1 frame_total]);
        end
    end
    
    speed_upby = 8;
    implay(uint8(full_frames),hs.fps*speed_upby);
    delete(hs.plotted_coords);
end

%%%
function Thresh_State_Algo(state_times,d_ratio,t_ratio,emg_avg,move_avg)
    
t = state_times;
d_to_b_ratio = d_ratio;
t_to_d_ratio = t_ratio;
move_average = move_avg;
emg_average = emg_avg;

%brain state classifier
%1 = awake
%2 = rem
%3 = nrem

nrem_thresh = get(hs.dthresh_box,'String');
hs.nrem_thresh = str2double(nrem_thresh);
rem_thresh = get(hs.tthresh_box,'String');
hs.rem_thresh = str2double(rem_thresh);
emg_thresh = get(hs.emgthresh_box,'String');
hs.emg_thresh = str2double(emg_thresh);
movt_thresh = get(hs.movtthresh_box,'String');
hs.movt_thresh = str2double(movt_thresh);

brain_state = hs.sleep_state_data;
if hs.first_run == 1
    
    %first classification point
    if (move_average(1) > hs.movt_thresh*1.2) || (emg_average(1) > hs.emg_thresh*1.2)
         brain_state(1) = 1;
    elseif t_to_d_ratio(1) > hs.rem_thresh*2
        brain_state(1) = 2;
    elseif d_to_b_ratio(1) > hs.nrem_thresh
        brain_state(1) = 3;
    else
        brain_state(1) = 1;
    end
    
    %second classification point
    if (move_average(2) > hs.movt_thresh*1.2) || (emg_average(1) > hs.emg_thresh*1.2)
         brain_state(2) = 1;
    elseif t_to_d_ratio(2) > hs.rem_thresh*2
        brain_state(2) = 2;
    elseif d_to_b_ratio(2) > hs.nrem_thresh
        brain_state(2) = 3;
    else
        brain_state(2) = 1;
    end
    
    start_loop = 3;
else
    start_loop = 1;
end


%classification loop
for min = start_loop:length(t)
    one_back = brain_state(end);
    two_back = brain_state(end-1);
    if two_back == 1
        if one_back == 1
            if (move_average(min) > hs.movt_thresh*0.8) || (emg_average(1) > hs.emg_thresh*0.8)
                brain_state(end+1) = 1;
            elseif t_to_d_ratio(min) > hs.rem_thresh*3
                brain_state(end+1) = 2;
            elseif d_to_b_ratio(min) > hs.nrem_thresh*1.1
                brain_state(end+1) = 3;
            else
                brain_state(end+1) = 1;
            end
        end 
        
        if one_back == 2
            if (move_average(min) > hs.movt_thresh*1.1) || (emg_average(1) > hs.emg_thresh*1.1)
                brain_state(end) = 1;
                brain_state(end+1) = 1;
            elseif t_to_d_ratio(min) > hs.rem_thresh
                brain_state(end+1) = 2;
            elseif d_to_b_ratio(min) > hs.nrem_thresh*1.1
                try
                if (move_average(min-1) > hs.movt_thresh*0.25) || (emg_average(1) > hs.emg_thresh*0.25)
                    brain_state(end) = 1;
                end                    
                brain_state(end+1) = 3;
                end
            else
                brain_state(end) = 1;
                brain_state(end+1) = 1;
            end
        end         
        
        if one_back == 3
            if (move_average(min) > hs.movt_thresh*1.2) || (emg_average(1) > hs.emg_thresh*1.2)
                brain_state(end+1) = 1;
            elseif t_to_d_ratio(min) > hs.rem_thresh*1.1
                brain_state(end+1) = 2;
            elseif d_to_b_ratio(min) > hs.nrem_thresh
                brain_state(end+1) = 3;
            else
                brain_state(end+1) = 1;
            end
        end                 
    end
    
    if two_back == 2
        if one_back == 1
            if (move_average(min) > hs.movt_thresh*1) || (emg_average(1) > hs.emg_thresh*1.2)
                brain_state(end+1) = 1;
            elseif t_to_d_ratio(min) > hs.rem_thresh*1
                try
                if (move_average(min-1) < hs.movt_thresh*0.9) || (emg_average(1) < hs.emg_thresh*1.2)
                    brain_state(end) = 2;
                end
                brain_state(end+1) = 2;
                end
            elseif d_to_b_ratio(min) > hs.nrem_thresh*1            
                brain_state(end+1) = 3;
            else
                brain_state(end+1) = 1;
            end
        end 
        
        if one_back == 2
            if (move_average(min) > hs.movt_thresh*1.2) || (emg_average(1) > hs.emg_thresh*1.2)
                brain_state(end+1) = 1;
            elseif t_to_d_ratio(min) > hs.rem_thresh*0.85
                brain_state(end+1) = 2;
            elseif d_to_b_ratio(min) > hs.nrem_thresh*0.85
                brain_state(end+1) = 3;
            else
                brain_state(end+1) = 1;
            end
        end         
        
        if one_back == 3
            if (move_average(min) > hs.movt_thresh*1.1) || (emg_average(1) > hs.emg_thresh*1.1)
                brain_state(end+1) = 1;
            elseif t_to_d_ratio(min) > hs.rem_thresh*1.1
                brain_state(end+1) = 2;
            elseif d_to_b_ratio(min) > hs.nrem_thresh
                brain_state(end+1) = 3;
            else
                brain_state(end+1) = 1;
            end
        end                 
    end     
        

    if two_back == 3
        if one_back == 1
            if (move_average(min) > hs.movt_thresh*1) || (emg_average(1) > hs.emg_thresh*1)
                brain_state(end+1) = 1;
            elseif t_to_d_ratio(min) > hs.rem_thresh*1
                rem_ratio = t_to_d_ratio(min-1)/hs.rem_thresh;
                nrem_ratio = d_to_b_ratio(min-1)/hs.nrem_thresh;
                
                try
                if (move_average(min-1) < hs.movt_thresh*0.9) || (emg_average(1) < hs.emg_thresh*0.9)
                    if rem_ratio > nrem_ratio
                        brain_state(end) = 2;
                    else
                        brain_state(end) = 3;
                    end
                end
                catch
                    if rem_ratio > nrem_ratio
                        brain_state(end) = 2;
                    else
                        brain_state(end) = 3;
                    end
                end
                brain_state(end+1) = 2;
            elseif d_to_b_ratio(min) > hs.nrem_thresh*1
                brain_state(end+1) = 3;
            else
                brain_state(end+1) = 1;
            end
        end 
        
        if one_back == 2
            if (move_average(min) > hs.movt_thresh*1.2) || (emg_average(1) > hs.emg_thresh*1.2)
                brain_state(end+1) = 1;
            elseif t_to_d_ratio(min) > hs.rem_thresh*0.9
                brain_state(end+1) = 2;
            elseif d_to_b_ratio(min) > hs.nrem_thresh*0.85
                brain_state(end+1) = 3;
            else
                brain_state(end+1) = 1;
            end
        end         
        
        if one_back == 3
            if (move_average(min) > hs.movt_thresh*1.25) || (emg_average(1) > hs.emg_thresh*1.25)
                brain_state(end+1) = 1;
            elseif t_to_d_ratio(min) > hs.rem_thresh*1
                brain_state(end+1) = 2;
            elseif d_to_b_ratio(min) > hs.nrem_thresh*0.9
                brain_state(end+1) = 3;
            else
                brain_state(end+1) = 1;
            end
        end                 
    end
    
end
    
brain_state = brain_state(end-length(t)+1:end);
hs.brain_state = repelem(brain_state,hs.upsamp_by);

end

%%%
end

    