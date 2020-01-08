function [thisData, frequency, thiscell]  = mini_detector(RAW_DATA,inputs)

% Template based miniEPSC/IPSC event detection.

try
%% Setting initial vars
% If there is not this item in the filter list, set it to this

readFilters = inputs; % list of filter parameters given to function

warning('off','signal:findpeaks:largeMinPeakHeight') % suppress warning text for findpeaks

% if = 1, will print when mini events are thrown out based on some
% threshold
verbose_on = 0;

if ~isfield(readFilters,'Smoothingfxn')
    readFilters.Smoothingfxn = 0;
end
if ~isfield(readFilters,'Templatethresh')
    readFilters.Templatethresh = 5.5;
end
if ~isfield(readFilters,'Ampthresh')
    readFilters.Ampthresh = -4;
end
if ~isfield(readFilters,'Risetime')
    readFilters.Risetime = 0.006;
end
if ~isfield(readFilters,'temp_length') % Later called template length
    readFilters.temp_length = 20;
end
if ~isfield(readFilters,'excisecoord')
    readFilters.excisecoord = [0 0];
end
if ~isfield(readFilters,'samplerate')
    readFilters.samplerate = 10000;
end
if ~isfield(readFilters,'pre_min_amp')
    readFilters.pre_min_amp = 5;
end
if ~isfield(readFilters,'matching_skip')
    readFilters.matching_skip = 1;
end
if ~isfield(readFilters,'epsc_yes')
    readFilters.epsc_yes = 1;
end

if readFilters.epsc_yes == 0
    sign_change = -1;
else
    sign_change = 1;
end

thisData = RAW_DATA; % create another variable for inputted data

thisData = sign_change*(thisData); % convert trace data from pA's to A's. Might be scaling issues.


%% Pre-Processing


% high pass filter
filt = designfilt('bandpassiir',...
    'FilterOrder',2, ...
    'HalfPowerFrequency1', 4,...
    'HalfPowerFrequency2',2000, ...
    'SampleRate',readFilters.samplerate);

% apply high pass filter to trace and mild smoothing
if readFilters.Smoothingfxn == 1
    try
        filt_data = filter(filt,detrend(thisData));
    catch
        filt_data = detrend(thisData);
    end
    temp_DATA = sgolayfilt(filt_data,3,15); % create temp smoothed data for template matching
    DATA = smooth(filt_data,3); % mild 3 pt smoothing for analysis
else
    temp_DATA = detrend(thisData);
    DATA = smooth(detrend(thisData),3);
end
    


% Exclusion Criteria
amp_cutoff_low = readFilters.pre_min_amp*10^-12;

refractory_per = 0.0045*readFilters.samplerate; % minumum number of points between events (REFRACTORY PERIOD)
S_E_L = .02*readFilters.samplerate; % Standard event length
template_thresh = readFilters.Templatethresh;
event_start_thresh = 0.014*readFilters.samplerate;% Maximum number of points before fitted event start index
glob.max_rise = 0.01*readFilters.samplerate; % maximum length between start and peak of mini
peak_window = 1:(0.008*readFilters.samplerate); %index window for findpeaks function

%% Template Matching
%build template

sample_rate        = 10000;
tRISE           = 15;
tDECAY           = 25;
rise_time = 15;
pleateau_t = 3;
temp_baseline_length = 60;
templatelength = readFilters.temp_length*sample_rate/1000;
template       = zeros(1,templatelength); 

% This creates the template shape in the list 'template' make rise
for t = 1:rise_time
    template(t) = (exp((t)/tRISE));
end
template = template - template(1);

%add peak plateau
template(rise_time:rise_time+pleateau_t) = template(rise_time);
%make decay
peak = template(rise_time);
for t = (rise_time+pleateau_t+1):length(template)
    template(t) = peak*(exp((-t+rise_time+pleateau_t)/tDECAY));
end

template = -template./(max(template));
template = [zeros(1,temp_baseline_length*sample_rate/10000) template]; %this forces a stable baseline prior to each mini. Adds 10 zeros.

%%%
SCALE=[];
OFFSET=[];
first_to_pk_amp=[];
SSE=zeros([1 length(temp_DATA) - length(template)]);
fitted_template_mat=zeros([length(temp_DATA) - length(template) length(template)]);
STANDARDERROR=[];
matching_skip = 1;

max_offset_dev = 3.5e-12;
indx_ahead_to_scale = 28e-4*readFilters.samplerate;
indx_behind_to_scale = 0;

SSE_tot = [];
SSE_pk = [];
SSE_tail = [];
SSE_pk2 = [];
temp_event_corrs = [];
pk_amps = [];

indx= 1;
start_indx = 1;
for ii = start_indx:matching_skip:(length(temp_DATA) - length(template))
    
%     if mod(ii,length(temp_DATA)/10) == 0
%         disp(['Finished ',num2str(ii/(length(temp_DATA)))])
%     end
    
    dat = temp_DATA(ii:ii+length(template)-1); % sliding window the size of template

    OFFSET = median(dat(1:temp_baseline_length-25));
    if OFFSET > max_offset_dev
        OFFSET = max_offset_dev;
    elseif OFFSET < -max_offset_dev
        OFFSET = -max_offset_dev;
    end
    
    %first_to_pk_amp(indx) = abs(fitted_template(1,1)-fitted_template(1,pk_ind));
    first_to_pk_amp = abs(min(dat(temp_baseline_length-indx_behind_to_scale:temp_baseline_length+indx_ahead_to_scale))-OFFSET);
    fitted_template = template*(first_to_pk_amp)*1.1 + OFFSET; % apply scaling to template and offset
    if first_to_pk_amp > amp_cutoff_low
%         SSE_tot(indx) = sum((offset_line - dat.').^2)*2;
%         SSE_pk(indx) = sum((fitted_template(temp_baseline_length-10:temp_baseline_length+70) - dat(10:90).').^2)*5;
        SSE_pk2(indx) = 1.5*sum((fitted_template(temp_baseline_length:temp_baseline_length+30) - dat(temp_baseline_length:temp_baseline_length+30).').^2);
        
        SSE_tail(indx) = sum((fitted_template(temp_baseline_length+30:temp_baseline_length+60) - dat(temp_baseline_length+30:temp_baseline_length+60).').^2)*0.75;
        
        pk_amps(indx) = first_to_pk_amp;
%         corr_coef = corrcoef(fitted_template(temp_baseline_length:temp_baseline_length+50).', dat(temp_baseline_length:temp_baseline_length+50));
%         temp_event_corrs(indx) = corr_coef(2);
%         event_area(indx) = sum(-OFFSET + dat(50:100).');
    else
%         event_area(indx) = 0;
        SSE_pk2(indx) = 0;
        SSE_tail(indx) = 0;
        pk_amps(indx) = 0;
    end
    
%     startindx = 7.220*readFilters.samplerate;
%     if indx > startindx        
%         fig10 = figure(10);
%         subplot(2,1,1)
%         plot(dat); hold on;
%         plot(fitted_template);
%         text(20,-first_to_pk_amp,num2str((SSE_pk2(indx)+SSE_tail(indx))),'FontSize',18)
%         hold off
%         subplot(2,1,2)
%         plot(indx-startindx,(SSE_pk2(indx)+SSE_tail(indx)),'ro')
%         xlim([0 1.5*(indx-startindx)])
%         hold on
%         pause(0.1)
% %       
%         if ~ishandle(fig10); keyboard; end
%     end
    indx = indx+1;
end

SSE_pk_and_tail = SSE_pk2 + SSE_tail;
x_axis = 1/(readFilters.samplerate):(1/readFilters.samplerate):length(SSE_pk_and_tail)/readFilters.samplerate;

% figure;plot(x_axis,SSE_pk_and_tail)
% 
% sm_SSE_pk = smooth(SSE_pk2, 15);
% norm_SSE_pk = sm_SSE_pk(25:end).'./SSE_pk2(1:end-24);
% norm_SSE_pk(norm_SSE_pk==inf) = 0;
% norm_SSE_pk(isnan(norm_SSE_pk)) = 0;

% make vector with 1's for nonzeros
SSE_event_hits = [0 diff((SSE_pk_and_tail > 0))];
SSE_event_hits(SSE_event_hits<0) = 0;
event_inds = find(SSE_event_hits);

% make vect with num of consecutive nonzero pts
i = find(diff((SSE_pk_and_tail > 0))) ;
n = [i numel(SSE_event_hits)] - [0 i];
c = arrayfun(@(X) X-1:-1:0, n , 'un',0);
num_consec = cat(2,c{:});

min_wind_st = round(7.5e-4*readFilters.samplerate); % first pt for finding SSE minimum
SSE_diff = [];
pk_inds = [];
for event_num = 1:length(event_inds)
        
    event_start = event_inds(event_num);
    event_sse = SSE_pk_and_tail(event_start:(event_start+num_consec(event_start)));
    event_pks = pk_amps(event_start:(event_start+num_consec(event_start)));
    if length(event_sse) > min_wind_st*1.5
        try
            min_window = min_wind_st:(indx_ahead_to_scale-rise_time)*2;
            [min_sse, min_ind] = min(event_sse(min_window));
            max_pk = max(event_pks(min_window));
        catch % if the event SSE interval is too short
            min_window = min_wind_st:length(event_sse);
            [min_sse, min_ind] = min(event_sse(min_window));
            max_pk = max(event_pks(min_window));
        end
        pk_inds(event_num) = min_ind+event_start;
        max_sse = max(event_sse);
        SSE_diff(event_num) = max_sse/min_sse; % ratio of max to min SSE
        pk_vals(event_num) = max_pk;
    else
        pk_inds(event_num) = NaN;
        max_sse = NaN;
        SSE_diff(event_num) = NaN;
        pk_vals(event_num) = NaN;
    end
end

%This part creates an adjusted threshold that is lower for larger events
%using a decay function. larger decay const means faster drop off --> more
%forgiving of SSE in detecting larger events.
timeindx = [];
if isempty(event_inds) == 0
    adj_decay_const = readFilters.AmpBiasCoeff;
%     adj_decay_const = 0.15;
    adj_const = (amp_cutoff_low*10^12+1)*exp(-adj_decay_const.*pk_vals*10^12);
    ymax = (amp_cutoff_low*10^12+1)*exp(-adj_decay_const.*(amp_cutoff_low*10^12+1));
    adj_const = adj_const./(ymax);
    size_adj_thresh = (template_thresh)*adj_const;

    timeindx = pk_inds(SSE_diff > size_adj_thresh);
    
    inds_base_before_pk = 40e-4*readFilters.samplerate;
    timeindx = timeindx - inds_base_before_pk + temp_baseline_length + rise_time;
end

% figure;hist(adj_const,100)
% figure;hist(size_adj_thresh,100)
% figure;plot(size_adj_thresh, SSE_diff,'o'); refline(1,0);


% timeindx = event_inds(SSE_diff > template_thresh)+temp_baseline_length+timeindx_shift;
disp(['num of hits = ', num2str(length(timeindx))])

% figure(17);plot(SSE_pk_and_tail);title('SSE_both');
% figure(19);plot(SSE_diff);title('SSE Diff')

%% Peak/Start Detection
% Loop through detected events and estimate where time of event start,
% peak is, and decay ends. Also extracts peak amplitude and rejects events
% based on rejection parameters

% Removes portion of trace in excisecoord filter
% Find a way to put this before the detection of mini's?
for coord_pair = 1:size(readFilters.excisecoord,1)
    timeindx(find(timeindx > readFilters.excisecoord(coord_pair,1) & ...
    timeindx < readFilters.excisecoord(coord_pair,2))) = [];
end
    
smoothed_events = [];

%%%
samp_rate = readFilters.samplerate;
rise_cutoff = readFilters.Risetime*samp_rate;
length_cutoff = 0.0025*samp_rate;

expand_event_by = 0.0000*samp_rate;
rise_time = nan(1,length(timeindx));
decay_ind = nan(1,length(timeindx));
sm_pk_val = nan(1,length(timeindx));
sm_pk_ind = nan(1,length(timeindx));
event_start_ind = nan(1,length(timeindx));
event_start_val = nan(1,length(timeindx));
event = nan(S_E_L+expand_event_by+1, length(timeindx));
excl_mini_timeindx = zeros([1 length(timeindx)]);            

for mini_num = 1:numel(timeindx)
    try
    
    max_rise = glob.max_rise;
    
    % boundary fix
    if timeindx(mini_num)+S_E_L < length(DATA)
        event(:,mini_num) = DATA(timeindx(mini_num)-expand_event_by:timeindx(mini_num)+S_E_L); % create 'event' that contains DATA 
    else
        event(:,mini_num) = DATA(timeindx(mini_num)-expand_event_by:end);
    end
    
    %%%% for testing
%     if mini_num == 74
%         figure;
%         plot(event(:,mini_num))
%         disp('stopped for adjustments')
% 
%     end
    %%%%
    
    % process mini event
    mini_threshold = (readFilters.pre_min_amp - 1)*10^-12; % threshold for detecting mini peak used here
    smooth_event = sgolayfilt(event(:,mini_num), 2e-4*samp_rate, 9e-4*samp_rate); % smooths event with polynomial 2, window size 13
    smoothed_events(:,mini_num) = smooth_event;
    offset_smooth_event = smooth_event - mean(smooth_event(1:15e-4*samp_rate));
    
    % find initial event peak (which is negative)
    int_st = 0.0005*samp_rate;
    int_end = 0.0150*samp_rate;
    try
        [min_pks, min_locs] = findpeaks(-offset_smooth_event(peak_window), 'MinPeakHeight', mini_threshold);
    catch
        [min_pks, min_locs] = min(offset_smooth_event(int_st:int_end));
    end
    
    if isempty(min_locs) == 1
        [min_pks, min_locs] = min(offset_smooth_event(1:int_end));
    end

    [largest_pk, ~] = max(min_pks);
    
    mini_pk_loc = min_locs(find(abs(min_pks) > 0.75*abs(largest_pk), 1)); % pick large enough first peak

    % make sure you don't index out of bounds
    if max_rise >= mini_pk_loc
        max_rise = mini_pk_loc - 1;
    end
    
    % find potential pre-rise pk before mini pk
    try
        [max_pks, max_locs] = findpeaks(event(mini_pk_loc-max_rise:mini_pk_loc, mini_num));
    catch
        %disp(['pre-rise mini start finding problem', ' mini ', num2str(mini_num)])
        [max_pks, max_locs] = max(smooth_event(mini_pk_loc-max_rise:mini_pk_loc));
    end
    
    if isempty(max_locs) == 1
        [max_pks, max_locs] = max(smooth_event(mini_pk_loc-max_rise:mini_pk_loc));
    end
        
    % using median of pre-peaks find a good pre-peak candidate
    pks_above_med = (max_pks >= median(smooth_event(mini_pk_loc-max_rise:mini_pk_loc - round(max_rise/2))));
    approx_mini_sizes = abs(event(mini_pk_loc,mini_num)) +  max_pks(pks_above_med);
    locs_above_med = max_locs(pks_above_med);
    pre_pk_loc = locs_above_med(find(approx_mini_sizes > mini_threshold, 1, 'last'));
    if isempty(pre_pk_loc) == 1
        try
            pre_pk_loc = max_locs(round(length(max_locs)/2));
        catch
            pre_pk_loc = round(mini_pk_loc/2);
        end
    end
    
    pre_pk_val = nanmean([nanmean(event(1:pre_pk_loc,mini_num)),...
        nanmedian(event(1:pre_pk_loc,mini_num))]);
    
    % specify ev start
    event_start_ind(mini_num) = (mini_pk_loc - max_rise + pre_pk_loc - 1); %pick last one and find time to beginning
    event_start_val(mini_num) = pre_pk_val;

    if mini_pk_loc-20e-4*samp_rate < 1 % prevent mini pks from being too early
        mini_pk_loc = 21e-4*samp_rate;
    end

    % another pass at mini peak
    [sm_pk_val(mini_num),sm_pk_ind(mini_num)] = min(event((mini_pk_loc-20e-4*samp_rate:mini_pk_loc+30e-4*samp_rate),mini_num)); 
    sm_pk_ind(mini_num) = mini_pk_loc - 20e-4*samp_rate + sm_pk_ind(mini_num) - 1; % correct placement
    sm_pk_val(mini_num) = abs(sm_pk_val(mini_num)) + event_start_val(mini_num);

    % eliminate mini's too close to one another
    if mini_num>1
        if timeindx(mini_num) - timeindx(mini_num-1) <= refractory_per   %%%%%%% minimum time between events (refractory period)
            excl_mini_timeindx(mini_num) = 1;            
            %timeindx(mini_num) = NaN;
            if verbose_on == 1
            disp(['too close ', num2str(mini_num),' indx: ',num2str(timeindx(mini_num))])
            end

        end
    end

    if isempty(event_start_ind(mini_num)) == 1
        event_start_ind(mini_num) = 1;
    end

    if event_start_ind(mini_num) >= event_start_thresh % event start index must be below or equal to threshold
        excl_mini_timeindx(mini_num) = 1;            
        %timeindx(mini_num) = NaN;
        if verbose_on == 1
        disp(['event start >= threshold ', num2str(mini_num),' indx: ',num2str(timeindx(mini_num))])
        end
    end

    if ~isnan(sm_pk_ind(mini_num))
        % find where mini has decreased by 90% from peak
        decay_loc = find((smooth_event(sm_pk_ind(mini_num):end) - event_start_val(mini_num))...
            > -0.2*sm_pk_val(mini_num), 1, 'first');

        if isempty(decay_loc) == 0
            decay_ind(mini_num) = decay_loc + mini_pk_loc - 1;
        else
            %Roughest error catch
            %fprintf(['decay error -> could not find decay', 'mini ', num2str(mini_num)])
            decay_ind(mini_num) = 0.008*samp_rate + mini_pk_loc;
        end

        %Find peak again after having found decay pt
        first_sm_pk_ind = sm_pk_ind(mini_num);

        [sm_pk_val(mini_num),sm_pk_ind(mini_num)] = min(smooth_event(first_sm_pk_ind-20e-4*samp_rate:...
            decay_ind(mini_num)));

        sm_pk_ind(mini_num) = first_sm_pk_ind - 20e-4*samp_rate + sm_pk_ind(mini_num) - 1; % correct placement
        sm_pk_val(mini_num) = abs(sm_pk_val(mini_num)) + event_start_val(mini_num);

        if sm_pk_val(mini_num) < abs(readFilters.Ampthresh*10^-12)
            if verbose_on == 1
                disp(['Mini too small ', num2str(mini_num),' indx: ',num2str(timeindx(mini_num))])
            end
            excl_mini_timeindx(mini_num) = 1;            
        end

        rise_time(mini_num) = sm_pk_ind(mini_num) - event_start_ind(mini_num);
        if rise_time(mini_num) > rise_cutoff
            excl_mini_timeindx(mini_num) = 1;
            if verbose_on == 1
                disp(['large rise ', num2str(mini_num),' indx: ',num2str(timeindx(mini_num))])
            end
        end

        decay_loc = find((smooth(event(sm_pk_ind(mini_num):end,mini_num),10) - event_start_val(mini_num))...
            > -0.15*sm_pk_val(mini_num), 1, 'first');

        if decay_loc > 200e-4*samp_rate
            decay_loc = find((smooth(event(sm_pk_ind(mini_num):end,mini_num),10) - event_start_val(mini_num))...
            > -0.3*sm_pk_val(mini_num), 1, 'first');
        end

        if isempty(decay_loc) == 0
            decay_ind(mini_num) = decay_loc + sm_pk_ind(mini_num);
        else
            % Rough catch if decay cannot be found
            decay_ind(mini_num) = 80e-4*samp_rate + mini_pk_loc;
        end

        event_duration = decay_ind(mini_num) - event_start_ind(mini_num);
        if event_duration < length_cutoff
            if verbose_on == 1
                disp(['too short ', num2str(mini_num),' indx: ',num2str(timeindx(mini_num))])
            end
            excl_mini_timeindx(mini_num) = 1;
        end
    end

%     figure;plot(event(:,mini_num)); NumTicks = 15;
%     L = get(gca,'YLim');
%     set(gca,'YTick',linspace(L(1),L(2),NumTicks))
%     
%     disp(['amp ', num2str(sm_pk_val(mini_num))])
%     disp(['peak ind ', num2str(sm_pk_ind(mini_num))])
%     disp(['event st ', num2str(event_start_ind(mini_num))])
%     disp(['event st value ', num2str(event_start_val(mini_num))])
%     disp(['decay ind ', num2str(decay_ind(mini_num))])

    if mini_num > 1
        prev_mini_decay_time = timeindx(mini_num-1) + decay_ind(mini_num-1);
        this_mini_event_start = timeindx(mini_num) + event_start_ind(mini_num);
        if prev_mini_decay_time >= this_mini_event_start
            decay_ind(mini_num-1) = this_mini_event_start - timeindx(mini_num-1) - 1;
        end
    end
    
    catch
        if verbose_on == 1
            disp(['mini characterization error ',num2str(mini_num),' indx: ',num2str(timeindx(mini_num))])
        end
        excl_mini_timeindx(mini_num) = 1;            
    end
    
%     subplot(1,2,1)
%     plot(event(:,mini_num))
%     subplot(1,2,2)
%     plot(smooth_event)
%     waitforbuttonpress
end

%% The Rest

excl_mini_timeindx = logical(excl_mini_timeindx);
deletethis = excl_mini_timeindx;

S_E_L = S_E_L + expand_event_by;
timeindx = timeindx - expand_event_by;

excl_mini_timeindx = timeindx(excl_mini_timeindx);

RAW_AMP = sm_pk_val;
RAW_AMP = 10^12*RAW_AMP;

timeindx(deletethis) = [];
RAW_AMP(deletethis) = [];
rise_time(deletethis) = [];
decay_ind(deletethis) = [];
sm_pk_val(deletethis) = [];
sm_pk_ind(deletethis) = [];
event_start_ind(deletethis) = [];
event_start_val(deletethis) = [];
event(:,deletethis) = [];

disp(['Events Detected ', num2str(length(timeindx))])

wavgs = NaN(S_E_L,numel(timeindx));

% Shift the events so that starts align
shift_axis = 40e-4*samp_rate;
shifted_events_trace = [];
for mini = 1:length(timeindx)

    start_ind = event_start_ind(mini);
    pk_ind = sm_pk_ind(mini);
    pk_val = RAW_AMP(mini)*10^-12;
    mini_raw = event(:,mini);  

    if pk_val > 5e-12
        smooth_mini = smooth(mini_raw(~isnan(mini_raw)),20e-4*samp_rate,'sgolay',3)';

        baseline_cur = nanmedian(smooth_mini(1:start_ind));
        smooth_mini = smooth_mini - baseline_cur;

        end_baseline = nanmean(mini_raw(50e-4*samp_rate:end));

        sm_mini_raw = smooth(mini_raw(~isnan(mini_raw)),13e-4*samp_rate,'sgolay',3)';
        sm_mini_raw = sm_mini_raw - baseline_cur;

        near_pk_ind = find(sm_mini_raw < (sm_mini_raw(start_ind)-pk_val)*.80,1);
        shift_thresh = (sm_mini_raw(start_ind)-pk_val)/2;
        [~, half_way_slope] = min(abs(smooth_mini(start_ind:near_pk_ind)-shift_thresh));
        half_way_slope = half_way_slope + start_ind -1;
        shift_ind_2 = shift_axis - half_way_slope -1;

        baseline_mini_raw = (mini_raw(~isnan(mini_raw)) - baseline_cur)';
        if shift_ind_2 > 0 
            shifted_event = [NaN(1,shift_ind_2) baseline_mini_raw];
        elseif shift_ind_2 < 0 
            shifted_event = baseline_mini_raw(-shift_ind_2:end);
        else
            shifted_event = baseline_mini_raw(shift_ind_2+1:end);
        end
%         if end_baseline>-6e-12
        shifted_events_trace = padmat(shifted_events_trace, shifted_event, 1);
%         end

    end
end

shifted_events = shifted_events_trace;

total_time = (numel(DATA)-numel(DATA(isnan(DATA))))*(1/readFilters.samplerate);
frequency = numel(~isnan(timeindx))/total_time;
% 

thiscell = [];
thiscell.timeindx = timeindx;
thiscell.sm_pk_ind = sm_pk_ind;
thiscell.event_start_ind = event_start_ind;
thiscell.event_start_val = event_start_val;
thiscell.decay_ind = decay_ind;
thiscell.wavgs = wavgs;
thiscell.template = template;
thiscell.DATA = DATA;
thiscell.n = readFilters.temp_length;
thiscell.detectcrit = SSE_pk_and_tail;
thiscell.rawdata = thisData;
thiscell.amp = RAW_AMP;
thiscell.risetime = rise_time;
thiscell.freq = frequency;
thiscell.events = event;
thiscell.shifted_events = shifted_events;
thiscell.excl_minis = excl_mini_timeindx;

catch
    disp(['error in mini analysis fxn']);    
end

end




