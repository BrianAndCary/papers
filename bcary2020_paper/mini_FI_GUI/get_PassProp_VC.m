function [Ra, Cp, Rin, Vr] = get_PassProp_VC(data, step_start, step_length, vstep, samprate)

% data - can either be a vector containing the volt trace with sealtest or
% 'choose', which will open up a browser to pick the file
% step start - this is where the sealtest starts, should be in seconds
% vstep - the voltage of the step
% the sampling rate




if strcmp(data,'choose')
    [hs.zFile, hs.zDir] = uigetfile('*.ibw','Select folder with mini data');
    data = IBWread([hs.zDir, hs.zFile]);
    data = data.y;
end

figures_on = 1;%% 0 prevents plotting


V_test = -abs(vstep); %test depolarization
V_hold = -7e-2; %holding potential


t=step_start; %starting point for fit
tracevals = data(t*samprate:(t*samprate+9000));

ss_Start = 400; %timepoints after the peak, used to define the steady state region
ss_End = 900; %timepoints after the peak, used to define the steady state region

pulse_end = step_length*samprate;
peak = max(tracevals(pulse_end:pulse_end+ss_End));
startfit = find(tracevals==peak,1,'first');% offset to avoid the effects of pipette capacitance 

I_test = nanmedian(tracevals(pulse_end - 600:pulse_end - 100)); %calculate the average baseline current
I_test_std = nanstd(tracevals(pulse_end - 600:pulse_end - 100));
I_baseline = nanmedian(tracevals(pulse_end + ss_Start:pulse_end + ss_End)); %calculate the average steady state current during the test
dI = I_test - I_baseline; 
% I_peak = tracevals(startfit(1)-10) - I_baseline;
% Rs = V_test/I_peak %calculate series resistance based on peak amplitude of depolarization, which is unreliable
Rin = abs(V_test/dI); %total resistance 
% Rc = V_test/dI - Rs; %input resistance

Vr = -(I_baseline * Rin) + V_hold; 


whole_seal_fit_x = (startfit:startfit + ss_End)';
ss_seal_fit_x = (startfit + ss_Start:startfit + ss_End)';  %%%%%%%steady state for BL subtract
baseline = nanmedian(tracevals(ss_seal_fit_x));


transient_vals = tracevals(whole_seal_fit_x) - baseline;
transient_vals = transient_vals*10^12;


[~,pk_ind] = max(tracevals);
charge_st = pk_ind - 1;
transient_seal_x = (charge_st:charge_st + ss_Start)';
transient_seal_vals = tracevals(transient_seal_x) - baseline;
Q1 = sum((abs(transient_seal_vals(1:200).*(1/samprate))));    %%%%%% divide by sam. rate

% figure;plot(abs(transient_seal_vals(1:200)));

tau_est = (find(transient_vals<(transient_vals(1)*0.37), 1) - 1)/(samprate);

transient_st = find(transient_vals<(transient_vals(1)*1.1),1);
transient_end = find(transient_vals<(transient_vals(1)*0.08),1);
transient_st = 1;
transient_end = 90;

transient_vals = transient_vals(transient_st:transient_end);

time = (0:(1/samprate):(length(transient_vals)/samprate)-(1/samprate)).';

transient_vals = transient_vals - mean(transient_vals(end-5:end));
s = fitoptions('Method', 'NonlinearLeastSquares', ...
    'StartPoint', [transient_vals(1), tau_est*0.5, transient_vals(10), tau_est*1.5],...
     'Lower', [mean(transient_vals(1:2)*.9), tau_est*0.001, mean(transient_vals(1:5))*.1, tau_est*0.01],...
     'Upper', [transient_vals(1)*1.1, tau_est*5, transient_vals(1)*1.25, tau_est*15]);
f = fittype('a*exp(-x/b) + c*exp(-x/d)','options',s);
% watchfithappen(time, transient_vals, f, 10)

[exp_fit,~] = fit(time,transient_vals,f);
cval = coeffvalues(exp_fit);

fast_exp = cval(1)*exp(-time./cval(2));
slow_exp = cval(3)*exp(-time./cval(4));
% 
% figure; plot(time, fast_exp); hold on
% plot(time, slow_exp);

Q1_fromfit = (sum(fast_exp.*(1/samprate)) + sum(slow_exp(1:10).*(1/samprate)))*10^-12;
Q1_fromfit = (sum(fast_exp.*(1/samprate)))*10^-12;

Q1_minus_slow = Q1 - sum(slow_exp.*(1/samprate))*10^-12;

fast_Tau = cval(2);
slow_Tau = cval(4);

Q2 = dI * fast_Tau;
Qt = Q1 + Q2;
Cp = abs(Qt/V_test);
Ra_from_Q = fast_Tau/Cp; %predicts model cell best but tends to underest.

Rs_corrected = fast_Tau / ((Q2+Q1_minus_slow)/V_test); %almost always overest.

num_pts_back_for_extrap = (12-find(tracevals(startfit-10:startfit)>I_test+5*I_test_std,1));
% disp(num_pts_back_for_extrap)
[vals] = find(tracevals==peak);
if length(vals) < 2
    if num_pts_back_for_extrap > 2
        num_pts_back_for_extrap = 2;
    end
end

baseline_correction = abs(abs(nanmedian(tracevals((startfit+transient_end-10):(startfit+transient_end))))...
    -abs(baseline));

Ra = abs(V_test)/((exp_fit((-(num_pts_back_for_extrap/samprate)))*10^-12) + baseline_correction); %usually closest off by around 2-3 Mohms for model

Cp_est = slow_Tau/Rin;




if figures_on == 1
%     figure(4);
%     plot(tracevals)
%     figure(6);
%     plot(exp_fit,time,transient_vals)
    disp(['Ra (MOhms): ',num2str(Ra*1e-6),'  Cp (pf): ',num2str(Cp*1e12),'  Rin (MOhms): ',num2str(Rin*1e-6),'  Vr (mV): ',num2str(Vr*1e3)])
%     text(10,50,num2str(Rs*1e-6))
%     axis([0 50 min(transient_seal_vals) 100])
end
