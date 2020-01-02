function [Cp_est, Rin] = get_PassProp_IC(data, step_start, step_length, istep, samprate)

% data - can either be a vector containing the volt trace with sealtest or
% 'choose', which will open up a browser to pick the file
% step start - this is where the sealtest starts, should be in seconds
% vstep - the voltage of the step
% the sampling rate

% %Example parameters
% data = raw_trace;
% step_start = 0.1;
% step_length = 0.4;
% istep = 100e-12;
% samprate = 10000;


% 
if strcmp(data,'choose')
    [hs.zFile, hs.zDir] = uigetfile('*.ibw','Select folder with mini data');
    data = IBWread([hs.zDir, hs.zFile]);
    data = data.y;
end
% 
figures_on = 1;%% 0 prevents plotting

% figure;plot(data)

iclamp_sealtest_current = istep;

seal_volt_ind = round((step_start+0.05)*samprate:(step_start+step_length-0.05)*samprate);
steady_state_ind = round((step_start+step_length+0.05)*samprate:(step_start+step_length+0.45)*samprate);
sealtest_volt = mean(data(seal_volt_ind));
base_volt = mean(data(steady_state_ind));
Rin = abs((sealtest_volt-base_volt)) ./ (iclamp_sealtest_current);


peak = step_start*samprate+1;
fit_end = (step_start + 0.05)*samprate;
fit_data = (data(round(peak):round(fit_end)));
fit_data = fit_data - mean(fit_data(round(0.04*samprate):end));


tau_est = (find(fit_data<(fit_data(1)*0.37), 1) - 1)/(samprate);

time = (0:(1/samprate):(length(fit_data)/samprate)-(1/samprate)).';
% s = fitoptions('Method', 'NonlinearLeastSquares', ...
%     'StartPoint', [transient_vals(1), tau_est*0.5, transient_vals(10), tau_est*1.5],...
%      'Lower', [mean(transient_vals(1:2)*.9), tau_est*0.001, mean(transient_vals(1:5))*.1, tau_est*0.01],...
%      'Upper', [transient_vals(1)*1.1, tau_est*5, transient_vals(1)*1.25, tau_est*15]);
% f = fittype('a*exp(-x/b) + c*exp(-x/d)','options',s);

s = fitoptions('Method', 'NonlinearLeastSquares', ...
    'StartPoint', [fit_data(1), tau_est*0.5],...
     'Lower', [mean(fit_data(1:2)*.85), tau_est*0.001],...
     'Upper', [fit_data(1)*1.15, tau_est*5]);
f = fittype('a*exp(-x/b)','options',s);

% watchfithappen(time, fit_data, f, 10)

[exp_fit,~] = fit(time,fit_data,f);
cval = coeffvalues(exp_fit);
tau = cval(2);

Cp_est = tau/Rin;

if figures_on == 1
%     figure(4);
%     plot(fit_data)
%     figure(5);
%     plot(exp_fit,time,fit_data)
%     disp(istep)
    disp(['Cp (pf): ',num2str(Cp_est*1e12),'  Rin (MOhms): ',num2str(Rin*1e-6)])
    %text(10,50,num2str(Rs*1e-6))
    %axis([0 50 min(transient_seal_vals) 100])
end


% disp(Rs_nate/10^6);
% disp(Rt/10^6);
% disp(Cm*10^12);
%disp(['Vm is ', num2str(Vm*10^3)]);