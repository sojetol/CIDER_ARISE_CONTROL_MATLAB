%% Initialize
add_all_paths_Sandia
load E3SM_params
pattern_T_base = pattern_T_base/288.5983*288;
total_years = 80;
start_year = 2020;
%% CO2
load CO2_concentrations.mat
CO2_ppm_growth = 3;
CO2_SSP245 = 400+(0:total_years+21)'*CO2_ppm_growth;
CO2_ref = CO2_SSP245(4);
CO2_forcing_SSP245 = 5.35*log((CO2_SSP245(5:end)+0)/CO2_ref);
CO2_forcing_SSP245_month = repeatElements(CO2_forcing_SSP245,12);
CO2_SSP245_smooth = linspace(CO2_SSP245(5)-CO2_ppm_growth/2,CO2_SSP245(end)+CO2_ppm_growth/2,length(CO2_forcing_SSP245_month))';
CO2_forcing_SSP245_month_smooth =  5.35*log((CO2_SSP245_smooth)/CO2_ref);

CO2_forcing_SSP370 =  5.35*log((CO2_SSP370(5:end)+0)/CO2_ref);
CO2_forcing_SSP370_month = repeatElements(CO2_forcing_SSP370,12);

% CO2_forcing_SSP245_month = CO2_forcing_SSP245_month_smooth;

%% Figure 1: Controlled SAI
refvals = [288.5, 0.8767, -5.89]; % Where controller is bringing climate to
%           T0       T1      T2

time = annualToMonthly(start_year:start_year+total_years-1);

tseries = 0:1/12:total_years-1/12;
offset = 150;
log_array = [];
log_array_370 = [];
q_array = zeros(4,total_years);
q_array_370 = zeros(4,total_years);
my_title = "Emulating with T0-T1-T2 (Different) PreInd Variability, Plotted annual";
years_range = 1+offset:12*total_years+offset;
load("Tools\PI_T0T1T2.mat")
preset_T_noise = 1*[1*(PI_T0(years_range)-mean(PI_T0(years_range))),...
                    1*(PI_T1(years_range)-mean(PI_T1(years_range))),...
                    1*(PI_T2(years_range)-mean(PI_T2(years_range)))];
% Just to initialize controls
init_year = 2035-start_year+1;
time_length = 12*init_year;
inj_array = [zeros(time_length,1) repeatElements(q_array(4,1:init_year)'/12,12) repeatElements(q_array(3,1:init_year)'/12,12)  zeros(time_length,1) repeatElements(q_array(2,1:init_year)'/12,12) repeatElements(q_array(1,1:init_year)'/12,12)  zeros(time_length,1)];
T_arise_ann_control = CIDER_pattern_from_all_injections_and_CO2([ inj_array CO2_forcing_SSP245_month(1:time_length)],param_AOD_all,param_T_all,pattern_T_all)+pattern_T_base;
T_arise_ann_control_370 = CIDER_pattern_from_all_injections_and_CO2([ inj_array CO2_forcing_SSP370_month(1:time_length)],param_AOD_all,param_T_all,pattern_T_all)+pattern_T_base;
[q_array(:,init_year),log_array] = ARISE_controller_translation(T_arise_ann_control,log_array,refvals,preset_T_noise(1:12,:));
[q_array_370(:,init_year),log_array_370] = ARISE_controller_translation(T_arise_ann_control_370,log_array_370,refvals,preset_T_noise(1:12,:));
% Simulate (annual controls)
tic
for i_year = init_year:total_years
time_length = 12*i_year;
inj_array = [zeros(time_length,1) repeatElements(q_array(4,1:i_year)'/12,12) repeatElements(q_array(3,1:i_year)'/12,12)  zeros(time_length,1) repeatElements(q_array(2,1:i_year)'/12,12) repeatElements(q_array(1,1:i_year)'/12,12)  zeros(time_length,1)];
inj_array_370 = [zeros(time_length,1) repeatElements(q_array_370(4,1:i_year)'/12,12) repeatElements(q_array_370(3,1:i_year)'/12,12)  zeros(time_length,1) repeatElements(q_array_370(2,1:i_year)'/12,12) repeatElements(q_array_370(1,1:i_year)'/12,12)  zeros(time_length,1)];

T_arise_ann_control = CIDER_pattern_from_all_injections_and_CO2([ inj_array CO2_forcing_SSP245_month(1:time_length)],param_AOD_all,param_T_all,pattern_T_all)+pattern_T_base;
T_arise_ann_control_370 = CIDER_pattern_from_all_injections_and_CO2([ inj_array_370 CO2_forcing_SSP370_month(1:time_length)],param_AOD_all,param_T_all,pattern_T_all)+pattern_T_base;

[q_array(:,i_year+1),log_array] = ARISE_controller_translation(T_arise_ann_control,log_array,refvals,preset_T_noise(1:time_length,:));
[q_array_370(:,i_year+1),log_array_370] = ARISE_controller_translation(T_arise_ann_control_370,log_array_370,refvals,preset_T_noise(1:time_length,:));

end
T_no_inj = CIDER_pattern_from_all_injections_and_CO2([inj_array*0 CO2_forcing_SSP245_month(1:time_length)],param_AOD_all,param_T_all,pattern_T_all)+pattern_T_base;
T_no_inj_worse = CIDER_pattern_from_all_injections_and_CO2([inj_array*0 CO2_forcing_SSP370_month(1:time_length)],param_AOD_all,param_T_all,pattern_T_all)+pattern_T_base;
toc
termination_year = 2060;
terminate = [ones((termination_year-start_year)*12,1);zeros(time_length-(termination_year-start_year)*12,1)];
T_term = CIDER_pattern_from_all_injections_and_CO2([inj_array.*terminate CO2_forcing_SSP245_month(1:time_length)],param_AOD_all,param_T_all,pattern_T_all)+pattern_T_base;
T_term_370 = CIDER_pattern_from_all_injections_and_CO2([inj_array_370.*terminate CO2_forcing_SSP370_month(1:time_length)],param_AOD_all,param_T_all,pattern_T_all)+pattern_T_base;
%%
CIDER_params = [];
CIDER_params.param_AOD_all = param_AOD_all;
CIDER_params.param_T_all = param_T_all;
CIDER_params.pattern_T_all = pattern_T_all;
CIDER_params.pattern_T_base = pattern_T_base;

controller_params = [];
controller_params.refvals =  [288.5, 0.8767, -5.89];
controller_params.kivals = [0.0183, 0.0753, 0.3120];
controller_params.kpvals = [0.0183, 0.0753, 0.3120];
controller_params.months_to_average = 12; % This is your rolling filter
controller_params.sampling_period = 12;
controller_params.sens = 4.1;
controller_params.ff_rate = 0.0273;
controller_params.x_ramp = 5;

[a,b,c]=CIDER_run_ARISE_controlled(CO2_forcing_SSP245_month_smooth,35*12,15*12+1,controller_params,CIDER_params,preset_T_noise);
% figure
% plot(globalMean(a))
%% Plot Figure 1
PIT = 287;
figure
hold on
box on
grid on
set(gca,'DefaultLineLineWidth',6)
set(gca,'FontSize',14)
plot(ann(time),ann(globalMean(T_no_inj_worse)+preset_T_noise(1:time_length,1))-PIT,"Color","k")

% plot(ann(time),ann(globalMean(T_no_inj)+preset_T_noise(1:time_length,1))-PIT,"Color","r")

% plot(ann(time),ann(globalMean(T_arise_ann_control)+preset_T_noise(1:time_length,1))-PIT,"Color","b")
% plot(ann(time),ann(globalMean(T_term_370)+preset_T_noise(1:time_length,1))-PIT,"Color",[168, 94, 50]/255)
line([2035 2100], [1.5 1.5],"LineStyle","--","Color",[0.5210    0.0860    0.8190])

% plot(ann(time),ann(globalMean(T_term)+preset_T_noise(1:time_length,1))-PIT,"Color","m")
% plot(ann(time),ann(globalMean(T_arise_ann_control)+preset_T_noise(1:time_length,1))-PIT,"Color","b")

% plot(ann(time),ann(globalMean(T_no_inj)+preset_T_noise(1:time_length,1))-PIT,"Color","r")
plot(ann(time),ann(globalMean(T_no_inj_worse)+preset_T_noise(1:time_length,1))-PIT,"Color","k")

% plot(ann(time),ann(globalMean(T_no_inj_worse)+preset_T_noise(1:time_length,1))-PIT)

ylabel("Temperature, °C above Pre-Industrial")
% ylabel("Temperature, K")
xlabel("Year")
% legend("Global Warming (Business-As-Usual)","Global Warming (Mitigated)","Stratopsheric Aerosol Injection (beginning in 2035)","Location","nw")
legend("Uncontrolled Greenhouse Gas Emissions","Mitigated Greenhouse Gas Emissions","Stratopsheric Aerosol Injection (beginning in 2035)","1.5°C Paris Agreement Target","Location","nw")
legend("Uncontrolled Greenhouse Gas Emissions","Mitigated Greenhouse Gas Emissions","1.5°C Paris Agreement Target","Location","nw")
legend("Uncontrolled Greenhouse Gas Emissions","1.5°C Paris Agreement Target","Location","nw")
title("Near-Surface Global Mean Temperature (CIDER trained on E3SM)")
% title("Temperature of a Nuclear-Fusion-Powered Ballistic Projectile")
% %% Plot Figure 2
% PIT = 287;
% figure
% hold on
% box on
% grid on
% set(gca,'DefaultLineLineWidth',4)
% % plot(ann(time),ann(globalMean(T_no_inj_worse)+preset_T_noise(1:time_length,1))-PIT,"Color","k")
% 
% plot(ann(time),ann(globalMean(T_no_inj)+preset_T_noise(1:time_length,1))-PIT,"Color","r")
% 
% plot(ann(time),ann(globalMean(T_arise_ann_control)+preset_T_noise(1:time_length,1))-PIT,"Color","b")
% plot(ann(time),ann(globalMean(T_arise_ann_control_370)+preset_T_noise(1:time_length,1))-PIT,"Color",[89, 161, 217]/255)
% plot(ann(time),ann(globalMean(T_term)+preset_T_noise(1:time_length,1))-PIT,"Color",[168, 94, 50]/255)
% plot(ann(time),ann(globalMean(T_term_370)+preset_T_noise(1:time_length,1))-PIT,"Color","m")
% plot(ann(time),ann(globalMean(T_arise_ann_control_370)+preset_T_noise(1:time_length,1))-PIT,"Color",[89, 161, 217]/255)
% 
% plot(ann(time),ann(globalMean(T_arise_ann_control)+preset_T_noise(1:time_length,1))-PIT,"Color","b")
% plot(ann(time),ann(globalMean(T_no_inj)+preset_T_noise(1:time_length,1))-PIT,"Color","r")
% 
% % 
% % plot(ann(time),ann(globalMean(T_no_inj)+preset_T_noise(1:time_length,1))-PIT,"Color","r")
% 
% % plot(ann(time),ann(globalMean(T_no_inj_worse)+preset_T_noise(1:time_length,1))-PIT)
% 
% ylabel({"Near-Surface Global Mean Temperature,", "°C above Pre-Industrial"})
% % ylabel("Temperature, K")
% xlabel("Year")
% % legend("Global Warming (Business-As-Usual)","Global Warming (Mitigated)","Stratopsheric Aerosol Injection (beginning in 2035)","Location","nw")
% legend("Mitigated Emissions, No SAI","SAI (Mitigated Emissions)","SAI (Uncontrolled Emissions)","Terminated SAI (Mitigated Emissions) (Bad)","Terminated SAI (Uncontrolled Emissions) (Very Bad)","Location","nw")
% % legend("Mitigated Emissions, No SAI","SAI (Mitigated Emissions)","SAI (Uncontrolled Emissions)","Location","nw")
% title("SAI under Different Emissions")
% % title("Temperature of a Nuclear-Fusion-Powered Ballistic Projectile")
