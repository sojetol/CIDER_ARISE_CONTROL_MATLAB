%% User Variables

% Target
refvals = [288.5, 0.8767, -5.89]; % Where controller is bringing climate to
%           T0       T1      T2

% Other User Variables
total_years = 35;
tseries = 0:1/12:total_years-1/12;
noise_on = 1;
offset = (randi(300)-1)*12;
offset = 100;
my_title = "Emulating with T0-T1-T2 (Different) PreInd Variability, Plotted annual";
years_range = 1+offset:12*total_years+offset;
load("Tools\PI_T0T1T2.mat")
preset_T_noise = 1*[1*(PI_T0(years_range)-mean(PI_T0(years_range))),...
                    1*(PI_T1(years_range)-mean(PI_T1(years_range))),...
                    1*(PI_T2(years_range)-mean(PI_T2(years_range)))];
                   
preset_T0_noise = 0*.35*(randn(total_years*12,1))+1*1.75*sin(2*pi*tseries)';
preset_T1_noise = 0*.35*(randn(total_years*12,1))+1*3.27*sin(2*pi*tseries)';
% preset_T_noise(:,1) = preset_T0_noise;
preset_T_noise = preset_T_noise*noise_on;
ann_on = 1;
CO2_ppm_growth = 3;

%% Initialize other parameters and arrays
add_all_paths_Sandia
load CESM_params
log_array = [];
q_array = zeros(4,total_years);

%% CO2
% load CO2_concentrations.mat
CO2_SSP245 = 400+(0:total_years+21)'*CO2_ppm_growth;
CO2_ref = CO2_SSP245(20);
CO2_forcing_SSP245 = 5.35*log((CO2_SSP245(21:end)+0)/CO2_ref);
CO2_forcing_SSP245_month = repeatElements(CO2_forcing_SSP245,12);
CO2_SSP245_smooth = linspace(CO2_SSP245(21)-CO2_ppm_growth/2,CO2_SSP245(end)+CO2_ppm_growth/2,length(CO2_forcing_SSP245_month))';
CO2_forcing_SSP245_month_smooth =  5.35*log((CO2_SSP245_smooth)/CO2_ref);
CO2_forcing_SSP245_month = CO2_forcing_SSP245_month_smooth;

%% Simulate
% Just to initialize controls
i_year = 1;
time_length = 12*i_year;
inj_array = [zeros(time_length,1) repeatElements(q_array(4,1:i_year)'/12,12) repeatElements(q_array(3,1:i_year)'/12,12)  zeros(time_length,1) repeatElements(q_array(2,1:i_year)'/12,12) repeatElements(q_array(1,1:i_year)'/12,12)  zeros(time_length,1)];
T_arise_ann_control = CIDER_pattern_from_all_injections_and_CO2([ inj_array CO2_forcing_SSP245_month(1:time_length)],param_AOD_all,param_T_all,pattern_T_all)+pattern_T_base;
[q_array(:,i_year),log_array] = ARISE_controller_translation(T_arise_ann_control,log_array,refvals,preset_T_noise(1:12,:));
% Simulate (annual controls)
for i_year = 1:total_years
time_length = 12*i_year;
inj_array = [zeros(time_length,1) repeatElements(q_array(4,1:i_year)'/12,12) repeatElements(q_array(3,1:i_year)'/12,12)  zeros(time_length,1) repeatElements(q_array(2,1:i_year)'/12,12) repeatElements(q_array(1,1:i_year)'/12,12)  zeros(time_length,1)];

T_arise_ann_control = CIDER_pattern_from_all_injections_and_CO2([ inj_array CO2_forcing_SSP245_month(1:time_length)],param_AOD_all,param_T_all,pattern_T_all)+pattern_T_base;

[q_array(:,i_year+1),log_array] = ARISE_controller_translation(T_arise_ann_control,log_array,refvals,preset_T_noise(1:time_length,:));

end
T_no_inj = CIDER_pattern_from_all_injections_and_CO2([inj_array*0 CO2_forcing_SSP245_month(1:time_length)],param_AOD_all,param_T_all,pattern_T_all)+pattern_T_base;
%% Simulate montly
% Just to initialize controls
q_array_month = zeros(4,total_years*12);
i_month = 1;
time_length = i_month;
log_array_monthly = [];
inj_array_monthly = [zeros(time_length,1) q_array_month(4,1:i_month)'/12 q_array_month(3,1:i_month)'/12  zeros(time_length,1) q_array_month(2,1:i_month)'/12 q_array_month(1,1:i_month)'/12  zeros(time_length,1)];
T_arise_monthly_control = CIDER_pattern_from_all_injections_and_CO2([ inj_array_monthly CO2_forcing_SSP245_month(1:time_length)],param_AOD_all,param_T_all,pattern_T_all)+pattern_T_base;
[q_array_month(:,i_month),log_array_monthly] = ARISE_controller_translation_monthly(T_arise_monthly_control,log_array_monthly,refvals,preset_T_noise(1:time_length,:));

q_array_month_LYM_filt = zeros(4,total_years*12);
i_month = 1;
time_length = i_month;
log_array_monthly_LYM_filt = [];
inj_array_lym_filt = [zeros(time_length,1) q_array_month_LYM_filt(4,1:i_month)'/12 q_array_month_LYM_filt(3,1:i_month)'/12  zeros(time_length,1) q_array_month_LYM_filt(2,1:i_month)'/12 q_array_month_LYM_filt(1,1:i_month)'/12  zeros(time_length,1)];
T_arise_monthly_control_LYM_filt = CIDER_pattern_from_all_injections_and_CO2([ inj_array_lym_filt CO2_forcing_SSP245_month(1:time_length)],param_AOD_all,param_T_all,pattern_T_all)+pattern_T_base;
[q_array_month_LYM_filt(:,i_month),log_array_monthly_LYM_filt] = ARISE_controller_translation_monthly(T_arise_monthly_control_LYM_filt,log_array_monthly_LYM_filt,refvals,preset_T_noise(1:time_length,:));

% Simulate (monthly controls)
for i_month = 1:total_years*12
time_length = i_month;


X = sprintf('Progress: %d / %d',time_length,total_years*12);
disp(X)

inj_array_monthly = [zeros(time_length,1) q_array_month(4,1:i_month)'/12 q_array_month(3,1:i_month)'/12 zeros(time_length,1) q_array_month(2,1:i_month)'/12 q_array_month(1,1:i_month)'/12 zeros(time_length,1)];

T_arise_monthly_control = CIDER_pattern_from_all_injections_and_CO2([ inj_array_monthly CO2_forcing_SSP245_month(1:time_length)],param_AOD_all,param_T_all,pattern_T_all)+pattern_T_base;

[q_array_month(:,i_month+1),log_array_monthly] = ARISE_controller_translation_monthly(T_arise_monthly_control,log_array_monthly,refvals,preset_T_noise(1:time_length,:));


inj_array_lym_filt = [zeros(time_length,1) q_array_month_LYM_filt(4,1:i_month)'/12 q_array_month_LYM_filt(3,1:i_month)'/12 zeros(time_length,1) q_array_month_LYM_filt(2,1:i_month)'/12 q_array_month_LYM_filt(1,1:i_month)'/12 zeros(time_length,1)];

T_arise_monthly_control_LYM_filt = CIDER_pattern_from_all_injections_and_CO2([ inj_array_lym_filt CO2_forcing_SSP245_month(1:time_length)],param_AOD_all,param_T_all,pattern_T_all)+pattern_T_base;

[q_array_month_LYM_filt(:,i_month+1),log_array_monthly_LYM_filt] = ARISE_controller_translation_monthly_w_filter(T_arise_monthly_control_LYM_filt,log_array_monthly_LYM_filt,refvals,preset_T_noise(1:time_length,:));

end
% T_no_inj = CIDER_pattern_from_all_injections_and_CO2([inj_array*0 CO2_forcing_SSP245_month(1:time_length)],param_AOD_all,param_T_all,pattern_T_all)+pattern_T_base;
%% Calculate AOD
AOD_annual = globalMean(CIDER_AOD_pattern_from_all_injections(param_AOD_all,inj_array,pattern_AOD_all)+pattern_AOD_base);
AOD_monthly = globalMean(CIDER_AOD_pattern_from_all_injections(param_AOD_all,inj_array_monthly,pattern_AOD_all)+pattern_AOD_base);
AOD_monthly_LYM = globalMean(CIDER_AOD_pattern_from_all_injections(param_AOD_all,inj_array_lym_filt,pattern_AOD_all)+pattern_AOD_base);

%% Simulate monthly with filter
% Just to initialize controls
% q_array_month_LYM_filt = zeros(4,total_years*12);
% i_month = 1;
% time_length = i_month;
% log_array_monthly_LYM_filt = [];
% inj_array_lym_filt = [zeros(time_length,1) q_array_month_LYM_filt(4,1:i_month)'/12 q_array_month_LYM_filt(3,1:i_month)'/12  zeros(time_length,1) q_array_month_LYM_filt(2,1:i_month)'/12 q_array_month_LYM_filt(1,1:i_month)'/12  zeros(time_length,1)];
% T_arise_monthly_control_LYM_filt = CIDER_pattern_from_all_injections_and_CO2([ inj_array_lym_filt CO2_forcing_SSP245_month(1:time_length)],param_AOD_all,param_T_all,pattern_T_all)+pattern_T_base;
% [q_array_month_LYM_filt(:,i_month),log_array_monthly_LYM_filt] = ARISE_controller_translation_monthly(T_arise_monthly_control_LYM_filt,log_array_monthly_LYM_filt,refvals,preset_noise(1:time_length));
% % Simulate (annual controls)
% for i_month = 1:total_years*12
% time_length = i_month
% inj_array_lym_filt = [zeros(time_length,1) q_array_month_LYM_filt(4,1:i_month)'/12 q_array_month_LYM_filt(3,1:i_month)'/12 zeros(time_length,1) q_array_month_LYM_filt(2,1:i_month)'/12 q_array_month_LYM_filt(1,1:i_month)'/12 zeros(time_length,1)];
% T_arise_monthly_control_LYM_filt = CIDER_pattern_from_all_injections_and_CO2([ inj_array_lym_filt CO2_forcing_SSP245_month(1:time_length)],param_AOD_all,param_T_all,pattern_T_all)+pattern_T_base;
% [q_array_month_LYM_filt(:,i_month+1),log_array_monthly_LYM_filt] = ARISE_controller_translation_monthly(T_arise_monthly_control_LYM_filt,log_array_monthly_LYM_filt,refvals,preset_noise(1:time_length));
% end
% T_no_inj = CIDER_pattern_from_all_injections_and_CO2([inj_array*0 CO2_forcing_SSP245_month(1:time_length)],param_AOD_all,param_T_all,pattern_T_all)+pattern_T_base;

%%
time = annualToMonthly(2035:2035+total_years-1);
figure
tiley = tiledlayout(3,2);
title(tiley,my_title)
nexttile(2,[1 1])
set(gca,'DefaultLineLineWidth',2)

hold on
box on
grid on
if ann_on == 1
plot(ann(time),ann(globalMean(T_no_inj)+preset_T_noise(:,1)),"Color","k")
plot(ann(time),ann(globalMean(T_arise_ann_control)+preset_T_noise(:,1)),"Color",'b')
plot(ann(time),ann(globalMean(T_arise_monthly_control)+preset_T_noise(:,1)),"Color",'g')
plot(ann(time),ann(globalMean(T_arise_monthly_control_LYM_filt)+preset_T_noise(:,1)),"Color",'c')
else
plot((time),(globalMean(T_no_inj)+preset_T_noise(:,1)),"Color","k")
plot((time),(globalMean(T_arise_ann_control)+preset_T_noise(:,1)),"Color",'b')
plot((time),(globalMean(T_arise_monthly_control)+preset_T_noise(:,1)),"Color",'g')
plot((time),(globalMean(T_arise_monthly_control_LYM_filt)+preset_T_noise(:,1)),"Color",'c')
end


line([time(1) time(end)],[refvals(1) refvals(1)],"Color","r")

set(gca,'DefaultLineLineWidth',4)
plot(time,globalMean(T_no_inj),"Color",'k')
plot(time,globalMean(T_arise_ann_control),"Color",'b')
plot(time,globalMean(T_arise_monthly_control),"Color",'g')
plot(time,globalMean(T_arise_monthly_control_LYM_filt),"Color",'c')
% line([time(1) time(end)],[refvals(1) refvals(1)],"Color","r")
title("Global Mean Temperature, CIDER (CESM2)")
ylabel("T0 (K)")

monthly_error = sqrt(mean((globalMean(T_arise_monthly_control)+preset_T_noise(:,1)-refvals(1)).^2))
annual_error = sqrt(mean((globalMean(T_arise_ann_control)+preset_T_noise(:,1)-refvals(1)).^2))
monthly_lym_error = sqrt(mean((globalMean(T_arise_monthly_control_LYM_filt)+preset_T_noise(:,1)-refvals(1)).^2))
legend("No SAI","ARISE (Annual)","ARISE (Monthly control)","ARISE (Monthly with LYM filter)","Target","Location","nw")
xlim([time(1) time(end)])

nexttile(4)
set(gca,'DefaultLineLineWidth',2)

hold on
box on
grid on

if ann_on == 1
plot(ann(time),ann(calculateT1(T_no_inj)+preset_T_noise(:,2)),"Color",'k')
plot(ann(time),ann(calculateT1(T_arise_ann_control)+preset_T_noise(:,2)),"Color",'b')
plot(ann(time),ann(calculateT1(T_arise_monthly_control)+preset_T_noise(:,2)),"Color",'g')
plot(ann(time),ann(calculateT1(T_arise_monthly_control_LYM_filt)+preset_T_noise(:,2)),"Color",'c')
else
plot(time,calculateT1(T_no_inj)+preset_T_noise(:,2),"Color",'k')
plot(time,calculateT1(T_arise_ann_control)+preset_T_noise(:,2),"Color",'b')
plot(time,calculateT1(T_arise_monthly_control)+preset_T_noise(:,2),"Color",'g')
plot(time,calculateT1(T_arise_monthly_control_LYM_filt)+preset_T_noise(:,2),"Color",'c')
end
set(gca,'DefaultLineLineWidth',4)

plot(time,calculateT1(T_no_inj),"Color",'k')
plot(time,calculateT1(T_arise_ann_control),"Color",'b')
plot(time,calculateT1(T_arise_monthly_control),"Color",'g')
plot(time,calculateT1(T_arise_monthly_control_LYM_filt),"Color",'c')
line([time(1) time(end)],[refvals(2) refvals(2)],"Color",'r')
xlim([time(1) time(end)])
ylabel("T1 (K/degree lat)")
title("N-S T Gradient")
nexttile(6)
hold on
box on
grid on
set(gca,'DefaultLineLineWidth',2)


if ann_on ==1
plot(ann(time),ann(calculateT2(T_no_inj)+preset_T_noise(:,3)),"Color",'k')
plot(ann(time),ann(calculateT2(T_arise_ann_control)+preset_T_noise(:,3)),"Color",'b')
plot(ann(time),ann(calculateT2(T_arise_monthly_control)+preset_T_noise(:,3)),"Color",'g')
plot(ann(time),ann(calculateT2(T_arise_monthly_control_LYM_filt)+preset_T_noise(:,3)),"Color",'c')
else
plot(time,calculateT2(T_no_inj)+preset_T_noise(:,3),"Color",'k')
plot(time,calculateT2(T_arise_ann_control)+preset_T_noise(:,3),"Color",'b')
plot(time,calculateT2(T_arise_monthly_control)+preset_T_noise(:,3),"Color",'g')
plot(time,calculateT2(T_arise_monthly_control_LYM_filt)+preset_T_noise(:,3),"Color",'c')
end
set(gca,'DefaultLineLineWidth',4)

plot(time,calculateT2(T_no_inj),"Color",'k')
plot(time,calculateT2(T_arise_ann_control),"Color",'b')
plot(time,calculateT2(T_arise_monthly_control),"Color",'g')
plot(time,calculateT2(T_arise_monthly_control_LYM_filt),"Color",'c')
line([time(1) time(end)],[refvals(3) refvals(3)],"Color",'r')
xlabel("Year")
ylabel("T2 (K/?)")
title("Equator-to-pole T Gradient")
xlim([time(1) time(end)])

%%
% figure
nexttile(1)
hold on
box on 
grid on
set(gca,'DefaultLineLineWidth',2)

plot((1:total_years+1)-1,q_array(1,1:end)','Color','r')
plot((1:total_years+1)-1,q_array(2,1:end)','Color','g')
plot((1:total_years+1)-1,q_array(3,1:end)','Color','b')
plot((1:total_years+1)-1,q_array(4,1:end)','Color','c')


if ann_on == 1

plot(ann(annualToMonthly(1:total_years)-1), ann(q_array_month(1,1:end-1)),'Color','r','LineStyle',":")
plot(ann(annualToMonthly(1:total_years)-1), ann(q_array_month(2,1:end-1)),'Color','g','LineStyle',":")
plot(ann(annualToMonthly(1:total_years)-1), ann(q_array_month(3,1:end-1)),'Color','b','LineStyle',":")
plot(ann(annualToMonthly(1:total_years)-1),ann( q_array_month(4,1:end-1)),'Color','c','LineStyle',":")

else
plot(annualToMonthly(1:total_years)-1, q_array_month(1,1:end-1),'Color','r','LineStyle',":")
plot(annualToMonthly(1:total_years)-1, q_array_month(2,1:end-1),'Color','g','LineStyle',":")
plot(annualToMonthly(1:total_years)-1, q_array_month(3,1:end-1),'Color','b','LineStyle',":")
plot(annualToMonthly(1:total_years)-1, q_array_month(4,1:end-1),'Color','c','LineStyle',":")
end

plot(annualToMonthly(1:total_years)-1, q_array_month_LYM_filt(1,1:end-1),'Color','r','LineStyle',"--")
plot(annualToMonthly(1:total_years)-1, q_array_month_LYM_filt(2,1:end-1),'Color','g','LineStyle',"--")
plot(annualToMonthly(1:total_years)-1, q_array_month_LYM_filt(3,1:end-1),'Color','b','LineStyle',"--")
plot(annualToMonthly(1:total_years)-1, q_array_month_LYM_filt(4,1:end-1),'Color','c','LineStyle',"--")

legend('30S','15S','15N','30N','Location','nw')
xlabel("Year")
ylabel("Injection Rate, TgSO2")
title("Controls (Annual: Solid, Monthly: Dotted, Monthly w/ LYM filter: Dashed)")
%%
% figure
% % plot(time,globalMean(T_no_inj),"Color",'k')
% hold on
% box on
% grid on
% plot(time,globalMean(T_arise_monthly_control)-globalMean(T_arise_ann_control),"Color",'g','LineWidth',4)
% plot(time,globalMean(T_arise_monthly_control_LYM_filt)-globalMean(T_arise_ann_control),"Color",'c','LineWidth',4)
% legend("Difference between Annual and Monthly","Difference between Annual and Monthly w/ LYM filter")
% ylabel("Temperature, K")
% xlabel("Year")
% xlim([time(1) time(end)])
% 
% title("Controlling ARISE")

%%
% figure
nexttile(5)
hold on
box on
grid on
set(gca,'DefaultLineLineWidth',2)
plot(time,globalMean(T_arise_ann_control)-refvals(1),"Color",'b')
plot(time,globalMean(T_arise_monthly_control)-refvals(1),"Color",'g')
plot(time,globalMean(T_arise_monthly_control_LYM_filt)-refvals(1),"Color",'c')
ylabel("Temperature Error from Target, K")
title("Forced Temperature - Temperature Target")
xlabel("Year")
xlim([time(1) time(end)])

legend("Annual","Monthly control","Monthly with LYM filter","Location","ne")

% line([time(1) time(end)],[refvals(1) refvals(1)],"Color","r")
%%
% figure
nexttile(3)
hold on
grid on
box on 
set(gca,'DefaultLineLineWidth',2)

plot(time,AOD_annual,"Color",'b')
plot(time,AOD_monthly,"Color",'g')
plot(time,AOD_monthly_LYM,"Color",'c')
xlim([time(1) time(end)])

title('AOD')
ylabel("AOD")
xlabel("Year")
legend("Annual","Monthly control","Monthly with LYM filter","Location","nw")
