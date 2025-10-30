add_all_paths_Sandia
load CESM_params.mat
AODVIS_30N_DJF = ncread("Climate_data\AODVIS_30N_DJF_201501_201912_readable.nc","AODVIS");
AODVIS_30N_SON = ncread("Climate_data\AODVIS_30N_SON_201501_201912_readable.nc","AODVIS");
AODVIS_30N_MAM = ncread("Climate_data\AODVIS_30N_MAM_201501_201912_readable.nc","AODVIS");
AODVIS_30N_JJA = ncread("Climate_data\AODVIS_30N_JJA_201501_201912_readable.nc","AODVIS");
AODVIS_30N_DJF_12Tg = ncread("Climate_data\AODVIS_30N_DJF_12Tg_201501_201912_readable.nc","AODVIS");
AODVIS_30N_SON_12Tg = ncread("Climate_data\AODVIS_30N_SON_12Tg_201501_201912_readable.nc","AODVIS");
AODVIS_30N_MAM_12Tg = ncread("Climate_data\AODVIS_30N_MAM_12Tg_201501_201912_readable.nc","AODVIS");
AODVIS_30N_JJA_12Tg = ncread("Climate_data\AODVIS_30N_JJA_12Tg_201501_201912_readable.nc","AODVIS");
AODVIS_SSP245 = ncread("Climate_data\AODVIS_SSP245_201501_201912_readable.nc","AODVIS");

TREFHT_30N_DJF = ncread("Climate_data\TREFHT_30N_DJF_201501_201912_readable.nc","TREFHT");
TREFHT_30N_SON = ncread("Climate_data\TREFHT_30N_SON_201501_201912_readable.nc","TREFHT");
TREFHT_30N_MAM = ncread("Climate_data\TREFHT_30N_MAM_201501_201912_readable.nc","TREFHT");
TREFHT_30N_JJA = ncread("Climate_data\TREFHT_30N_JJA_201501_201912_readable.nc","TREFHT");
TREFHT_30N_DJF_12Tg = ncread("Climate_data\TREFHT_30N_DJF_12Tg_201501_201912_readable.nc","TREFHT");
TREFHT_30N_SON_12Tg = ncread("Climate_data\TREFHT_30N_SON_12Tg_201501_201912_readable.nc","TREFHT");
TREFHT_30N_MAM_12Tg = ncread("Climate_data\TREFHT_30N_MAM_12Tg_201501_201912_readable.nc","TREFHT");
TREFHT_30N_JJA_12Tg = ncread("Climate_data\TREFHT_30N_JJA_12Tg_201501_201912_readable.nc","TREFHT");
TREFHT_SSP245 = ncread("Climate_data\TREFHT_SSP245_201501_201912_readable.nc","TREFHT");


AODVIS_15N_DJF = ncread("Climate_data\AODVIS_15N_DJF_201501_201912_readable.nc","AODVIS");
AODVIS_15N_SON = ncread("Climate_data\AODVIS_15N_SON_201501_201912_readable.nc","AODVIS");
AODVIS_15N_MAM = ncread("Climate_data\AODVIS_15N_MAM_201501_201912_readable.nc","AODVIS");
AODVIS_15N_JJA = ncread("Climate_data\AODVIS_15N_JJA_201501_201912_readable.nc","AODVIS");

TREFHT_15N_DJF = ncread("Climate_data\TREFHT_15N_DJF_201501_201912_readable.nc","TREFHT");
TREFHT_15N_SON = ncread("Climate_data\TREFHT_15N_SON_201501_201912_readable.nc","TREFHT");
TREFHT_15N_MAM = ncread("Climate_data\TREFHT_15N_MAM_201501_201912_readable.nc","TREFHT");
TREFHT_15N_JJA = ncread("Climate_data\TREFHT_15N_JJA_201501_201912_readable.nc","TREFHT");
index = 1;



%%
MAM_dif_30 = AODVIS_30N_MAM(index,3:end)-AODVIS_SSP245(index,3:end);
JJA_dif_30 = AODVIS_30N_JJA(index,6:end)-AODVIS_SSP245(index,6:end);
SON_dif_30 = AODVIS_30N_SON(index,9:end)-AODVIS_SSP245(index,9:end);
DJF_dif_30 = AODVIS_30N_DJF(index,12:end)-AODVIS_SSP245(index,12:end);
AOD_VIS_30_all = double((MAM_dif_30(1:49)+JJA_dif_30(1:49)+SON_dif_30(1:49)+DJF_dif_30)/4)';
MAM_dif_15 = AODVIS_15N_MAM(index,3:end)-AODVIS_SSP245(index,3:end);
JJA_dif_15 = AODVIS_15N_JJA(index,6:end)-AODVIS_SSP245(index,6:end);
SON_dif_15 = AODVIS_15N_SON(index,9:end)-AODVIS_SSP245(index,9:end);
DJF_dif_15 = AODVIS_15N_DJF(index,12:end)-AODVIS_SSP245(index,12:end);
AOD_VIS_15_all = double((MAM_dif_15(1:49)+JJA_dif_15(1:49)+SON_dif_15(1:49)+DJF_dif_15)/4)';% plot(AODVIS_30N_JJA(index,6:end)-AODVIS_SSP245(index,6:end))
% plot(AODVIS_30N_SON(index,9:end)-AODVIS_SSP245(index,9:end))
% plot(AODVIS_30N_DJF(index,12:end)-AODVIS_SSP245(index,12:end))
% param_T_all(:,1) = 200*param_T_all(:,1)
% param_T_all(:,2) = param_T_all(:,2)/20

amount_per_month = 2;
all_pulse_injections_OG = double(AOD_VIS_30_all*0);
all_pulse_injections_OG(1:3) = all_pulse_injections_OG(1:3)+amount_per_month;
all_pulse_injections_longer  = double(AOD_VIS_30_all*0);
all_pulse_injections_longer(3:5) = all_pulse_injections_longer(3:5)+amount_per_month;
SO2_params = [0.2689 0.2405  0];
SO4_in =  CIDER_AOD_from_injection(SO2_params,all_pulse_injections_OG);
% SO4_in =  CIDER_AOD_from_injection([0.5 0.5  0],all_pulse_injections_OG);
all_pulse_injections = SO4_in;
all_pulse_injections_longer = CIDER_AOD_from_injection(SO2_params,all_pulse_injections_longer);
all_pulse_injections_taller = all_pulse_injections*2;

figure 
plot(SO4_in)
weightings = [1 0];
all_pulse_responses = AOD_VIS_30_all;
all_pulse_responses2 = AOD_VIS_15_all;
paramsE3SM = CIDER_train_AOD_params_seasonal(all_pulse_injections,all_pulse_responses,[],[],weightings);
paramsE3SM2 = CIDER_train_AOD_params_seasonal(all_pulse_injections,all_pulse_responses2,[],[],weightings);
emu_AOD = CIDER_AOD_from_injection(paramsE3SM,all_pulse_injections);
emu_AOD2 = CIDER_AOD_from_injection(paramsE3SM2,all_pulse_injections);
emu_AOD_taller = CIDER_AOD_from_injection(paramsE3SM,all_pulse_injections_taller);
emu_AOD_longer = CIDER_AOD_from_injection(paramsE3SM,all_pulse_injections_longer);

emu_AOD_CESM = CIDER_AOD_from_injection(param_AOD_all(2,:),all_pulse_injections);
emu_AOD_step = CIDER_AOD_from_injection(paramsE3SM,CIDER_AOD_from_injection([0.2689 0.2405  0],all_pulse_injections*0+1));
emu_AOD_CESM_step = CIDER_AOD_from_injection(param_AOD_all(2,:),all_pulse_injections*0+1);
%%
tseries = annualToMonthly(2015:2019);
index = 1;
figure
tiley = tiledlayout(2,2);
nexttile(1)
set(gca,'DefaultLineLineWidth',2)
hold on
box on
grid on
plot(tseries,AODVIS_15N_MAM(index,:)-AODVIS_SSP245(index,:))
plot(tseries,AODVIS_15N_JJA(index,:)-AODVIS_SSP245(index,:))
plot(tseries,AODVIS_15N_SON(index,:)-AODVIS_SSP245(index,:))
plot(tseries,AODVIS_15N_DJF(index,:)-AODVIS_SSP245(index,:))
title("AODVIS, 6Tg@15N - SSP245")
legend("MAM","JJA","SON","DJF")
ylabel("AOD Difference")
xlabel("Year")
nexttile(3)
set(gca,'DefaultLineLineWidth',2)
hold on
box on
grid on
plot(tseries,TREFHT_15N_MAM(index,:)-TREFHT_SSP245(index,:))
plot(tseries,TREFHT_15N_JJA(index,:)-TREFHT_SSP245(index,:))
plot(tseries,TREFHT_15N_SON(index,:)-TREFHT_SSP245(index,:))
plot(tseries,TREFHT_15N_DJF(index,:)-TREFHT_SSP245(index,:))
ylabel("Temperature Difference (K)")
ylim([-.8 .5])
xlabel("Year")
title("TREFHT, 6Tg@15N - SSP245")
% legend("MAM","JJA","SON","DJF")

nexttile(2)
set(gca,'DefaultLineLineWidth',2)
hold on
box on
grid on
plot(tseries,AODVIS_30N_MAM(index,:)-AODVIS_SSP245(index,:))
plot(tseries,AODVIS_30N_JJA(index,:)-AODVIS_SSP245(index,:))
plot(tseries,AODVIS_30N_SON(index,:)-AODVIS_SSP245(index,:))
plot(tseries,AODVIS_30N_DJF(index,:)-AODVIS_SSP245(index,:))
ylabel("AOD Difference")
xlabel("Year")
title("AODVIS, 6Tg@30N - SSP245")
% legend("MAM","JJA","SON","DJF")

nexttile(4)
set(gca,'DefaultLineLineWidth',2)
hold on
box on
grid on
plot(tseries,TREFHT_30N_MAM(index,:)-TREFHT_SSP245(index,:))
plot(tseries,TREFHT_30N_JJA(index,:)-TREFHT_SSP245(index,:))
plot(tseries,TREFHT_30N_SON(index,:)-TREFHT_SSP245(index,:))
plot(tseries,TREFHT_30N_DJF(index,:)-TREFHT_SSP245(index,:))
ylabel("Temperature Difference (K)")
xlabel("Year")
ylim([-.8 .5])
title("TREFHT, 6Tg@30N - SSP245")
% legend("MAM","JJA","SON","DJF")
%% Plot rezeroed

index = 1;
tseries = annualToMonthly(0:4);
figure
tiley = tiledlayout(2,2);
nexttile(1)
set(gca,'DefaultLineLineWidth',2)
hold on
box on
grid on
MAM_dif = AODVIS_15N_MAM(index,3:end)-AODVIS_SSP245(index,3:end);
JJA_dif = AODVIS_15N_JJA(index,6:end)-AODVIS_SSP245(index,6:end);
SON_dif = AODVIS_15N_SON(index,9:end)-AODVIS_SSP245(index,9:end);
DJF_dif = AODVIS_15N_DJF(index,12:end)-AODVIS_SSP245(index,12:end);
plot(AODVIS_15N_MAM(index,3:end)-AODVIS_SSP245(index,3:end))
plot(AODVIS_15N_JJA(index,6:end)-AODVIS_SSP245(index,6:end))
plot(AODVIS_15N_SON(index,9:end)-AODVIS_SSP245(index,9:end))
plot(AODVIS_15N_DJF(index,12:end)-AODVIS_SSP245(index,12:end))
plot(emu_AOD2,"LineWidth",4)
plot((MAM_dif(1:49)+JJA_dif(1:49)+SON_dif(1:49)+DJF_dif)/4, "LineWidth",4)

title("AODVIS, 6Tg@15N - SSP245")
legend("MAM","JJA","SON","DJF","CIDER*","Mean")
ylabel("AOD Difference")
xlabel("Months after initiation")
nexttile(3)
set(gca,'DefaultLineLineWidth',2)
hold on
box on
grid on
MAM_dif = TREFHT_15N_MAM(index,3:end)-TREFHT_SSP245(index,3:end);
JJA_dif = TREFHT_15N_JJA(index,6:end)-TREFHT_SSP245(index,6:end);
SON_dif = TREFHT_15N_SON(index,9:end)-TREFHT_SSP245(index,9:end);
DJF_dif = TREFHT_15N_DJF(index,12:end)-TREFHT_SSP245(index,12:end);
plot(TREFHT_15N_MAM(index,3:end)-TREFHT_SSP245(index,3:end))
plot(TREFHT_15N_JJA(index,6:end)-TREFHT_SSP245(index,6:end))
plot(TREFHT_15N_SON(index,9:end)-TREFHT_SSP245(index,9:end))
plot(TREFHT_15N_DJF(index,12:end)-TREFHT_SSP245(index,12:end))
emu_CIDER_15 = CIDER_response_from_1_forcing(param_T_all(3,:),emu_AOD2);
plot(emu_CIDER_15,"LineWidth",4)
plot((MAM_dif(1:49)+JJA_dif(1:49)+SON_dif(1:49)+DJF_dif)/4, "LineWidth",4)

ylabel("Temperature Difference (K)")
ylim([-.8 .5])
xlabel("Months after initiation")
title("TREFHT, 6Tg@15N - SSP245")
% legend("MAM","JJA","SON","DJF")

nexttile(2)
set(gca,'DefaultLineLineWidth',2)
hold on
box on
grid on
MAM_dif = AODVIS_30N_MAM(index,3:end)-AODVIS_SSP245(index,3:end);
JJA_dif = AODVIS_30N_JJA(index,6:end)-AODVIS_SSP245(index,6:end);
SON_dif = AODVIS_30N_SON(index,9:end)-AODVIS_SSP245(index,9:end);
DJF_dif = AODVIS_30N_DJF(index,12:end)-AODVIS_SSP245(index,12:end);
plot(AODVIS_30N_MAM(index,3:end)-AODVIS_SSP245(index,3:end))
plot(AODVIS_30N_JJA(index,6:end)-AODVIS_SSP245(index,6:end))
plot(AODVIS_30N_SON(index,9:end)-AODVIS_SSP245(index,9:end))
plot(AODVIS_30N_DJF(index,12:end)-AODVIS_SSP245(index,12:end))
plot(emu_AOD,"LineWidth",4)
plot((MAM_dif(1:49)+JJA_dif(1:49)+SON_dif(1:49)+DJF_dif)/4, "LineWidth",4)
% plot(emu_AOD_taller,"LineWidth",4)
% % plot(all_pulse_injections/2/10)
% % plot(emu_AOD_CESM,"LineWidth",4,"LineStyle",":")
% plot(emu_AOD_step,"LineWidth",4)
% plot(emu_AOD_CESM_step,"LineWidth",4,"LineStyle",":")
ylabel("AOD Difference")
xlabel("Months after initiation")
title("AODVIS, 6Tg@30N - SSP245")
% legend("MAM","JJA","SON","DJF")

nexttile(4)
set(gca,'DefaultLineLineWidth',2)
hold on
box on
grid on
MAM_dif = TREFHT_30N_MAM(index,3:end)-TREFHT_SSP245(index,3:end);
JJA_dif = TREFHT_30N_JJA(index,6:end)-TREFHT_SSP245(index,6:end);
SON_dif = TREFHT_30N_SON(index,9:end)-TREFHT_SSP245(index,9:end);
DJF_dif = TREFHT_30N_DJF(index,12:end)-TREFHT_SSP245(index,12:end);
plot(TREFHT_30N_MAM(index,3:end)-TREFHT_SSP245(index,3:end))
plot(TREFHT_30N_JJA(index,6:end)-TREFHT_SSP245(index,6:end))
plot(TREFHT_30N_SON(index,9:end)-TREFHT_SSP245(index,9:end))
plot(TREFHT_30N_DJF(index,12:end)-TREFHT_SSP245(index,12:end))

emu_CIDER = CIDER_response_from_1_forcing(param_T_all(2,:),emu_AOD);
emu_CIDER2 = CIDER_response_from_1_forcing(param_T_all(2,:),CIDER_AOD_from_injection(param_AOD_all(2,:),all_pulse_injections_OG));
plot(emu_CIDER,"LineWidth",4)
% plot(emu_CIDER2,"LineWidth",4)
plot((MAM_dif(1:49)+JJA_dif(1:49)+SON_dif(1:49)+DJF_dif)/4, "LineWidth",4)


ylabel("Temperature Difference (K)")
xlabel("Months after initiation")
ylim([-.8 .5])
title("TREFHT, 6Tg@30N - SSP245")

%% Plot rezeroed

index = 1;
tseries = annualToMonthly(0:4);
figure
tiley = tiledlayout(2,2);
nexttile(1)
set(gca,'DefaultLineLineWidth',2)
hold on
box on
grid on
MAM_dif = AODVIS_15N_MAM(index,3:end)-AODVIS_SSP245(index,3:end);
JJA_dif = AODVIS_15N_JJA(index,6:end)-AODVIS_SSP245(index,6:end);
SON_dif = AODVIS_15N_SON(index,9:end)-AODVIS_SSP245(index,9:end);
DJF_dif = AODVIS_15N_DJF(index,12:end)-AODVIS_SSP245(index,12:end);
plot(AODVIS_15N_MAM(index,3:end)-AODVIS_SSP245(index,3:end))
plot(AODVIS_15N_JJA(index,6:end)-AODVIS_SSP245(index,6:end))
plot(AODVIS_15N_SON(index,9:end)-AODVIS_SSP245(index,9:end))
plot(AODVIS_15N_DJF(index,12:end)-AODVIS_SSP245(index,12:end))
plot(emu_AOD2,"LineWidth",4)
plot((MAM_dif(1:49)+JJA_dif(1:49)+SON_dif(1:49)+DJF_dif)/4, "LineWidth",4)

title("AODVIS, 6Tg@15N - SSP245")
legend("MAM","JJA","SON","DJF","CIDER*","Mean")
ylabel("AOD Difference")
xlabel("Months after initiation")
nexttile(3)
set(gca,'DefaultLineLineWidth',2)
hold on
box on
grid on
MAM_dif = TREFHT_15N_MAM(index,3:end)-TREFHT_SSP245(index,3:end);
JJA_dif = TREFHT_15N_JJA(index,6:end)-TREFHT_SSP245(index,6:end);
SON_dif = TREFHT_15N_SON(index,9:end)-TREFHT_SSP245(index,9:end);
DJF_dif = TREFHT_15N_DJF(index,12:end)-TREFHT_SSP245(index,12:end);
plot(TREFHT_15N_MAM(index,3:end)-TREFHT_SSP245(index,3:end))
plot(TREFHT_15N_JJA(index,6:end)-TREFHT_SSP245(index,6:end))
plot(TREFHT_15N_SON(index,9:end)-TREFHT_SSP245(index,9:end))
plot(TREFHT_15N_DJF(index,12:end)-TREFHT_SSP245(index,12:end))
emu_CIDER_15 = CIDER_response_from_1_forcing(param_T_all(3,:),emu_AOD2);
plot(emu_CIDER_15,"LineWidth",4)
plot((MAM_dif(1:49)+JJA_dif(1:49)+SON_dif(1:49)+DJF_dif)/4, "LineWidth",4)

ylabel("Temperature Difference (K)")
ylim([-.8 .5])
xlabel("Months after initiation")
title("TREFHT, 6Tg@15N - SSP245")
% legend("MAM","JJA","SON","DJF")

nexttile(2)
set(gca,'DefaultLineLineWidth',2)
hold on
box on
grid on
MAM_dif = AODVIS_30N_MAM(index,3:end)-AODVIS_SSP245(index,3:end);
JJA_dif = AODVIS_30N_JJA(index,6:end)-AODVIS_SSP245(index,6:end);
SON_dif = AODVIS_30N_SON(index,9:end)-AODVIS_SSP245(index,9:end);
DJF_dif = AODVIS_30N_DJF(index,12:end)-AODVIS_SSP245(index,12:end);
plot(AODVIS_30N_MAM_12Tg(index,3:end)-AODVIS_SSP245(index,3:end))
plot(AODVIS_30N_JJA_12Tg(index,6:end)-AODVIS_SSP245(index,6:end))
plot(AODVIS_30N_SON_12Tg(index,9:end)-AODVIS_SSP245(index,9:end))
plot(AODVIS_30N_DJF_12Tg(index,12:end)-AODVIS_SSP245(index,12:end))
plot(emu_AOD,"LineWidth",4)
plot((MAM_dif(1:49)+JJA_dif(1:49)+SON_dif(1:49)+DJF_dif)/4, "LineWidth",4)
% plot(emu_AOD_taller,"LineWidth",4)
% % plot(all_pulse_injections/2/10)
% % plot(emu_AOD_CESM,"LineWidth",4,"LineStyle",":")
% plot(emu_AOD_step,"LineWidth",4)
% plot(emu_AOD_CESM_step,"LineWidth",4,"LineStyle",":")
ylabel("AOD Difference 12Tg")
xlabel("Months after initiation")
title("AODVIS, 6Tg@30N - SSP245")
% legend("MAM","JJA","SON","DJF")

nexttile(4)
set(gca,'DefaultLineLineWidth',2)
hold on
box on
grid on
MAM_dif = TREFHT_30N_MAM(index,3:end)-TREFHT_SSP245(index,3:end);
JJA_dif = TREFHT_30N_JJA(index,6:end)-TREFHT_SSP245(index,6:end);
SON_dif = TREFHT_30N_SON(index,9:end)-TREFHT_SSP245(index,9:end);
DJF_dif = TREFHT_30N_DJF(index,12:end)-TREFHT_SSP245(index,12:end);
plot(TREFHT_30N_MAM_12Tg(index,3:end)-TREFHT_SSP245(index,3:end))
plot(TREFHT_30N_JJA_12Tg(index,6:end)-TREFHT_SSP245(index,6:end))
plot(TREFHT_30N_SON_12Tg(index,9:end)-TREFHT_SSP245(index,9:end))
plot(TREFHT_30N_DJF_12Tg(index,12:end)-TREFHT_SSP245(index,12:end))

emu_CIDER = CIDER_response_from_1_forcing(param_T_all(2,:),emu_AOD);
emu_CIDER2 = CIDER_response_from_1_forcing(param_T_all(2,:),CIDER_AOD_from_injection(param_AOD_all(2,:),all_pulse_injections_OG));
plot(emu_CIDER,"LineWidth",4)
% plot(emu_CIDER2,"LineWidth",4)
plot((MAM_dif(1:49)+JJA_dif(1:49)+SON_dif(1:49)+DJF_dif)/4, "LineWidth",4)


ylabel("Temperature Difference (K)")
xlabel("Months after initiation")
ylim([-.8 .5])
title("TREFHT, 6Tg@30N - SSP245")


%%
figure
tiley2 = tiledlayout(2,1);

nexttile
AOD_VIS_30_all = smooth(double(AODVIS_30N_MAM(index,3:end)-AODVIS_SSP245(index,3:end))');

AOD_VIS_30_all = smooth(double(AODVIS_30N_MAM(index,:)-AODVIS_SSP245(index,:))');

% plot(AODVIS_30N_JJA(index,6:end)-AODVIS_SSP245(index,6:end))
% plot(AODVIS_30N_SON(index,9:end)-AODVIS_SSP245(index,9:end))
% plot(AODVIS_30N_DJF(index,12:end)-AODVIS_SSP245(index,12:end))
hold on


plot((1:24)-1,AOD_VIS_30_all(1:24),"LineWidth",4)
plot((2:5)-1,AOD_VIS_30_all(2:5),"LineWidth",4)
plot((1:24)-1,emu_AOD_longer(1:24),"LineWidth",4)
plot((2:5)-1,emu_AOD_longer(2:5),"LineWidth",4)

nexttile
hold on
slope = diff(AOD_VIS_30_all);
emu_slope = diff(emu_AOD_longer);
plot((1:23)+.5-1,slope(1:23),"LineWidth",4)
plot((2:4)+.5-1,slope(2:4),"LineWidth",4)
plot((1:23)+.5-1,diff(emu_AOD_longer(1:24)),"LineWidth",4)
plot((2:4)+.5-1,emu_slope(2:4),"LineWidth",4)


%% Figure
figure
hold on
params1 = paramsE3SM;
coeff = 33.0907;
params1(1)=paramsE3SM(1);
params2 = paramsE3SM;

params2(1) = paramsE3SM(1)*28*.75;
params2(2) = params2(2)*3;
% params2(2) = params2(2);
params2(3) = 0;
aa = [params1;params2];
aa(:,2).^-1
% dualconvolved_bad = CIDER_AOD_from_injection(paramsE3SM,CIDER_AOD_from_injection(paramsE3SM,all_pulse_injections_longer));
dualconvolved = CIDER_AOD_from_injection(params1,CIDER_AOD_from_injection(params2,all_pulse_injections_longer));
% plot(dualconvolved)
% plot(dualconvolved*off_by)
plot((AOD_VIS_30_all))
plot(emu_AOD_longer)
off_by = max(AOD_VIS_30_all)/max(emu_AOD_longer)


% line([0 2],[max(diff(AOD_VIS_30_all)) max(diff(AOD_VIS_30_all))])
%%
% figure
% hold on
% 
% plot(all_pulse_injections_longer)
% plot(SO4_in)
best_param_1 = 6/trapz(all_pulse_injections_longer)*SO2_params(1)
% plot(offby2*CIDER_AOD_from_injection(params2,all_pulse_injections_longer))