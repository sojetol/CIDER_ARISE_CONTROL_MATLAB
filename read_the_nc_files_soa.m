add_all_paths_Sandia
load CESM_params.mat
AODVIS_15N_2015_2019 = ncread("Climate_data\soa_15N_AODVIS_201501_201912_readable.nc","AODVIS");
AODVIS_15N_2020_2024 = ncread("Climate_data\soa_15N_AODVIS_202001_202412_readable.nc","AODVIS");
AODVIS_15N_2025_2029 = ncread("Climate_data\soa_15N_AODVIS_202501_202912_readable.nc","AODVIS");
TREFHT_15N_2015_2019 = ncread("Climate_data\soa_15N_TREFHT_201501_201912_readable.nc","TREFHT");
TREFHT_15N_2020_2024 = ncread("Climate_data\soa_15N_TREFHT_202001_202412_readable.nc","TREFHT");
TREFHT_15N_2025_2029 = ncread("Climate_data\soa_15N_TREFHT_202501_202912_readable.nc","TREFHT");
AODVIS_15N_all = [AODVIS_15N_2015_2019,AODVIS_15N_2020_2024,AODVIS_15N_2025_2029];
TREFHT_15N_all = [TREFHT_15N_2015_2019,TREFHT_15N_2020_2024,TREFHT_15N_2025_2029];

AODVIS_15S_2015_2019 = ncread("Climate_data\soa_15S_AODVIS_201501_201912_readable.nc","AODVIS");
AODVIS_15S_2020_2024 = ncread("Climate_data\soa_15S_AODVIS_202001_202412_readable.nc","AODVIS");
AODVIS_15S_2025_2029 = ncread("Climate_data\soa_15S_AODVIS_202501_202912_readable.nc","AODVIS");
TREFHT_15S_2015_2019 = ncread("Climate_data\soa_15S_TREFHT_201501_201912_readable.nc","TREFHT");
TREFHT_15S_2020_2024 = ncread("Climate_data\soa_15S_TREFHT_202001_202412_readable.nc","TREFHT");
TREFHT_15S_2025_2029 = ncread("Climate_data\soa_15S_TREFHT_202501_202912_readable.nc","TREFHT");
AODVIS_15S_all = [AODVIS_15S_2015_2019,AODVIS_15S_2020_2024,AODVIS_15S_2025_2029];
TREFHT_15S_all = [TREFHT_15S_2015_2019,TREFHT_15S_2020_2024,TREFHT_15S_2025_2029];

AODVIS_30N_2015_2019 = ncread("Climate_data\soa_30N_AODVIS_201501_201912_readable.nc","AODVIS");
AODVIS_30N_2020_2024 = ncread("Climate_data\soa_30N_AODVIS_202001_202412_readable.nc","AODVIS");
AODVIS_30N_2025_2029 = ncread("Climate_data\soa_30N_AODVIS_202501_202912_readable.nc","AODVIS");
TREFHT_30N_2015_2019 = ncread("Climate_data\soa_30N_TREFHT_201501_201912_readable.nc","TREFHT");
TREFHT_30N_2020_2024 = ncread("Climate_data\soa_30N_TREFHT_202001_202412_readable.nc","TREFHT");
TREFHT_30N_2025_2029 = ncread("Climate_data\soa_30N_TREFHT_202501_202912_readable.nc","TREFHT");
AODVIS_30N_all = [AODVIS_30N_2015_2019,AODVIS_30N_2020_2024,AODVIS_30N_2025_2029];
TREFHT_30N_all = [TREFHT_30N_2015_2019,TREFHT_30N_2020_2024,TREFHT_30N_2025_2029];

AODVIS_30S_2015_2019 = ncread("Climate_data\soa_30S_AODVIS_201501_201912_readable.nc","AODVIS");
AODVIS_30S_2020_2024 = ncread("Climate_data\soa_30S_AODVIS_202001_202412_readable.nc","AODVIS");
AODVIS_30S_2025_2029 = ncread("Climate_data\soa_30S_AODVIS_202501_202912_readable.nc","AODVIS");
TREFHT_30S_2015_2019 = ncread("Climate_data\soa_30S_TREFHT_201501_201912_readable.nc","TREFHT");
TREFHT_30S_2020_2024 = ncread("Climate_data\soa_30S_TREFHT_202001_202412_readable.nc","TREFHT");
TREFHT_30S_2025_2029 = ncread("Climate_data\soa_30S_TREFHT_202501_202912_readable.nc","TREFHT");
AODVIS_30S_all = [AODVIS_30S_2015_2019,AODVIS_30S_2020_2024,AODVIS_30S_2025_2029];
TREFHT_30S_all = [TREFHT_30S_2015_2019,TREFHT_30S_2020_2024,TREFHT_30S_2025_2029];

AODVIS_SSP245 = [ncread("Climate_data\AODVIS_SSP245_201501_201912_readable.nc","AODVIS")...
                 ncread("Climate_data\AODVIS_SSP245_202001_202412_readable.nc","AODVIS")...
                 ncread("Climate_data\AODVIS_SSP245_202501_202912_readable.nc","AODVIS")];

TREFHT_SSP245 = [ncread("Climate_data\TREFHT_SSP245_201501_201912_readable.nc","TREFHT")...
                 ncread("Climate_data\TREFHT_SSP245_202001_202412_readable.nc","TREFHT")...
                 ncread("Climate_data\TREFHT_SSP245_202501_202912_readable.nc","TREFHT")];


index = 1;

TREFHT_SSP245_pattern = cat(3,...
                            ncread("Climate_data\TREFHT_SSP245_202001_202412_readable_180x360_aave.nc","TREFHT"),...
                            ncread("Climate_data\TREFHT_SSP245_202501_202912_readable_180x360_aave.nc","TREFHT"));

TREFHT_15N_pattern_ts = cat(3,ncread("Climate_data\soa_15N_180x360_aave_TREFHT_202001_202412_readable.nc","TREFHT"),...
                           ncread("Climate_data\soa_15N_180x360_aave_TREFHT_202501_202912_readable.nc","TREFHT"))-TREFHT_SSP245_pattern;

TREFHT_15S_pattern_ts = cat(3,ncread("Climate_data\soa_15S_180x360_aave_TREFHT_202001_202412_readable.nc","TREFHT"),...
                           ncread("Climate_data\soa_15S_180x360_aave_TREFHT_202501_202912_readable.nc","TREFHT"))-TREFHT_SSP245_pattern;

TREFHT_30N_pattern_ts = cat(3,ncread("Climate_data\soa_30N_180x360_aave_TREFHT_202001_202412_readable.nc","TREFHT"),...
                           ncread("Climate_data\soa_30N_180x360_aave_TREFHT_202501_202912_readable.nc","TREFHT"))-TREFHT_SSP245_pattern;
TREFHT_30S_pattern_ts = cat(3,ncread("Climate_data\soa_30S_180x360_aave_TREFHT_202001_202412_readable.nc","TREFHT"),...
                           ncread("Climate_data\soa_30S_180x360_aave_TREFHT_202501_202912_readable.nc","TREFHT"))-TREFHT_SSP245_pattern;

TREFHT_SSP245_pattern_base = mean(ncread("Climate_data\TREFHT_SSP245_201501_201912_readable_180x360_aave.nc","TREFHT"),3);

pattern_T_15N = mean(TREFHT_15N_pattern_ts(:,:,25:end),3)/mean(globalMean(TREFHT_15N_pattern_ts(:,:,25:end)));
pattern_T_15S = mean(TREFHT_15S_pattern_ts(:,:,25:end),3)/mean(globalMean(TREFHT_15S_pattern_ts(:,:,25:end)));
pattern_T_30N = mean(TREFHT_30N_pattern_ts(:,:,25:end),3)/mean(globalMean(TREFHT_30N_pattern_ts(:,:,25:end)));
pattern_T_30S = mean(TREFHT_30S_pattern_ts(:,:,25:end),3)/mean(globalMean(TREFHT_30S_pattern_ts(:,:,25:end)));
pattern_T_ssp245 = mean(TREFHT_SSP245_pattern(:,:,:)-TREFHT_SSP245_pattern_base,3)/mean(globalMean(TREFHT_SSP245_pattern(:,:,:)-TREFHT_SSP245_pattern_base));

%%
pattern_T_all(:,:,2) = regrideE3SMtoCESM(pattern_T_30N);
pattern_T_all(:,:,3) = regrideE3SMtoCESM(pattern_T_15N);
pattern_T_all(:,:,5) = regrideE3SMtoCESM(pattern_T_15S);
pattern_T_all(:,:,6) = regrideE3SMtoCESM(pattern_T_30S);
pattern_T_all(:,:,8) = regrideE3SMtoCESM(pattern_T_ssp245);
%%
AODVIS_30N_diff = AODVIS_30N_all(index,61:end) - AODVIS_SSP245(index,61:end);
AODVIS_15N_diff = AODVIS_15N_all(index,61:end) - AODVIS_SSP245(index,61:end);
AODVIS_15S_diff = AODVIS_15S_all(index,61:end) - AODVIS_SSP245(index,61:end);
AODVIS_30S_diff = AODVIS_30S_all(index,61:end) - AODVIS_SSP245(index,61:end);

TREFHT_30N_diff = TREFHT_30N_all(index,61:end) - TREFHT_SSP245(index,61:end);
TREFHT_15N_diff = TREFHT_15N_all(index,61:end) - TREFHT_SSP245(index,61:end);
TREFHT_15S_diff = TREFHT_15S_all(index,61:end) - TREFHT_SSP245(index,61:end);
TREFHT_30S_diff = TREFHT_30S_all(index,61:end) - TREFHT_SSP245(index,61:end);

%%
% plot(AODVIS_30N_SON(index,9:end)-AODVIS_SSP245(index,9:end))
% plot(AODVIS_30N_DJF(index,12:end)-AODVIS_SSP245(index,12:end))
% param_T_all(:,1) = 200*param_T_all(:,1)
% param_T_all(:,2) = param_T_all(:,2)/20


all_step_injections = .5*ones(120,4);
all_step_responses_AOD = double([AODVIS_30N_diff' AODVIS_15N_diff' AODVIS_15S_diff' AODVIS_30S_diff']);
all_step_responses_T = double([TREFHT_30N_diff' TREFHT_15N_diff' TREFHT_15S_diff' TREFHT_30S_diff']);

param_AOD_all_E3SM = CIDER_train_AOD_params_seasonal(all_step_injections,all_step_responses_AOD,[],[],[1 0])
param_T_all_E3SM = CIDER_train_climate_params(all_step_responses_AOD,all_step_responses_T)


emu_AOD_30N_CESM = CIDER_AOD_from_injection(param_AOD_all(2,:),.5*ones(120,1));
emu_AOD_15N_CESM = CIDER_AOD_from_injection(param_AOD_all(3,:),.5*ones(120,1));
emu_AOD_15S_CESM = CIDER_AOD_from_injection(param_AOD_all(5,:),.5*ones(120,1));
emu_AOD_30S_CESM = CIDER_AOD_from_injection(param_AOD_all(6,:),.5*ones(120,1));

emu_T_30N_CESM = CIDER_response_from_1_forcing(param_T_all(2,:),emu_AOD_30N_CESM);
emu_T_15N_CESM = CIDER_response_from_1_forcing(param_T_all(3,:),emu_AOD_15N_CESM);
emu_T_15S_CESM = CIDER_response_from_1_forcing(param_T_all(5,:),emu_AOD_15S_CESM);
emu_T_30S_CESM = CIDER_response_from_1_forcing(param_T_all(6,:),emu_AOD_30S_CESM);

param_T_all_old = param_T_all;
param_AOD_all_old = param_AOD_all;
param_T_all(2:3,:) = param_T_all_E3SM(1:2,:);
param_T_all(5:6,:) = param_T_all_E3SM(3:4,:);
param_AOD_all(2:3,:) = param_AOD_all_E3SM(1:2,:);
param_AOD_all(5:6,:) = param_AOD_all_E3SM(3:4,:);

param_AOD_comp = [param_AOD_all param_AOD_all_old]
param_T_comp = [param_T_all param_T_all_old]
%%
save("E3SM_params.mat","param_AOD_all","param_P_all","param_T_all","param_Q_all","param_RH_all","param_U10_all","pattern_P_all","pattern_T_all","pattern_Q_all","pattern_RH_all","pattern_U10_all","pattern_AOD_all","pattern_AOD_base","pattern_T_base","pattern_P_base","pattern_Q_base","pattern_U10_base","pattern_RH_base");
%%
emu_AOD_30N = CIDER_AOD_from_injection(param_AOD_all_E3SM(1,:),.5*ones(120,1));
emu_AOD_15N = CIDER_AOD_from_injection(param_AOD_all_E3SM(2,:),.5*ones(120,1));
emu_AOD_15S = CIDER_AOD_from_injection(param_AOD_all_E3SM(3,:),.5*ones(120,1));
emu_AOD_30S = CIDER_AOD_from_injection(param_AOD_all_E3SM(4,:),.5*ones(120,1));

emu_T_30N = CIDER_response_from_1_forcing(param_T_all_E3SM(1,:),emu_AOD_30N);
emu_T_15N = CIDER_response_from_1_forcing(param_T_all_E3SM(2,:),emu_AOD_15N);
emu_T_15S = CIDER_response_from_1_forcing(param_T_all_E3SM(3,:),emu_AOD_15S);
emu_T_30S = CIDER_response_from_1_forcing(param_T_all_E3SM(4,:),emu_AOD_30S);


tseries = annualToMonthly(0:9);
figure
tiley= tiledlayout(2,1);
nexttile
hold on
box on
grid on
set(gca,'DefaultLineLineWidth',4)
set(gca,'FontSize',12)
plot(tseries,AODVIS_30N_diff,"Color",[230, 53, 65]/255)
plot(tseries,AODVIS_15N_diff,"Color",[0, 168, 31]/255)
plot(tseries,AODVIS_15S_diff,"Color",[66, 135, 245]/255)
plot(tseries,AODVIS_30S_diff,"Color",[247, 165, 49]/255)
plot(tseries,emu_AOD_30N,"Color",[230, 53, 65]/255,"LineStyle",":")
plot(tseries,emu_AOD_15N,"Color",[0, 168, 31]/255,"LineStyle",":")
plot(tseries,emu_AOD_15S,"Color",[66, 135, 245]/255,"LineStyle",":")
plot(tseries,emu_AOD_30S,"Color",[247, 165, 49]/255,"LineStyle",":")
ylabel("Change in Aerosol Optical Depth")
title("Global Mean Responses to Single-Latitude 6 Tg SO_2 Injections")
legend("30N, E3SM","15N, E3SM","15S, E3SM","30S, E3SM","30N, CIDER","15N, CIDER","15S, CIDER", "30S, CIDER","Location","se","NumColumns",2)
nexttile
hold on
box on 
grid on
set(gca,'DefaultLineLineWidth',4)
set(gca,'FontSize',12)
plot(tseries,TREFHT_30N_diff,"Color",[230, 53, 65]/255)
plot(tseries,TREFHT_15N_diff,"Color",[0, 168, 31]/255)
plot(tseries,TREFHT_15S_diff,"Color",[66, 135, 245]/255)
plot(tseries,TREFHT_30S_diff,"Color",[247, 165, 49]/255)
plot(tseries,emu_T_30N,"Color",[230, 53, 65]/255,"LineStyle",":")
plot(tseries,emu_T_15N,"Color",[0, 168, 31]/255,"LineStyle",":")
plot(tseries,emu_T_15S,"Color",[66, 135, 245]/255,"LineStyle",":")
plot(tseries,emu_T_30S,"Color",[247, 165, 49]/255,"LineStyle",":")
ylabel("Near-Surface Cooling, K")
xlabel("Years of Injection")
% plot(AODVIS_15N_diff)
% plot(AODVIS_15S_diff)
% plot(AODVIS_30S_diff)
% plot(tseries,emu_AOD_30S)
% plot(emu_AOD_30S_CESM)
% plot(emu_AOD_15N)
% plot(emu_AOD_15S)
% plot(emu_AOD_30S)
% 
% figure
% hold on
% plot(TREFHT_30S_diff)
% 
% % plot(TREFHT_15N_diff)
% % plot(TREFHT_15S_diff)
% % plot(TREFHT_30S_diff)
% plot(emu_T_30S)
% plot(emu_T_30S_CESM)
% % plot(emu_T_15N)
% % plot(emu_T_15S)
% % plot(emu_T_30S)