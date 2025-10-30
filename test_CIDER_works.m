add_all_paths_Sandia
load CESM_params.mat
% CESM_AOD_param = param_AOD_all
% load UKESM_params.mat
% UKESM_AOD_param = param_AOD_all
% param_AOD_all = CESM_AOD_param

%%
time_length = 60;
inj_test_0N = [zeros(time_length,3) ones(time_length,1) zeros(time_length,3) zeros(time_length,1)];
myseasonalinj = [zeros(11,1);
                  0;2;2;
                  2;0;0;
                  0;0;0;
                  0;0;0;
                  zeros(time_length-24+1,1)];
inj_test_0N = [zeros(time_length,3) myseasonalinj zeros(time_length,3) zeros(time_length,1)];
inj_test_30N = [zeros(time_length,1) myseasonalinj zeros(time_length,5) zeros(time_length,1)];
inj_test_30S = [zeros(time_length,5) myseasonalinj zeros(time_length,1) zeros(time_length,1)];
T_tested_0N = CIDER_response_from_all_injections_and_CO2(inj_test_0N,param_AOD_all,param_T_all);
T_tested_30N = CIDER_response_from_all_injections_and_CO2(inj_test_30N,param_AOD_all,param_T_all);
T_tested_30S = CIDER_response_from_all_injections_and_CO2(inj_test_30S,param_AOD_all,param_T_all);
time_years = linspace(-1+1/12,time_length/12-1,time_length);

figure
all_plot_on
plot(time_years,T_tested_30N,"LineWidth",4)
plot(time_years,T_tested_0N,"LineWidth",4)
plot(time_years,T_tested_30S,"LineWidth",4)
xlabel("Years after Injection")
ylabel("Cooling (K)")
title("Response to 6Tg, 3-Month Injection in CESM (CIDER's Guess)")
legend("30N","0N","30S")