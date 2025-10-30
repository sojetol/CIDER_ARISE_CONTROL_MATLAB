%% Initialize
add_all_paths_Sandia % Add paths for directories
controller_params = []; % Empty, will fill later

%% Change things here
plot_title = "CIDER (E3SM-Trained): Various Controllers, Forced Response";

% Fill the parameters (these are the defaults)
controller_params.refvals =  [288.5, 0.8767, -5.89]; % [T0 T1 T2] Targets
controller_params.kivals = [0.0183, 0.0753, 0.3120]; % Integral gains
controller_params.kpvals = [0.0183, 0.0753, 0.3120]; % Proportional Gains
controller_params.months_to_average = 12; % This is your rolling filter. Keep at 12 months unless you're testing something
controller_params.sampling_period = 12; % How often to run the controller and recalculate the injection rate. 1 is monthly, 12 is annual, etc.
controller_params.sens = 4.1; % Something with the feed-forward - high is less? 
controller_params.ff_rate = 0.0273; % Scales feedforward up over time by this amount
controller_params.x_ramp = 5; % Amount of time to scale the feedback in over
controller_params.M = [0, 30, 30, 0;  0, 0, 45, 20; 20, 45, 0, 0; 40, 0, 0, 40];

total_months = 50*12; % Total length of simulation
starting_month = 1+1*12; % When to start controlling (I haven't tried changing it, probably don't touch this)
ensemble_size =1; % Number of ensemble members (each will have different variability)
noise_on = 1; % Turn on natural variability [0/1]
plot_noise = 1; % Plot the natural variability if it's on [0/1]
plot_annual = 1; % Plot the annual mean values [0/1]

controller_params_new_M = controller_params;
controller_params_new_M.M = [0 35.5 35.7 0; 0 10.8 45.5 18.8; 17.3 46.9 10.5 0; 39.1 6.6 4.8 39.7];

% Once all this is set up ^^^ you can easily create new controllers by
% copying and then modifying parameters e.g.:
controller_params_harder = controller_params;
controller_params_harder.kivals = 2*[0.0183, 0.0753, 0.3120];
controller_params_harder.kpvals = 2*[0.0183, 0.0753, 0.3120];

controller_params_even_harder = controller_params;
controller_params_even_harder.kivals = 10*[0.0183, 0.0753, 0.3120];
controller_params_even_harder.kpvals = 10*[0.0183, 0.0753, 0.3120];

controller_params_even_harder_new_M = controller_params_even_harder;
controller_params_even_harder_new_M.M = [0 35.5 35.7 0; 0 10.8 45.5 18.8; 17.3 46.9 10.5 0; 39.1 6.6 4.8 39.7];

controller_params_softer = controller_params;
controller_params_softer.kivals = 0.5*[0.0183, 0.0753, 0.3120];
controller_params_softer.kpvals = 0.5*[0.0183, 0.0753, 0.3120];

controller_params_even_softer = controller_params;
controller_params_even_softer.kivals = 0.1*[0.0183, 0.0753, 0.3120];
controller_params_even_softer.kpvals = 0.1*[0.0183, 0.0753, 0.3120];



uncontrolled = controller_params;
uncontrolled.kivals = [0 0 0];
uncontrolled.kpvals = [0 0 0];
uncontrolled.sens = inf;

% Then you put all your controllers into an array, with names and
% colorschemes for plotting
% all_controllers = [uncontrolled controller_params   controller_params_harder  controller_params_even_harder];
% all_controllers = [uncontrolled controller_params   controller_params_new_M  controller_params_even_harder controller_params_even_harder_new_M];
all_controllers = [uncontrolled controller_params   controller_params_harder  controller_params_softer controller_params_even_softer];
controller_names = [ "No SAI" "ARISE (Normal)" "ARISE (Double Gains)" "ARISE (10x Gains)"];
% controller_names = [ "No SAI" "ARISE (CESM2 M)" "ARISE (E3SM M)" "ARISE (CESM2 M) (10x Gains)" "ARISE (E3SM M) (10x Gains)"];
controller_names = [ "No SAI" "ARISE (Normal)" "ARISE (Double Gains)" "ARISE (half Gains)" "ARISE (10th Gains)"];
colorschemes = ["Reds" "Blues"  "Greens" "Purples" "Oranges"];

%% CIDER Parameters
% Load the AOD and Temperature parameters to make CIDER mimic E3SM
load E3SM_params.mat 
CIDER_params = [];
CIDER_params.param_AOD_all = param_AOD_all;
CIDER_params.param_T_all = param_T_all;
CIDER_params.pattern_T_all = pattern_T_all;
CIDER_params.pattern_T_base = pattern_T_base;
%% CO2
% Set up the background CO2 emissions - Doesn't really need to be changed.
load CO2_concentrations.mat
CO2_ppm_growth = 3;
CO2_ref = CO2_SSP245(20);
CO2_forcing_SSP245 = 5.35*log((CO2_SSP245(21:end)+0)/CO2_ref);
CO2_forcing_SSP245_month = repeatElements(CO2_forcing_SSP245,12);
CO2_SSP245_smooth = linspace(CO2_SSP245(21)-CO2_ppm_growth/2,CO2_SSP245(end)+CO2_ppm_growth/2,length(CO2_forcing_SSP245_month))';
CO2_forcing_SSP245_month_smooth =  5.35*log((CO2_SSP245_smooth)/CO2_ref);

CO2_forcing_SSP370 =  5.35*log((CO2_SSP370(21:end)+0)/CO2_ref);
CO2_forcing_SSP370_month = repeatElements(CO2_forcing_SSP370,12);
%% Make legend
legend_to_plot = ["Target"];
for j = 1:length(all_controllers)
    for i = 1:ensemble_size
        if i == ensemble_size
            legend_to_plot = [legend_to_plot controller_names(j)];
        else
            legend_to_plot = [legend_to_plot ""];
        end
    end
end

%% Run CIDER
load("Tools\PI_T0T1T2.mat") % Load Pre-industrial Variability

% Pre-allocate arrays
ens_T0s = zeros(total_months,ensemble_size,length(all_controllers));
ens_T1s = zeros(total_months,ensemble_size,length(all_controllers));
ens_T2s = zeros(total_months,ensemble_size,length(all_controllers));
tic


for j = 1:length(all_controllers) % For Every Controller

for i = 1:ensemble_size % For all the ensemble members

    % Generate a series of internal variability
    offset = round(500/ensemble_size)*(i-1)*12;
    years_range = 1+offset:total_months+offset;
    preset_T_noise = noise_on*[1*(PI_T0(years_range)-mean(PI_T0(years_range))),...
                        1*(PI_T1(years_range)-mean(PI_T1(years_range))),...
                        1*(PI_T2(years_range)-mean(PI_T2(years_range)))];

    % Run CIDER with controller
    [Temperature,Controls,Log]=CIDER_run_ARISE_controlled(CO2_forcing_SSP245_month_smooth,total_months,starting_month,all_controllers(j),CIDER_params,preset_T_noise);
    
    % Store T0, T1, T2 for plotting
    ens_T0s(:,i,j) =  globalMean(Temperature);
    ens_T1s(:,i,j) =  calculateT1(Temperature);
    ens_T2s(:,i,j) =  calculateT2(Temperature);
    disp("Controller "+j+ ", Ensemble member "+i+" calculated...")

end
end
toc
%% Plot the responses of the controllers
figure
tiley = tiledlayout(3,1);
PIT = 287;
ens_T0s = ens_T0s-PIT;

time_array = annualToMonthly(2035:total_months/12+2034);
for j = 1:length(all_controllers)
    this_controller_params = all_controllers(j);
    controller_colors = brewermap(ensemble_size+2,convertStringsToChars(colorschemes(j)));
for i = 1:ensemble_size
    offset = round(500/ensemble_size)*(i-1)*12;
    years_range = 1+offset:total_months+offset;
    preset_T_noise = noise_on*[1*(PI_T0(years_range)-mean(PI_T0(years_range))),...
                        1*(PI_T1(years_range)-mean(PI_T1(years_range))),...
                        1*(PI_T2(years_range)-mean(PI_T2(years_range)))];
    nexttile(tiley,1)
    set(gca,'DefaultLineLineWidth',4)
    title(plot_title)
    box on 
    grid on
    hold on
    if i == 1 && j == 1
    line([time_array(starting_month) time_array(end)], [this_controller_params.refvals(1)-PIT this_controller_params.refvals(1)-PIT],"LineStyle","--","Color","k","LineWidth",2)
    end
    if plot_annual == 1
        plot(ann(time_array),ann(ens_T0s(:,i,j)+plot_noise*preset_T_noise(:,1)),"Color",controller_colors(i+2,:))
    else
        plot(time_array,ens_T0s(:,i,j)+plot_noise*preset_T_noise(:,1),"Color",controller_colors(i+2,:))
    end
    ylabel("Global Mean T, K above PI")
    xlim([time_array(1) time_array(end)])
    if i == ensemble_size && j == length(all_controllers)
        legend(legend_to_plot,"Location","nw")
    end
    nexttile(tiley,2)
    % set(gca,"DefaultAxesFontSize",16)
    set(gca,'DefaultLineLineWidth',4)
    fontsize(12,"points")
    box on 
    grid on
    hold on
    line([time_array(starting_month) time_array(end)], [this_controller_params.refvals(2) this_controller_params.refvals(2)],"LineStyle","--","Color","k","LineWidth",2)
    if plot_annual == 1
        plot(ann(time_array),ann(ens_T1s(:,i,j)+plot_noise*preset_T_noise(:,2)),"Color",controller_colors(i+2,:))
    else
        plot(time_array,ens_T1s(:,i,j)+plot_noise*preset_T_noise(:,2),"Color",controller_colors(i+2,:))
    end
    ylabel("N-S T Gradient, K")
    xlim([time_array(1) time_array(end)])
    ylim([.7 1.2])
    nexttile(tiley,3)
    set(gca,'DefaultLineLineWidth',4)
    box on 
    grid on
    hold on
    line([time_array(starting_month) time_array(end)], [this_controller_params.refvals(3) this_controller_params.refvals(3)],"LineStyle","--","Color","k","LineWidth",2)
    if plot_annual == 1
        plot(ann(time_array),ann(ens_T2s(:,i,j)+plot_noise*preset_T_noise(:,3)),"Color",controller_colors(i+2,:))
    else
        plot(time_array,ens_T2s(:,i,j)+plot_noise*preset_T_noise(:,3),"Color",controller_colors(i+2,:))
    end
    ylabel("Equator-Pole T Gradient, K")
    xlabel("Year")
    xlim([time_array(1) time_array(end)])
    ylim([-6 -5.6])
end
end