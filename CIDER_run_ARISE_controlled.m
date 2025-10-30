function [T_array,q_array,log_array] = CIDER_run_ARISE_controlled(CO2_forcing,total_months,init_month,controller_params,CIDER_params,preset_T_noise)
% CIDER_run_ARISE_controlled Run CIDER with the ARISE controller in a
% closed loop.
% CO2_forcing: the background GHG forcing
% total_months: length of simulation, in months
% init_month: what month you start controlling
% controller_params: aspects of the controller (e.g. gains) 
% CIDER_params: govern the dynamics of the aerosols and climate
% preset_T_noise: T0-T1-T2 variability that you want to feed in

% Begin loading in the controller parameters
months_to_average=controller_params.months_to_average; % How many months to average over to get the "current" value (usually 12)
sampling_period = controller_params.sampling_period; % How often to run the controller (1 is monthly, 12 is annual)
refvals = controller_params.refvals; % T0-T1-T2 Target

% Initialize controller and log arrays
q_array = zeros(4,total_months);
log_array = [];

% Load in CIDER parameters
param_AOD_all = CIDER_params.param_AOD_all;
param_T_all = CIDER_params.param_T_all;
pattern_T_all = CIDER_params.pattern_T_all;
pattern_T_base = CIDER_params.pattern_T_base;

% Start up injection -- not sure how to really do this without having a
% little bit of emulation beforehand. More or less a repeat of code further
% down.
if init_month == 1
    time_length = months_to_average;
    zero_array = zeros(time_length,1);
    inj_array = [zero_array q_array(4,1:time_length)'/12 q_array(3,1:time_length)'/12  zero_array q_array(2,1:time_length)'/12 q_array(1,1:time_length)'/12,12  zero_array];
    % for i = 1:4
    %     inj_array(:,i) = CIDER_AOD_from_injection([0.2689 0.2405  0],inj_array(:,i));
    % end
    
    T_array = CIDER_pattern_from_all_injections_and_CO2([ inj_array CO2_forcing(1:time_length)],param_AOD_all,param_T_all,pattern_T_all)+pattern_T_base;
    T_into_controller = mean(T_array(:,:,end-months_to_average+1:end),3);
    T_into_controller = T_array;
    noise_into_controller = mean(preset_T_noise(time_length-months_to_average+1:time_length,:));
    noise_into_controller = preset_T_noise(1:time_length,:);
    [q_out,log_array] = ARISE_controller_for_CIDER(T_into_controller,log_array,controller_params,noise_into_controller);
    q_array(:,init_month:init_month+sampling_period-1)= ones(4,sampling_period).*q_out;
else
    time_length = init_month-1;
    zero_array = zeros(time_length,1);

    inj_array = [zero_array q_array(4,1:(init_month-1))'/12 q_array(3,1:(init_month-1))'/12  zero_array q_array(2,1:(init_month-1))'/12 q_array(1,1:(init_month-1))'/12  zero_array];
    % for i = 1:4
    %     inj_array(:,i) = CIDER_AOD_from_injection([0.2689 0.2405  0],inj_array(:,i));
    % end
    T_array = CIDER_pattern_from_all_injections_and_CO2([ inj_array CO2_forcing(1:time_length)],param_AOD_all,param_T_all,pattern_T_all)+pattern_T_base;
    
    
    if mod(time_length,sampling_period) == 0
        T_into_controller = T_array;
        noise_into_controller = preset_T_noise(1:time_length,:);
        [q_out,log_array] = ARISE_controller_for_CIDER(T_into_controller,log_array,controller_params,noise_into_controller);
        q_array(:,init_month:init_month+sampling_period-1)= ones(4,sampling_period).*q_out;
    end
end


% Emulate through the entire timeseries
for i_month = init_month:total_months

    % Emulation only nessisary when controller is sampled, every
    % sampling_period months
    if mod(i_month,sampling_period) == 0
        time_length = i_month; % length of simulation - this bit could have just been i_month
        zero_array = zeros(time_length,1); % zeros for all places where injection isn't happening
    
        % Create the injection array - SO2 per month
        inj_array = [zero_array q_array(4,1:i_month)'/12 q_array(3,1:i_month)'/12  zero_array q_array(2,1:i_month)'/12 q_array(1,1:i_month)'/12  zero_array];
        
        % This vvv might be something to implement when we go seasonal
        % for i = 1:4
        %     inj_array(:,i) = CIDER_AOD_from_injection([0.2689 0.2405  0],inj_array(:,i));
        % end

        % Emulate using the injection (T_array will have lat, lon, and
        % time)
        T_array = CIDER_pattern_from_all_injections_and_CO2([inj_array CO2_forcing(1:time_length)],param_AOD_all,param_T_all,pattern_T_all)+pattern_T_base;
        
        % Feed T into the controller and get next injection command
        if i_month~=total_months && mod(time_length,sampling_period) == 0
            T_into_controller=T_array; % Unnessisary
            noise_into_controller = preset_T_noise(1:time_length,:); % Feed in only variability that's "happened" so far
            [q_out,log_array] = ARISE_controller_for_CIDER(T_into_controller,log_array,controller_params,noise_into_controller); % Controller
            q_array(:,i_month+1:i_month+sampling_period)= ones(4,sampling_period).*q_out; % Zero-order-hold on the controller's command
        end
    end
end


end