function all_params = CIDER_train_AOD_params_seasonal(all_pulse_injections,all_pulse_responses,all_feedback_injections,all_feedback_responses,weightings,varargin)
pulse_injections_size = size(all_pulse_injections);
pulse_injections_count = pulse_injections_size(2);
pulse_injections_length = pulse_injections_size(1);

all_params = zeros(pulse_injections_count,3);
if nargin>5
    all_short_step_injections = varargin{1};
    all_simulated_short_step_response = varargin{2};
    for i = 1:pulse_injections_count
        step_injection = all_pulse_injections(:,i);
        %feedback_injection = all_feedback_injections(:,i);
        simulated_pulse_response = all_pulse_responses(:,i);
        %simulated_feedback_response = all_feedback_responses(:,i);
        short_step_injection = all_short_step_injections(:,i);
        simulated_short_step_response = all_simulated_short_step_response(:,i);
        short_step_error = @(params) sum((simulated_short_step_response-CIDER_AOD_from_injection(params,short_step_injection)).^2) / length(simulated_short_step_response);
        pulse_error = @(params) sum((simulated_pulse_response-CIDER_AOD_from_injection(params,step_injection)).^2) / length(simulated_pulse_response);
        %feedback_error = @(params) sum((simulated_feedback_response-CIDER_AOD_from_injection(params,feedback_injection)).^2) / length(simulated_feedback_response);
        emulator_error = @(params) weightings(1)*pulse_error(params)+weightings(1)*short_step_error(params);%+weightings(2)*feedback_error(params);
        x0 = [0.02; 0.01; 0.0002];
        % optimal_params = fmincon(@(x) emulator_error(x),x0,[],[],[],[],[0;0;0],[.1;.2;.1]);
        optimal_params = fmincon(@(x) emulator_error(x),x0,[],[],[],[],[0;0.08;0],[.1;.2;1]);
        all_params(i,:) = optimal_params;
    end
else
for i = 1:pulse_injections_count
    pulse_injection = all_pulse_injections(:,i);
    %feedback_injection = all_feedback_injections(:,i);
    simulated_pulse_response = all_pulse_responses(:,i);
    %simulated_feedback_response = all_feedback_responses(:,i);
    pulse_error = @(params) sum((simulated_pulse_response-CIDER_AOD_from_injection(params,pulse_injection)).^2) / length(simulated_pulse_response);
    feedback_error = @(params) sum((simulated_feedback_response-CIDER_AOD_from_injection(params,feedback_injection)).^2) / length(simulated_feedback_response);
    emulator_error = @(params) weightings(1)*pulse_error(params);%+weightings(2)*feedback_error(params);
    x0 = [0.02; 0.1; 0.0002];
    % optimal_params = fmincon(@(x) emulator_error(x),x0,[],[],[],[],[0;0;0],[.1;.2;.1]);
    optimal_params = fmincon(@(x) emulator_error(x),x0,[],[],[],[],[0;0.08;0],[1;.3;.1]);
    all_params(i,:) = optimal_params;
end
end

end