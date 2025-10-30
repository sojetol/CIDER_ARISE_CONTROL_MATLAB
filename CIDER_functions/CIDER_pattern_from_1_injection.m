function response_pattern = CIDER_pattern_from_1_injection(injection,param_AOD,param_climate,pattern_to_scale)
% clc
% tic
AOD_emulated = CIDER_AOD_from_injection(param_AOD,injection);
% toc
% tic
emulated_response = CIDER_response_from_1_forcing(param_climate,AOD_emulated);
% toc
% tic
response_pattern = CIDER_pattern_scale(emulated_response,pattern_to_scale);
% toc
if length(emulated_response)>400
    dummy = 1;
end

end