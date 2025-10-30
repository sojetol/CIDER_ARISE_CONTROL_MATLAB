function AOD_emulated = CIDER_AOD_from_injection(param_AOD,injection)
% Load in parameters for AOD
beta = param_AOD(1);
alpha = param_AOD(2);
gamma = param_AOD(3);

AOD_emulated = zeros(length(injection),1);
if mean(injection)==0
    return
end

% impulse_shape = zeros(length(injection),1);
% 
% for i = 1:length(AOD_emulated)
% impulse_shape(i) = impulse_semiInfDiff(i,impulse_p_SAI);
% end






% tic
% for k = 1:length(AOD_emulated)
%     for j = 1:k
%         AOD_emulated(k) = AOD_emulated(k) + impulse_firstOdiff_wGammaInj(beta, alpha,gamma,injection(k-j+1), j)*injection(k-j+1);
% 
% 
%         % if ydata(k-j+1,j)-impulse_firstOdiff_wGammaInj(beta, alpha,gamma,injection(k-j+1), j) ~=0
%         %     aaa = 1;
%         % end
%     end
% end
% toc

% tic

t = (1:length(AOD_emulated));

ydata = beta*exp(-(alpha+gamma*injection)*t);
for k = 1:length(AOD_emulated)
    for j = 1:k
        AOD_emulated(k) = AOD_emulated(k) + ydata(k-j+1,j)*injection(k-j+1);
        
        
        % if ydata(k-j+1,j)-impulse_firstOdiff_wGammaInj(beta, alpha,gamma,injection(k-j+1), j) ~=0
        %     aaa = 1;
        % end
    end
end
% toc
end
