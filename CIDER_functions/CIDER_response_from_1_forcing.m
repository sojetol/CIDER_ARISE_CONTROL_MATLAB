function emulated_response = CIDER_response_from_1_forcing(params,forcing)
impulse_p_SAI = [];
impulse_p_SAI.tau = params(1);
impulse_p_SAI.mu = params(2);
emulated_response = zeros(length(forcing),1);
impulse_shape = zeros(length(forcing),1);

for i = 1:length(emulated_response)
impulse_shape(i) = impulse_semiInfDiff(i,impulse_p_SAI);
end

% tic
em_test = conv(impulse_shape,forcing);
% toc

% tic
% % Convolve the impulse response
% for k = 1:length(emulated_response)
%     for j = 1:k
%         emulated_response(k) = emulated_response(k) + impulse_semiInfDiff(j,impulse_p_SAI)*forcing(k-j+1);
%     end
% end
% toc

new_em = em_test(1:length(emulated_response));


emulated_response = new_em;

end