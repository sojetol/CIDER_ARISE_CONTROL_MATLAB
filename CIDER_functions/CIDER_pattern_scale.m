function pattern = CIDER_pattern_scale(mean_response,pattern_to_scale)
    % tic
    reshaped_response = reshape(mean_response,1,1,[]);
    % toc
    % tic
    % b = pattern_to_scale.*a;
    % toc
    % tic
    pattern = pagemtimes(pattern_to_scale,reshaped_response);
    % toc

    % pattern = pattern_to_scale.*reshape(mean_response,1,1,[]);
end