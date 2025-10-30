% Function to calculate global mean
function outfield = gmean(infield, weights)
    a = size(infield);
    b = size(weights);
    otherinds = -ones(1, length(a));
    
    for k = 1:length(b)
        otherinds(k) = find(a == b(k));
    end
    otherinds(otherinds == -1) = [];
    
    timeind = find(~ismember(a, b));
    
    % Want time to be last
    c = perms(a);
    d = [b, a(timeind)];
    i = find(ismember(c, d, 'rows'));
    e = perms(1:length(a));
    infield2 = permute(infield, e(i, :));
    w3 = reshape(repmat(weights, [1, a(timeind)]), d);
    theprod = infield2 .* w3;
    
    while ndims(theprod) > 1
        theprod = sum(theprod, 1);
        w3 = sum(w3, 1);
    end
    
    outfield = theprod ./ w3;
end
