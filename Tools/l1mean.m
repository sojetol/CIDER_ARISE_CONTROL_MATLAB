% Function to calculate L1 mean
function outfield = l1mean(infield, weights, lats)
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
    lats2 = permute(repmat(lats, [12, b(2), 1]), [3, 2, 1]);
    sinlats = sin(deg2rad(lats2));
    theprod = infield2 .* w3 .* sinlats;
    
    while ndims(theprod) > 1
        theprod = sum(theprod, 1);
        w3 = sum(w3, 1);
    end
    
    outfield = theprod ./ w3;
end