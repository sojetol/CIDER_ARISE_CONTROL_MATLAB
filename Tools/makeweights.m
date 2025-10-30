% Function to create an array of weights for area weighting
function weights = makeweights(lats, lons)
    nlats = length(lats);
    nlons = length(lons);
    latdiff = diff(lats);
    temparray = -90 + latdiff(1)/2 : latdiff(1) : 90;
    latedges = [-90, temparray, 90];
    colatedges = pi/2 - latedges * (pi / 180);
    dphi = (lons(3) - lons(2)) * (pi / 180);
    areaout = abs(diff(cos(colatedges))) * dphi;
    biggest = max(areaout);
    mapcols = areaout / biggest;
    weights = repmat(mapcols, nlons, 1);
end