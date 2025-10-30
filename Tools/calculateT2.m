function T2 = calculateT2(T_matrix)

size_of_T_matrix = size(T_matrix);
if length(size_of_T_matrix) == 2
    size_of_T_matrix = [size_of_T_matrix 2];
    T_matrix = cat(3,T_matrix,T_matrix);
end

T_zonal = squeeze(mean(T_matrix,1));

T2 = zeros(size_of_T_matrix(3),1);
% load get_lat_and_lon.mat
lat = linspace(-90,90,size_of_T_matrix(2))';


A = 2*pi*earthRadius^2;

fun = @(j) T_zonal(j,i)*cosd(lat(j))*sind(lat(j));
j = 1:size_of_T_matrix(2)';
for i = 1:size_of_T_matrix(3)
   T2(i) = sum(T_zonal(j,i).*cosd(lat(j)).*(3*sind(lat(j)).^2-1)/2)/sum(cosd(lat(j)));
end

end