function output = ann(input)
size_input = size(input);
if size_input(1)==1
output = averageEvery2d(12,1,input')';

else

output = averageEvery2d(12,1,input);
end

end