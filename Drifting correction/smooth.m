function [output] = smooth(input,smooth_length)
length_data = length(input);
output = zeros(size(input));
for i = 1:length_data
   ini_idx = max(1,i-smooth_length/2);
   end_idx = min(length_data,i+smooth_length/2);
   output(i) = mean(input(ini_idx:end_idx));
end
end