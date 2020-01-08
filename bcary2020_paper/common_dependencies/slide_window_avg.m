function [output] = slide_window_avg(in_vec, wind_avg)


%non overlapping sliding window average
%cuts off last elements in vector if vector is not divisible by wind_avg


output = nanmean(reshape(in_vec(1:wind_avg * floor(numel(in_vec) / wind_avg)), wind_avg, []), 1);



end
