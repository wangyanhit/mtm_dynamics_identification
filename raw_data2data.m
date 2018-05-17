function [data] = raw_data2data(data_raw, f)
d_t = 1/f;
dq_f = data_raw(3:end, 2);
dq_b = data_raw(1:end-2, 2);

data = zeros(size(data_raw, 1)-2, size(data_raw, 2)+1);
data(:,1:2) = data_raw(2:end-1,1:2);
data(:,3) = (dq_f -dq_b)/(2*d_t);
data(:,4) = data_raw(2:end-1,3);
end

