function [data_filted] = filt_data(data, b, a)
data_filted = data;
data_filted(:,1) = filtfilt(b, a, data(:,1));
data_filted(:,2) = filtfilt(b, a, data(:,2));
data_filted(:,3) = filtfilt(b, a, data(:,3));
data_filted(:,4) = filtfilt(b, a, data(:,4));
end

