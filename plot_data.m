function [] = plot_data(q1_data, q2_data, q3_data, q1_data_filted, q2_data_filted, q3_data_filted, f)
figure
d_t = 1/f;
t = ((1:size(q1_data))-1)*d_t;

subplot(4,3,1)
plot(t, q1_data(:,1), t, q1_data_filted(:,1))
xlabel("t (s)");
ylabel("Angle (rad)");
subplot(4,3,4)
plot(t, q1_data(:,2), t, q1_data_filted(:,2))
xlabel("t (s)");
ylabel("Velocity (rad/s)");
subplot(4,3,7)
plot(t, q1_data(:,3), t, q1_data_filted(:,3))
xlabel("t (s)");
ylabel("Acceleration (rad/s^2)");
subplot(4,3,10)
plot(t, q1_data(:,4), t, q1_data_filted(:,4))
xlabel("t (s)");
ylabel("Effort (N*m)");


subplot(4,3,2)
plot(t, q2_data(:,1), t, q2_data_filted(:,1))
xlabel("t (s)");
ylabel("Angle (rad)");
subplot(4,3,5)
plot(t, q2_data(:,2), t, q2_data_filted(:,2))
xlabel("t (s)");
ylabel("Velocity (rad/s)");
subplot(4,3,8)
plot(t, q2_data(:,3), t, q2_data_filted(:,3))
xlabel("t (s)");
ylabel("Acceleration (rad/s^2)");
subplot(4,3,11)
plot(t, q2_data(:,4), t, q2_data_filted(:,4))
xlabel("t (s)");
ylabel("Effort (N*m)");


subplot(4,3,3)
plot(t, q3_data(:,1), t, q3_data_filted(:,1))
xlabel("t (s)");
ylabel("Angle (rad)");
subplot(4,3,6)
plot(t, q3_data(:,2), t, q3_data_filted(:,2))
xlabel("t (s)");
ylabel("Velocity (rad/s)");
subplot(4,3,9)
plot(t, q3_data(:,3), t, q3_data_filted(:,3))
xlabel("t (s)");
ylabel("Acceleration (rad/s^2)");
subplot(4,3,12)
plot(t, q3_data(:,4), t, q3_data_filted(:,4))
xlabel("t (s)");
ylabel("Effort (N*m)");

end

