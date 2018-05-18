function [] = plot_excitation_traj(tr)
%close all;


T = 2*pi/tr.w_f;
num = 100;
interval = T/num;
vt = 0:interval:T;

q1 = [];
q2 = [];
q3 = [];
dq1 = [];
dq2 = [];
dq3 = [];
ddq1 = [];
ddq2 = [];
ddq3 = [];

for i = 1:length(vt)
    [q1(i), dq1(i), ddq1(i)] = fourier_func(tr.q01, tr.a1, tr.b1, tr.w_f, vt(i));
end
for i = 1:length(vt)
    [q2(i), dq2(i), ddq2(i)] = fourier_func(tr.q02, tr.a2, tr.b2, tr.w_f, vt(i));
end
for i = 1:length(vt)
    [q3(i), dq3(i), ddq3(i)] = fourier_func(tr.q03, tr.a3, tr.b3, tr.w_f, vt(i));
end

figure
subplot(3,1,1)
plot(vt, q1, '-r', 'LineWidth', 2)
hold on
plot(vt, q2, '--g', 'LineWidth', 2)
plot(vt, q3, '-.b', 'LineWidth', 2)
xlabel('t (s)')
ylabel('Joint angle (rad)')
legend('q1', 'q2', 'q3')

subplot(3,1,2)
plot(vt, dq1, '-r', 'LineWidth', 2)
hold on
plot(vt, dq2, '--g', 'LineWidth', 2)
plot(vt, dq3, '-.b', 'LineWidth', 2)
xlabel('t (s)')
ylabel('Joint velocity (rad/s)')
legend('dq1', 'dq2', 'dq3')

subplot(3,1,3)
plot(vt, ddq1, '-r', 'LineWidth', 2)
hold on
plot(vt, ddq2, '--g', 'LineWidth', 2)
plot(vt, ddq3, '-.b', 'LineWidth', 2)
xlabel('t (s)')
ylabel('Joint acceleration (rad/s^2)')
legend('ddq1', 'ddq2', 'ddq3')

end

