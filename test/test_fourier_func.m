close all;
q0 = 2;
a = [1 0.5 0.25 0.125 0.0625];
b = [1 0.5 0.25 0.125 0.0625];
w_f = 2*pi*0.5;
T = 2*pi/w_f;
vt = 0:0.01:T*2;
q = [];
dq = [];
ddq = [];
for i = 1:length(vt)
    [q(i), dq(i), ddq(i)] = fourier_func(q0, a, b, w_f, vt(i));
end

figure
plot(vt, q, 'r')
hold on
plot(vt, dq , 'g')
plot(vt, ddq, 'b')
xlabel('t(s)')
legend('q', 'velocity', 'acceleration')