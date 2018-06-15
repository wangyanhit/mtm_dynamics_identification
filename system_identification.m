%% Define geometric parameters of robot
clc;clear;
close all

dof = 7;
pi = sym('pi');
%%
% Generate joint variables

%%
% Serial link method
syms q1 q2 q3 q4 q5 q6 q7 real;
syms dq1 dq2 dq3 dq4 dq5 dq6 dq7 real;
syms ddq1 ddq2 ddq3 ddq4 ddq5 ddq6 ddq7 real;
q = [q1 q2 q3 q4 q5 q6 q7];
q2p = q2 + q3;
q2pp = -q3;
dq = [dq1 dq2 dq3 dq4 dq5 dq6 dq7];
dq2p = dq2 + dq3;
dq2pp = -dq3;
dq = [ddq1 ddq2 ddq3 ddq4 ddq5 ddq6 ddq7];
ddq2p = ddq2 + ddq3;
ddq2pp = -ddq3;

%%
% paralell link method
% syms q1 q2 q2p q4 q5 q6 q7 t real;
% syms dq1 dq2 dq2p dq4 dq5 dq6 dq7 t real;
% q = [q1 q2 q2p q4 q5 q6 q7];
% q3 = -q2 + q2p;
% q2pp = -q2p + q2;
% dq = [dq1 dq2 dq2p dq4 dq5 dq6 dq7];
% dq3 = -dq2 + dq2p;
% dq2pp = -dq2p + dq2; 

%%
% Define DH parameters
% Order: $a, \alpha, d \theta$
dh_link1 = [0 pi/2 0 q1-pi/2]; % 1
dh_link2 = [0.279 0 0 q2-pi/2]; % 1 to 2
dh_link2p = [0.1 0 0 q2p]; % 1 to 2p
dh_link2pp = [0.279 0 0 q2pp-pi/2]; % 2p to 2pp
dh_link3 = [0.365 -pi/2 0 q3+pi/2]; % 2 to 3
dh_link4 = [0 pi/2 0.151 q4]; % 3 to 4
dh_link5 = [0 -pi/2 0 q5]; % 4 to 5
dh_link6 = [0 pi/2 0 q6-pi/2]; % 5 to 6
dh_link7 = [0 0 0 q7]; % 6 to 7

T01 = dhparam2matrix(dh_link1);
T12 = dhparam2matrix(dh_link2);
T23 = dhparam2matrix(dh_link3);
T34 = dhparam2matrix(dh_link4);
T45 = dhparam2matrix(dh_link5);
T56 = dhparam2matrix(dh_link6);
T67 = dhparam2matrix(dh_link7);
T12p = dhparam2matrix(dh_link2p);
T2p2pp = dhparam2matrix(dh_link2pp);

T02 = T01*T12;
T03 = T02*T23;
T04 = T03*T34;
T05 = T04*T45;
T06 = T05*T56;
T07 = T06*T67;
T02p = T01*T12p;
T2pp = T02p*T2p2pp;

T(:,:,1) = T01;
T(:,:,2) = T02;
T(:,:,3) = T03;
T(:,:,4) = T04;
T(:,:,5) = T05;
T(:,:,6) = T06;
T(:,:,7) = T07;
T(:,:,8) = T02p;
T(:,:,9) = T2pp;

%%
% Visualization to see if the transformations are right
visualize_transform = true;
if visualize_transform
    figure
    hold on
    axis equal
    for idx = 1:9
        %T_num = subs(T(:,:,idx), {q1 q2 q2p q4 q5 q6 q7}, {0 0 0 0 0.2 0.0 0.4});
        T_num = subs(T(:,:,idx), {q1 q2 q3 q4 q5 q6 q7}, {0 -0.2 0.2 0 0.2 0.0 0.4});
        p = T_num*([0 0 0 1]');
        x = T_num*([0.1 0 0 1]');
        y = T_num*([0 0.1 0 1]');
        z = T_num*([0 0 0.1 1]');
        plot3([p(1) x(1)], [p(2) x(2)], [p(3) x(3)], 'r')
        plot3([p(1) y(1)], [p(2) y(2)], [p(3) y(3)], 'g')
        plot3([p(1) z(1)], [p(2) z(2)], [p(3) z(3)], 'b')
        text(p(1), p(2), p(3), num2str(idx));
    end
    legend('x', 'y', 'z')
    xlabel('x (m)')
    ylabel('y (m)')
    zlabel('z (m)')
    hold off
end

%%
% Define inertia parameters
syms Lxx1 Lxy1 Lxz1 Lyy1 Lyz1 Lzz1 lx1 ly1 lz1 m1 real;
syms Fv1 Fc1 Fo1 real;
[delta_L1, r1, I1] = inertia_bary2std(Lxx1, Lxy1, Lxz1, Lyy1, Lyz1, Lzz1, lx1, ly1, lz1, m1);
delta_A1 = [Fv1;];% Fc1; Fo1];
%delta_A1 = [];
syms Lxx2 Lxy2 Lxz2 Lyy2 Lyz2 Lzz2 lx2 ly2 lz2 m2 real;
syms Fv2 Fc2 Fo2 real;
[delta_L2, r2, I2] = inertia_bary2std(Lxx2, Lxy2, Lxz2, Lyy2, Lyz2, Lzz2, lx2, ly2, lz2, m2);
delta_A2 = [Fv2;];% Fc2; Fo2];
%delta_A2 = [];
syms Lxx3 Lxy3 Lxz3 Lyy3 Lyz3 Lzz3 lx3 ly3 lz3 m3 real;
syms Fv3 Fc3 Fo3 real;
[delta_L3, r3, I3] = inertia_bary2std(Lxx3, Lxy3, Lxz3, Lyy3, Lyz3, Lzz3, lx3, ly3, lz3, m3);
delta_A3 = [Fv3;];% Fc3; Fo3];
%delta_A3 = [];
%%
% Linear and rotational velocities of link mass centers 
% Tranformations for mass centers
T01_mc = simplify(T01*trans_mat(r1));
T02_mc = simplify(T02*trans_mat(r2));
T03_mc = simplify(T03*trans_mat(r3));

p01_mc = T01_mc(1:3,4);
p02_mc = T02_mc(1:3,4);
p03_mc = T03_mc(1:3,4);

syms t real;
syms q1t(t) q2t(t) q3t(t)
T01_mc_t = subs(T01_mc, {q1, q2, q3}, {q1t, q2t, q3t});
T02_mc_t = subs(T02_mc, {q1, q2, q3}, {q1t, q2t, q3t});
T03_mc_t = subs(T03_mc, {q1, q2, q3}, {q1t, q2t, q3t});

dT01_mc_t = simplify(diff(T01_mc_t, t));
dT02_mc_t = simplify(diff(T02_mc_t, t));
dT03_mc_t = simplify(diff(T03_mc_t, t));

dT01_mc = subs(dT01_mc_t, {diff(q1t(t), t), diff(q2t(t), t), diff(q3t(t), t)}, {dq1, dq2, dq3});
dT02_mc = subs(dT02_mc_t, {diff(q1t(t), t), diff(q2t(t), t), diff(q3t(t), t)}, {dq1, dq2, dq3});
dT03_mc = subs(dT03_mc_t, {diff(q1t(t), t), diff(q2t(t), t), diff(q3t(t), t)}, {dq1, dq2, dq3});
dT01_mc = subs(dT01_mc, {q1t, q2t, q3t}, {q1, q2, q3});
dT02_mc = subs(dT02_mc, {q1t, q2t, q3t}, {q1, q2, q3});
dT03_mc = subs(dT03_mc, {q1t, q2t, q3t}, {q1, q2, q3});
w01_mc = simplify(so3ToVec(dT01_mc(1:3, 1:3)*(T01_mc(1:3,1:3).')));
%dR01_mc = dT01_mc(1:3, 1:3);
v01_mc = dT01_mc(1:3, 4);
w02_mc = simplify(so3ToVec(dT02_mc(1:3, 1:3)*(T02_mc(1:3,1:3).')));
%dR02_mc = dT02_mc(1:3, 1:3);
v02_mc = dT02_mc(1:3, 4);
%dR03_mc = dT03_mc(1:3, 1:3);
w03_mc = simplify(so3ToVec(dT03_mc(1:3, 1:3)*(T03_mc(1:3,1:3).')));
v03_mc = dT03_mc(1:3, 4);

%%
% Kinetic energy
Ke = 1/2*m1*(v01_mc.')*v01_mc + 1/2*m2*(v02_mc.')*v02_mc + 1/2*m3*(v03_mc.')*v03_mc;
Ke = Ke + 1/2*(w01_mc.')*inertia_tensor2world(T01_mc, I1)*w01_mc +...
    1/2*(w02_mc.')*inertia_tensor2world(T02_mc, I2)*w02_mc +...
    1/2*(w03_mc.')*inertia_tensor2world(T03_mc, I3)*w03_mc;
Ke = simplify(Ke);

%%
% Potential energy
g = [0 0 -9.81];
Pe = simplify(dot(p01_mc, -g)*m1 + dot(p02_mc, -g)*m2 + dot(p03_mc, -g)*m3);

%%
% Lagrangian
L = Ke - Pe;

tau1 = subs(diff(subs(diff(L, dq1), {q1, q2, q3, dq1, dq2 ,dq3}, {q1t, q2t, q3t, diff(q1t, t), diff(q2t, t), diff(q3t, t)}), t),...
    {diff(q1t, t, 2), diff(q2t, t, 2), diff(q3t, t, 2), diff(q1t, t), diff(q2t, t), diff(q3t, t), q1t, q2t, q3t}, {ddq1, ddq2 ,ddq3, dq1, dq2 ,dq3, q1, q2, q3})...
    - diff(L, q1)...
    + Fv1*dq1;% + Fc1*sign(dq1) + Fo1;
% tau1 = simplify(tau1);
tau2 = subs(diff(subs(diff(L, dq2), {q1, q2, q3, dq1, dq2 ,dq3}, {q1t, q2t, q3t, diff(q1t, t), diff(q2t, t), diff(q3t, t)}), t),...
    {diff(q1t, t, 2), diff(q2t, t, 2), diff(q3t, t, 2), diff(q1t, t), diff(q2t, t), diff(q3t, t), q1t, q2t, q3t}, {ddq1, ddq2 ,ddq3, dq1, dq2 ,dq3, q1, q2, q3})...
    - diff(L, q2)...
    + Fv2*dq2;% + Fc2*sign(dq2) + Fo2;
% tau2 = simplify(tau2);
tau3 = subs(diff(subs(diff(L, dq3), {q1, q2, q3, dq1, dq2 ,dq3}, {q1t, q2t, q3t, diff(q1t, t), diff(q2t, t), diff(q3t, t)}), t),...
    {diff(q1t, t, 2), diff(q2t, t, 2), diff(q3t, t, 2), diff(q1t, t), diff(q2t, t), diff(q3t, t), q1t, q2t, q3t}, {ddq1, ddq2 ,ddq3, dq1, dq2 ,dq3, q1, q2, q3})...
    - diff(L, q3)...
    + Fv3*dq3;% + Fc3*sign(dq3) + Fo3;
% tau3 = simplify(tau3);

Tau = [tau1; tau2; tau3];

%%
% Write energy function in linear equation of inertia parameters
%X = [delta_L1; delta_L2; delta_L3];
X = [delta_L1; delta_A1; delta_L2; delta_A2; delta_L3; delta_A3];
% Energy
%h = equationsToMatrix(L, X);
% Torque
h = equationsToMatrix(Tau, X);
%%
% Calculate base parameters
% This method is refered to the following papar
% Gautier, Maxime. "Numerical calculation of the base inertial parameters of robots." Journal of Field Robotics 8.4 (1991): 485-506.

rand_var_file_name = "data/rand_var.mat";
rand_var.rand_num = length(h)+5;

if 2 == exist(rand_var_file_name)
    load(rand_var_file_name)
else
    % Generate random data
    rand_var.rand_num = length(h)+5;
    rand_var.q1_rand = (rand(rand_var.rand_num,1)-0.5)*6.28;
    rand_var.q2_rand = (rand(rand_var.rand_num,1)-0.5)*6.28;
    rand_var.q3_rand = (rand(rand_var.rand_num,1)-0.5)*6.28;
    rand_var.dq1_rand = (rand(rand_var.rand_num,1)-0.5)*6.28;
    rand_var.dq2_rand = (rand(rand_var.rand_num,1)-0.5)*6.28;
    rand_var.dq3_rand = (rand(rand_var.rand_num,1)-0.5)*6.28;
    rand_var.ddq1_rand = (rand(rand_var.rand_num,1)-0.5)*6.28;
    rand_var.ddq2_rand = (rand(rand_var.rand_num,1)-0.5)*6.28;
    rand_var.ddq3_rand = (rand(rand_var.rand_num,1)-0.5)*6.28;

    save(rand_var_file_name, "rand_var")
end

W = [];
for i=1:rand_var.rand_num
    W(i*3-2:i*3,:) = subs(h, {q1, q2, q3, dq1, dq2, dq3, ddq1, ddq2, ddq3},...
        {rand_var.q1_rand(i), rand_var.q2_rand(i), rand_var.q3_rand(i),...
        rand_var.dq1_rand(i), rand_var.dq2_rand(i), rand_var.dq3_rand(i),...
        rand_var.ddq1_rand(i), rand_var.ddq2_rand(i), rand_var.ddq3_rand(i)});
end

b = rank(W);
c = length(X);
% H*P = Q*R
[Q R P] = qr(W);
Xp = P'*X;
X1 = Xp(1:b,:);
X2 = Xp(b+1:end,:);
[Q, R] = qr(W*P);
R1 = R(1:b,1:b);
R2 = R(1:b, b+1:end);
% remove zero terms caused by computational precision
R1invR2 = inv(R1)*R2;
small_indices = find(abs(R1invR2) < 0.000001);
R1invR2(small_indices) = 0;
XB1 = X1 + R1invR2*X2;
% W1*XB1=K
% Numerical validation
% Or symbolic validation
%W = h;

Wp = W*P;
W1 = Wp(:,1:b);
W2 = Wp(:,b+1:end);
%%
% Validation of this reduction
K_err = simplify(W1*XB1-W*X);
[h_err] = equationsToMatrix(K_err, X);
vpa(h_err, 2)
disp("The number of errors which are larger than 0.000001 is:");
find(abs(h_err)>0.000001)


%% Optimal Trajctory Generation
h_b = h*P(:,1:b);
tr.h_b = h_b;
% Fundamental frequency
PI = 3.1415926;
w_f = 2*PI*0.1;
% Number of harmonics
n_H = 6;
%tr_file_name = "data/tr.mat"; % large workspace
tr_file_name = "data/tr_s2.mat"; % small workspace
regenerate_trajectory = 0;
% if file exists, load it; otherwise, compute one.
if 2 == exist(tr_file_name) && regenerate_trajectory == 0
    load(tr_file_name);
else
    tr = optimal_exciting_traj(h_b, n_H, w_f);
    save(tr_file_name, "tr")
end

plot_excitation_traj(tr);

%% Experiment data processing

%%
% Loading data
q1_file_name = "data/experiment_data/50s/outer_yaw_joint_states.csv";
q1_data_raw = csvread(q1_file_name);
q2_file_name = "data/experiment_data/50s/shoulder_pitch_joint_states.csv";
q2_data_raw = csvread(q2_file_name);
q3_file_name = "data/experiment_data/50s/elbow_pitch_joint_states.csv";
q3_data_raw = csvread(q3_file_name);

%%
% Get Acceleration by differentiation
sampling_freq = 200;
% d_t = 1/sampling_freq;
% dq_f = q1_data_raw(3:end, 2);
% dq_b = q1_data_raw(1:end-2, 2);
% 
% q1_data = zeros(size(q1_data_raw, 1)-2, size(q1_data_raw, 2)+1);
% q1_data(:,1:2) = q1_data_raw(2:end-1,1:2);
% q1_data(:,3) = (dq_f -dq_b)/(2*d_t);
% q1_data(:,4) = q1_data_raw(2:end-1,3);
q1_data = raw_data2data(q1_data_raw, sampling_freq);
q2_data = raw_data2data(q2_data_raw, sampling_freq);
q3_data = raw_data2data(q3_data_raw, sampling_freq);

%%
% remove abnormal acceleration data
% max_acc = 2.5;
% max_vel_change = 0.2;
% [q1_data, q2_data, q3_data] = remove_abnormal_acc_data(q1_data, q2_data, q3_data, max_acc, max_vel_change);

%%
% filter design
fc = 3;
fs = sampling_freq;

[b_f,a_f] = butter(10,fc/(fs/2));
freqz(b_f,a_f)

%%
% filt data
q1_data_filted = filt_data(q1_data, b_f, a_f);
q2_data_filted = filt_data(q2_data, b_f, a_f);
q3_data_filted = filt_data(q3_data, b_f, a_f);

plot_data(q1_data, q2_data, q3_data, q1_data_filted, q2_data_filted, q3_data_filted, sampling_freq);

%%
% remove near zero velocity data and outlier
vel_threshold = 0.00;
[q1_data_no_zero, q2_data_no_zero, q3_data_no_zero] = remove_near_zero_vel_data(q1_data_filted,...
    q2_data_filted, q3_data_filted, vel_threshold);

% q1_data_no_zero = q1_data_filted;
% q2_data_no_zero = q2_data_filted;
% q3_data_no_zero = q3_data_filted;
%%
% Generate regression matrix
[W_data, b_data] = generate_regression_mat(q1_data_no_zero, q2_data_no_zero, q3_data_no_zero, h_b);

%%
% least square
XB1_ols = W_data\b_data;
%%
% predict torque
predicted_vtau = W_data*XB1_ols;
l = size(predicted_vtau,1)/3;
predicted_tau1 = zeros(l,1);
predicted_tau2 = zeros(l,1);
predicted_tau3 = zeros(l,1);
for i = 1:l
    predicted_tau1(i) = predicted_vtau(3*(i-1)+1);
    predicted_tau2(i) = predicted_vtau(3*(i-1)+2);
    predicted_tau3(i) = predicted_vtau(3*(i-1)+3);
end
it =(1:l)/200.0;
figure
subplot(3,1,1)
plot(it, q1_data_no_zero(:,4), it, predicted_tau1);
xlabel("t (s)");
ylabel("Joint1 Torque (N*m)");
subplot(3,1,2)
plot(it, q2_data_no_zero(:,4), it, predicted_tau2);
xlabel("t (s)");
ylabel("Joint2 Torque (N*m)");
subplot(3,1,3)
plot(it, q3_data_no_zero(:,4), it, predicted_tau3);
xlabel("t (s)");
ylabel("Joint3 Torque (N*m)");
clear it l;
% %%
% % Generate regression matrix
% [W_raw_data, b_raw_data] = generate_regression_mat(q1_data, q2_data, q2_data, h_b);
% predicted_vtau_raw = W_raw_data*XB1_ols;
% l_raw = size(predicted_vtau_raw,1)/3;
% predicted_tau1_raw = zeros(l_raw,1);
% predicted_tau2_raw = zeros(l_raw,1);
% predicted_tau3_raw = zeros(l_raw,1);
% for i = 1:l_raw
%     predicted_tau1(i) = predicted_vtau(3*(i-1)+1);
%     predicted_tau2(i) = predicted_vtau(3*(i-1)+2);
%     predicted_tau3(i) = predicted_vtau(3*(i-1)+3);
% end
% it =1:l_raw;
% figure
% subplot(3,1,1)
% plot(it, q1_data_raw(:,4), it, predicted_tau1_raw);
% subplot(3,1,2)
% plot(it, q2_data_raw(:,4), it, predicted_tau2_raw);
% subplot(3,1,3)
% plot(it, q3_data_raw(:,4), it, predicted_tau3_raw);
% clear it l_raw;
%%
%
% variance of the regression error
var_reg_error_ols = norm(b_data - W_data*XB1_ols)/(length(b_data) - b);
% standard deviation of XB1_ols
std_XB1_ols = sqrt(diag(var_reg_error_ols*inv(W_data.'*W_data)));
% percentage of standard deviation of XB1_ols
pecent_std_XB1_ols = std_XB1_ols./abs(XB1_ols);

%%
% weighted least square

% variance of regression error of each joint
joint1_mask = 1:3:size(W_data,1);
joint2_mask = (1:3:size(W_data,1))+1;
joint3_mask = (1:3:size(W_data,1))+2;
var_error_ols_joint1 = norm(b_data(joint1_mask,:) - W_data(joint1_mask,:)*XB1_ols)...
    /(length(joint1_mask) - b);
var_error_ols_joint2 = norm(b_data(joint2_mask,:) - W_data(joint2_mask,:)*XB1_ols)...
    /(length(joint2_mask) - b);
var_error_ols_joint3 = norm(b_data(joint3_mask,:) - W_data(joint3_mask,:)*XB1_ols)...
    /(length(joint3_mask) - b);

W_data_weight = W_data;
W_data_weight(joint1_mask,:) = W_data(joint1_mask,:)/sqrt(var_error_ols_joint1);
W_data_weight(joint2_mask,:) = W_data(joint2_mask,:)/sqrt(var_error_ols_joint2);
W_data_weight(joint3_mask,:) = W_data(joint3_mask,:)/sqrt(var_error_ols_joint3);
b_data_weight = b_data;
b_data_weight(joint1_mask) = b_data(joint1_mask)/sqrt(var_error_ols_joint1);
b_data_weight(joint2_mask) = b_data(joint2_mask)/sqrt(var_error_ols_joint2);
b_data_weight(joint3_mask) = b_data(joint3_mask)/sqrt(var_error_ols_joint3);
clear joint1_mask joint2_mask joint3_mask;

XB1_wls = W_data_weight\b_data_weight;

%%
% predict torque
wls_predicted_vtau = W_data*XB1_wls;
l = size(wls_predicted_vtau,1)/3;
wls_predicted_tau1 = zeros(l,1);
wls_predicted_tau2 = zeros(l,1);
wls_predicted_tau3 = zeros(l,1);
for i = 1:l
    wls_predicted_tau1(i) = wls_predicted_vtau(3*(i-1)+1);
    wls_predicted_tau2(i) = wls_predicted_vtau(3*(i-1)+2);
    wls_predicted_tau3(i) = wls_predicted_vtau(3*(i-1)+3);
end
it =1:l;
figure
subplot(3,1,1)
plot(it, q1_data_no_zero(:,4), it, wls_predicted_tau1);
subplot(3,1,2)
plot(it, q2_data_no_zero(:,4), it, wls_predicted_tau2);
subplot(3,1,3)
plot(it, q3_data_no_zero(:,4), it, wls_predicted_tau3);
clear it l;

