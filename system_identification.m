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
        T_num = subs(T(:,:,idx), {q1 q2 q3 q4 q5 q6 q7}, {0 0.2 0 0 0.2 0.0 0.4});
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
[delta_L1, r1, I1] = inertia_bary2std(Lxx1, Lxy1, Lxz1, Lyy1, Lyz1, Lzz1, lx1, ly1, lz1, m1);
syms Lxx2 Lxy2 Lxz2 Lyy2 Lyz2 Lzz2 lx2 ly2 lz2 m2 real;
[delta_L2, r2, I2] = inertia_bary2std(Lxx2, Lxy2, Lxz2, Lyy2, Lyz2, Lzz2, lx2, ly2, lz2, m2);
syms Lxx3 Lxy3 Lxz3 Lyy3 Lyz3 Lzz3 lx3 ly3 lz3 m3 real;
[delta_L3, r3, I3] = inertia_bary2std(Lxx3, Lxy3, Lxz3, Lyy3, Lyz3, Lzz3, lx3, ly3, lz3, m3);

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
w01_mc = simplify(so3ToVec(dT01_mc(1:3, 1:3)*T01_mc(1:3,1:3).'));
%dR01_mc = dT01_mc(1:3, 1:3);
v01_mc = dT01_mc(1:3, 4);
w02_mc = simplify(so3ToVec(dT02_mc(1:3, 1:3)*T02_mc(1:3,1:3).'));
%dR02_mc = dT02_mc(1:3, 1:3);
v02_mc = dT02_mc(1:3, 4);
%dR03_mc = dT03_mc(1:3, 1:3);
w03_mc = simplify(so3ToVec(dT03_mc(1:3, 1:3)*T03_mc(1:3,1:3).'));
v03_mc = dT03_mc(1:3, 4);

%%
% Kinetic energy
Ke = 1/2*m1*v01_mc.'*v01_mc + 1/2*m2*v02_mc.'*v02_mc + 1/2*m3*v03_mc.'*v03_mc;
Ke = Ke + 1/2*w01_mc.'*inertia_tensor2world(T01_mc, I1)*w01_mc +...
    1/2*w02_mc.'*inertia_tensor2world(T02_mc, I2)*w02_mc +...
    1/2*w03_mc.'*inertia_tensor2world(T03_mc, I3)*w03_mc;
Ke = simplify(Ke);

%%
% Potential energy
g = [0 0 -0.981];
Pe = simplify(dot(p01_mc, -g)*m1 + dot(p02_mc, -g)*m2 + dot(p03_mc, -g)*m3);

%%
% Lagrangian
L = Ke - Pe;

tau1 = subs(diff(subs(diff(L, dq1), {q1, q2, q3, dq1, dq2 ,dq3}, {q1t, q2t, q3t, diff(q1t, t), diff(q2t, t), diff(q3t, t)}), t),...
    {diff(q1t, t, 2), diff(q2t, t, 2), diff(q3t, t, 2), diff(q1t, t), diff(q2t, t), diff(q3t, t), q1t, q2t, q3t}, {ddq1, ddq2 ,ddq3, dq1, dq2 ,dq3, q1, q2, q3})...
    - diff(L, q1);
% tau1 = simplify(tau1);
tau2 = subs(diff(subs(diff(L, dq2), {q1, q2, q3, dq1, dq2 ,dq3}, {q1t, q2t, q3t, diff(q1t, t), diff(q2t, t), diff(q3t, t)}), t),...
    {diff(q1t, t, 2), diff(q2t, t, 2), diff(q3t, t, 2), diff(q1t, t), diff(q2t, t), diff(q3t, t), q1t, q2t, q3t}, {ddq1, ddq2 ,ddq3, dq1, dq2 ,dq3, q1, q2, q3})...
    - diff(L, q2);
% assume(m2 ~= 0);
% assume(m3 ~= 0);
% tau2 = simplify(tau2);
tau3 = subs(diff(subs(diff(L, dq3), {q1, q2, q3, dq1, dq2 ,dq3}, {q1t, q2t, q3t, diff(q1t, t), diff(q2t, t), diff(q3t, t)}), t),...
    {diff(q1t, t, 2), diff(q2t, t, 2), diff(q3t, t, 2), diff(q1t, t), diff(q2t, t), diff(q3t, t), q1t, q2t, q3t}, {ddq1, ddq2 ,ddq3, dq1, dq2 ,dq3, q1, q2, q3})...
    - diff(L, q3);
% assume(m2 ~= 0);
% assume(m3 ~= 0);
% tau3 = simplify(tau3);

Tau = [tau1; tau2; tau3];
%%
% Write energy function in linear equation of inertia parameters
X = [delta_L1; delta_L2; delta_L3];
% Energy
%h = equationsToMatrix(L, X);
% Torque
h = equationsToMatrix(Tau, X);
%%
% Calculate base parameters
% This method is refered to the following papar
% Gautier, Maxime. "Numerical calculation of the base inertial parameters of robots." Journal of Field Robotics 8.4 (1991): 485-506.
% Generate random data
rand_num = length(h)+5;
q1_rand = (rand(rand_num,1)-0.5)*6.28;
q2_rand = (rand(rand_num,1)-0.5)*6.28;
q3_rand = (rand(rand_num,1)-0.5)*6.28;
dq1_rand = (rand(rand_num,1)-0.5)*6.28;
dq2_rand = (rand(rand_num,1)-0.5)*6.28;
dq3_rand = (rand(rand_num,1)-0.5)*6.28;
ddq1_rand = (rand(rand_num,1)-0.5)*6.28;
ddq2_rand = (rand(rand_num,1)-0.5)*6.28;
ddq3_rand = (rand(rand_num,1)-0.5)*6.28;
W = [];
for i=1:rand_num
    W(i*3-2:i*3,:) = subs(h, {q1, q2, q3, dq1, dq2, dq3, ddq1, ddq2, ddq3}, {q1_rand(i), q2_rand(i), q3_rand(i), dq1_rand(i), dq2_rand(i), dq3_rand(i), ddq1_rand(i), ddq2_rand(i), ddq3_rand(i)});
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
w_f = 2*pi*0.1;
% Number of harmonics
n_H = 4;
tr = optimal_exciting_traj(h_b, n_H, w_f);



%%
% Visualization to see if the transformations are right
% figure
% hold on
% for link = fieldnames(dh_params) 
%     for idx = 1:length(link)
%         T = dh_params.(link{idx}).T_base2current;
%         T = subs(T, {q1 q2 q3 q4 q5 q6 q7}, {0 0 0 0 0 0 0});
%         p = T*([0 0 0 1]');
%         x = T*([0.1 0 0 1]');
%         y = T*([0 0.1 0 1]');
%         z = T*([0 0 0.1 1]');
%         plot3([p(1) x(1)], [p(2) x(2)], [p(3) x(3)], 'r')
%         plot3([p(1) y(1)], [p(2) y(2)], [p(3) y(3)], 'g')
%         plot3([p(1) z(1)], [p(2) z(2)], [p(3) z(3)], 'b')
%         text(p(1), p(2), p(3), num2str(idx));
%     end
% end




% frame_names = ["1" "2" "3" "4" "5" "6" "7" "2p" "2pp"];

% dh_params.link1.param = [0 pi/2 0 q1];
% dh_params.link1.prev = '0';
% dh_params.link1.prev_link = 0;
% dh_params.link2.param = [0.279 0 0 q2];
% dh_params.link2.prev = '1';
% dh_params.link2.prev_link = dh_params.link1;
% dh_params.link3.param = [0.365 -pi/2 0 q3];
% dh_params.link3.prev = '2';
% dh_params.link3.prev_link = dh_params.link2;
% dh_params.link4.param = [0 pi/2 0.151 q4];
% dh_params.link4.prev = '3';
% dh_params.link4.prev_link = dh_params.link3;
% dh_params.link5.param = [0 -pi/2 0 q5];
% dh_params.link5.prev = '4';
% dh_params.link5.prev_link = dh_params.link4;
% dh_params.link6.param = [0 pi/2 0 q6];
% dh_params.link6.prev = '5';
% dh_params.link6.prev_link = dh_params.link5;
% dh_params.link7.param = [0 0 0 q7];
% dh_params.link7.prev = '6';
% dh_params.link7.prev_link = dh_params.link6;
% dh_params.link2p.param = [0.1 0 0 q2p];
% dh_params.link2p.prev = '1';
% dh_params.link2p.prev_link = dh_params.link1;
% dh_params.link2pp.param = [0.279 0 0 q2pp];
% dh_params.link2pp.prev = '2p';
% dh_params.link2pp.prev_link = dh_params.link2p;
% dh_params = [['1' '0' 0 pi/2 0 q1];...
%             ['2' '1' 0.279 0 0 q2];...
%             ['2p' '1' 0.1 0 0 q2p];...
%             ['2pp' '2p' 0.279 0 0 q2pp];...
%             ['3' '2' 0.365 -pi/2 0 q3];...
%             ['4' '3' 0 pi/2 0.151 q4];...
%             ['5' '4' 0 -pi/2 0 q5];...
%             ['6' '5' 0 pi/2 0 q6];...
%             ['7' '6' 0 0 0 q7]];

%% 
% Transformations from frame i-1 to frame i and from base to i
% for link = fieldnames(dh_params)
%     for idx = 1:length(link)
%         dh_params.(link{idx}).T_prev2current = dhparam2matrix(dh_params.(link{idx}).param);
%         if idx == 1
%             dh_params.(link{idx}).T_base2current = dh_params.(link{idx}).T_prev2current;
%         else
%             prev_idx = find(frame_names == dh_params.(link{idx}).prev);
%             dh_params.(link{idx}).T_base2current = dh_params.(link{prev_idx}).T_base2current*dh_params.(link{idx}).T_prev2current;
%         end
%     end
% end
%%
% Visualization to see if the transformations are right
% figure
% hold on
% for link = fieldnames(dh_params) 
%     for idx = 1:length(link)
%         T = dh_params.(link{idx}).T_base2current;
%         T = subs(T, {q1 q2 q3 q4 q5 q6 q7}, {0 0 0 0 0 0 0});
%         p = T*([0 0 0 1]');
%         x = T*([0.1 0 0 1]');
%         y = T*([0 0.1 0 1]');
%         z = T*([0 0 0.1 1]');
%         plot3([p(1) x(1)], [p(2) x(2)], [p(3) x(3)], 'r')
%         plot3([p(1) y(1)], [p(2) y(2)], [p(3) y(3)], 'g')
%         plot3([p(1) z(1)], [p(2) z(2)], [p(3) z(3)], 'b')
%         text(p(1), p(2), p(3), num2str(idx));
%     end
% end

% T01 = dhparam2matrix(dh_link1);
% T12 = dhparam2matrix(dh_link2);
% T23 = dhparam2matrix(dh_link3);
% T34 = dhparam2matrix(dh_link4);
% T45 = dhparam2matrix(dh_link5);
% T56 = dhparam2matrix(dh_link6);
% T67 = dhparam2matrix(dh_link7);
% T12p = dhparam2matrix(dh_link2p);
% T2p2pp = dhparam2matrix(dh_link2pp);
% %%
% % Transformation from 
% T02 = T01*T12;
% T03 = T02*T23;
% T04 = T03*T34;
% T05 = T04*T45;
% T06 = T05*T56;
% T07 = T06*T67;
% T02p = T01*T12p;
% T02pp = T02p*T2p2pp;
%%
% Define friction types

%%


%% Dynamic model generation


%% Calclualte base parameters

% QR factorization with colomn pivoting
%[Q R P] = qr()

%% Data processing

%% Model identification

