Fv1_sim = 0.1;
Fv2_sim = 0.1;
Fv3_sim = 0.1;


Lyy1_sim = 0.005;


m2_sim = 2;
Lxx2_sim = 0.0002;
Lyy2_sim = 0.02;
Lzz2_sim = 0.02;
Lxz2_sim = 0;
Lxy2_sim = 0;
Lyz2_sim = 0;

lx2_sim = -0.1;
ly2_sim = 0;
lz2_sim = 0;

m3_sim = 1;
Lxx3_sim = 0.0001;
Lyy3_sim = 0.01;
Lzz3_sim = 0.01;
Lxz3_sim = 0;
Lxy3_sim = 0;
Lyz3_sim = 0;

lx3_sim = -0.15;
ly3_sim = 0;
lz3_sim = 0;

XB1_num = subs(XB1, {Fv1, Fv2, Fv3, Lyy1,...
Lxx2, Lxy2, Lxz2, Lyy2, Lyz2, Lzz2, m2, lx2, ly2, lz2,...
Lxx3, Lxy3, Lxz3, Lyy3, Lyz3, Lzz3, m3, lx3, ly3, lz3},...
{Fv1_sim, Fv2_sim, Fv3_sim, Lyy1_sim,...
Lxx2_sim, Lxy2_sim, Lxz2_sim, Lyy2_sim, Lyz2_sim, Lzz2_sim, m2_sim, lx2_sim, ly2_sim, lz2_sim,...
Lxx3_sim, Lxy3_sim, Lxz3_sim, Lyy3_sim, Lyz3_sim, Lzz3_sim, m3_sim, lx3_sim, ly3_sim, lz3_sim});
XB1_num = double(XB1_num)

sim.tf = 10;
sim.freq_samp = 100;
sim.dt = 1/sim.freq_samp;
sim.q = zeros(3,1);
sim.dq = zeros(3,1);
sim.last_dq = zeros(3,1);
sim.ddq = zeros(3,1);
sim.tau = zeros(3,1);

sim.Kp = ones(3)*5;
sim.Kd = ones(3)*1;
sim.v_l = sim.tf*sim.freq_samp;
sim.v_q = zeros(sim.v_l, 3);
sim.v_qd = zeros(sim.v_l, 3);
sim.v_tau = zeros(sim.v_l, 3);
sim.cnt = 0;
h_b_xb1 = h_b*XB1_num;
for t = 0:sim.dt:sim.tf
    sim.cnt = sim.cnt+1
    
    sim.qd = fourier_func(tr.q01, tr.a1, tr.b1, tr.w_f, t);
    sim.qd(2,1) = fourier_func(tr.q02, tr.a2, tr.b2, tr.w_f, t);
    sim.qd(3,1) = fourier_func(tr.q03, tr.a3, tr.b3, tr.w_f, t);
    
    
    sim.tau = sim.Kp*(sim.qd - sim.q) + sim.Kd*(sim.last_dq - sim.dq);
    sim.last_dq = sim.dq;
    
    ddq = solve(subs(h_b_xb1, {q1, q2, q3, dq1, dq2, dq3},...
        {sim.q(1), sim.q(2), sim.q(3), sim.dq(1), sim.dq(2), sim.dq(3)}) == sim.tau,...
        [ddq1, ddq2, ddq3]);
    sim.ddq = double([ddq.ddq1; ddq.ddq2; ddq.ddq3]);
    sim.dq = sim.dq + sim.ddq*sim.dt;
    sim.q = sim.q + sim.dq*sim.dt;
    
    sim.v_q(sim.cnt, :) = sim.q.';
    sim.v_qd(sim.cnt, :) = sim.qd.';
    sim.v_tau(sim.cnt, :) = sim.tau.';
end

