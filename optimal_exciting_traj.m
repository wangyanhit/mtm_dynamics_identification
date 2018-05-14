function [tr] = optimal_exciting_traj(h, n_H, w_f)

tr.n_H = n_H;
tr.w_f = w_f;
tr.period = 2*pi/tr.w_f;
tr.num = 2^(n_H+1);
tr.dof = size(h,1);

tr.max_q = 0.5;
tr.min_q = -0.5;
tr.max_dq = 2;
tr.min_dq = -2;
tr.max_ddq = 5;
tr.min_ddq = -5;

q_init_max = 0.1;
ab_init_max = 0.05;

q01 = rand()*q_init_max-q_init_max/2;
a1 = rand(1,tr.n_H)*ab_init_max-ab_init_max/2;
b1 = rand(1,tr.n_H)*ab_init_max-ab_init_max/2;
q02 = rand()*q_init_max-q_init_max/2;
a2 = rand(1,tr.n_H)*ab_init_max-ab_init_max/2;
b2 = rand(1,tr.n_H)*ab_init_max-ab_init_max/2;
q03 = rand()*q_init_max-q_init_max/2;
a3 = rand(1,tr.n_H)*ab_init_max-ab_init_max/2;
b3 = rand(1,tr.n_H)*ab_init_max-ab_init_max/2;


x0 = [q01 a1 b1 q02 a2 b2 q03 a3 b3];
A = [];
b = [];
Aeq = [];
beq = [];
ub = ones(1,length(x0));
lb = -ub;



objfun = @(x)cond_traj(x, h, tr);
nonlcon = @(x)nonlcon_traj(x, h, tr);
% options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
options = optimoptions(@fmincon,'Display','iter','Algorithm','active-set','MaxFunctionEvaluations',20000, 'MaxIterations', 15);
[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(objfun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);

tr.q01 = x(1);
tr.a1 = x(2:1+tr.n_H);
tr.b1 = x(2+tr.n_H:1+2*tr.n_H);
b_c = 1+2*tr.n_H;
tr.q02 = x(1+b_c);
tr.a2 = x(2+b_c:1+tr.n_H+b_c);
tr.b2 = x(2+tr.n_H+b_c:1+2*tr.n_H+b_c);
b_c = (1+2*tr.n_H)*2;
tr.q03 = x(1+b_c);
tr.a3 = x(2+b_c:1+tr.n_H+b_c);
tr.b3 = x(2+tr.n_H+b_c:1+2*tr.n_H+b_c);
%cond_traj(x0, h, tr)


end

