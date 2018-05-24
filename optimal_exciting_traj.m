function [tr] = optimal_exciting_traj(h, n_H, w_f)

tr.n_H = n_H;
tr.w_f = w_f;
tr.period = 2*3.1415926/tr.w_f;
tr.num = 8*n_H;
tr.dof = size(h,1);
tr.base_num = size(h,2);
tr.h = h;

tr.max_q = 0.5;
tr.min_q = -0.5;
q_delta = 0.25;
tr.max_q1 = 0.785398 - q_delta;
tr.min_q1 = -1.309 + q_delta;
tr.max_q2 = 0.785398 - q_delta;
tr.min_q2 = -0.785398 + q_delta;
tr.max_q3 = 0.785398 - q_delta;
tr.min_q3 = -0.785398 + q_delta;


tr.max_dq = 2;
tr.min_dq = -2;
% The acceleration constraints are removed
tr.max_ddq = 10;
tr.min_ddq = -10;

q_init_max = 0.1;
ab_init_max = 0.1;

q01 = rand()*q_init_max-q_init_max/2;
a1 = rand(1,tr.n_H)*ab_init_max-ab_init_max/2;
b1 = rand(1,tr.n_H)*ab_init_max-ab_init_max/2;
q02 = rand()*q_init_max-q_init_max/2;
a2 = rand(1,tr.n_H)*ab_init_max-ab_init_max/2;
b2 = rand(1,tr.n_H)*ab_init_max-ab_init_max/2;
q03 = rand()*q_init_max-q_init_max/2;
a3 = rand(1,tr.n_H)*ab_init_max-ab_init_max/2;
b3 = rand(1,tr.n_H)*ab_init_max-ab_init_max/2;

% q01 = -0.2721;
% a1 = [0.0694 -0.3526 0.7608 0.1043 0.1654 0.0560]*0.99;
% b1 = [0.0138 0.3271 -0.0243 -0.3257 -0.1315 0.4086]*0.99;
% q02 = -0.0186;
% a2 = [0.0536 0.1499 -0.1590 -0.0349 -0.0703 -0.9585]*0.99;
% b2 = [0.0458 -0.2422 0.0560 -0.1084 0.0322 -0.1254]*0.99;
% q03 = 0.0379;
% a3 = [0.0642 0.1795 -0.3514 0.3555 0.1985 0.8235]*0.99;
% b3 = [0.0835 0.0456 -0.3276 0.0906 -0.1214 -0.3116]*0.99;


x0 = [q01 a1 b1 q02 a2 b2 q03 a3 b3];
[A, b] = A_traj(tr);
Aeq = [];
beq = [];
ub = ones(1,length(x0));
lb = -ub;



objfun = @(x)cond_traj(x, tr);
nonlcon = [];
%nonlcon = @(x)nonlcon_traj(x, h, tr);
% options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
options = optimoptions(@fmincon,'Display','iter','Algorithm','active-set','MaxFunctionEvaluations',20000, 'MaxIterations', 10);
[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(objfun,x0,A,b,Aeq,beq,lb,ub, nonlcon, options);

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

