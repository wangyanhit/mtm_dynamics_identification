function [tr] = optimal_exciting_traj(h, n_H, w_f)

tr.n_H = n_H;
tr.w_f = w_f;
tr.period = 2*3.141592653/tr.w_f;
tr.num = 8*n_H;
tr.dof = size(h,1);
tr.base_num = size(h,2);
tr.h = h;

% tr.max_q = 0.5;
% tr.min_q = -0.5;
q_delta = 0.2;
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
ab_init_max = 0.2;

% q01 = rand()*q_init_max-q_init_max/2;
% a1 = rand(1,tr.n_H)*ab_init_max-ab_init_max/2;
% b1 = rand(1,tr.n_H)*ab_init_max-ab_init_max/2;
% q02 = rand()*q_init_max-q_init_max/2;
% a2 = rand(1,tr.n_H)*ab_init_max-ab_init_max/2;
% b2 = rand(1,tr.n_H)*ab_init_max-ab_init_max/2;
% q03 = rand()*q_init_max-q_init_max/2;
% a3 = rand(1,tr.n_H)*ab_init_max-ab_init_max/2;
% b3 = rand(1,tr.n_H)*ab_init_max-ab_init_max/2;

q01 = -0.1590;
a1 = [0.4335 -0.0801 -0.3560 -0.0943 -0.1195 0.0041];
b1 = [-0.1509 -0.0845 0.0869 0.1800 0.1564 0.1114];
q02 = 0.0183;
a2 = [-0.0670 0.2883 -0.0050 -0.0592 -0.1569 -0.0273];
b2 = [0.3047 0.1555 0.0578 -0.0698 0.0652 -0.1081];
q03 = -0.0708;
a3 = [0.0996 -0.1643 -0.0212 0.0376 0.0936 -0.0192];
b3 = [0.3465 0.0534 -0.1563 -0.0571 -0.0031 0.1028];



x0 = [q01 a1 b1 q02 a2 b2 q03 a3 b3];
[A, b] = A_traj(tr);
Aeq = [];
beq = [];
ub = ones(1,length(x0))*5;
lb = -ub;



objfun = @(x)cond_traj(x, tr);
nonlcon = [];
%nonlcon = @(x)nonlcon_traj(x, h, tr);
% options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
options = optimoptions(@fmincon,'Display','iter','Algorithm','active-set','MaxFunctionEvaluations',20000, 'MaxIterations', 50);
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

