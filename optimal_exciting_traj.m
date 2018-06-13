function [tr] = optimal_exciting_traj(h, n_H, w_f)

tr.n_H = n_H;
tr.w_f = w_f;
tr.period = 2*3.141592653/tr.w_f;
tr.num = 10*n_H;
tr.dof = size(h,1);
tr.base_num = size(h,2);
tr.h = h;

% large workspace
% q_delta = 0.25;
% tr.max_q1 = 0.785398 - q_delta;
% tr.min_q1 = -1.309 + q_delta;
% tr.max_q2 = 0.785398 - q_delta;
% tr.min_q2 = -1.785398 + q_delta;
% tr.max_q3 = 0.785398 - q_delta;
% tr.min_q3 = -1.785398 + q_delta;
% small workspace due to dvrk paper
q_delta = deg2rad(5);
tr.max_q1 = deg2rad(65) - q_delta;
tr.min_q1 = deg2rad(-40) + q_delta;
tr.max_q2 = deg2rad(50) - q_delta;
tr.min_q2 = deg2rad(-15) + q_delta;
tr.max_q3 = deg2rad(35) - q_delta;
tr.min_q3 = deg2rad(-50) + q_delta;


tr.max_dq = 3;
tr.min_dq = -3;
% The acceleration constraints are removed
tr.max_ddq = 10;
tr.min_ddq = -10;

q_init_max = 0.1;
ab_init_max = 0.2;

q01 = rand()*q_init_max-q_init_max/2;
a1 = rand(1,tr.n_H)*ab_init_max-ab_init_max/2;
b1 = rand(1,tr.n_H)*ab_init_max-ab_init_max/2;
q02 = rand()*q_init_max-q_init_max/2;
a2 = rand(1,tr.n_H)*ab_init_max-ab_init_max/2;
b2 = rand(1,tr.n_H)*ab_init_max-ab_init_max/2;
q03 = rand()*q_init_max-q_init_max/2;
a3 = rand(1,tr.n_H)*ab_init_max-ab_init_max/2;
b3 = rand(1,tr.n_H)*ab_init_max-ab_init_max/2;

% q01= -0.3909;
% a1= [0.0415 0.3569 0.1357 -0.0863 0.0283 -0.0791];
% b1= [-0.1051 0.7755 -0.1038 -0.6518 0.0793 0.2084];
% q02= 0.0426;
% a2= [-0.3960 -0.0669 -0.1852 -0.1514 0.0587 -0.1562];
% b2= [0.1303 -0.0565 0.2649 -0.0548 0.2778 0.0169];
% q03= -0.0544;
% a3= [-0.0232 0.0270 0.4741 0.3429 -0.1960 0.2062];
% b3= [-0.2055 0.0100 0.0572 -0.0024 -0.4890 -0.0768];
% q01 = -0.4151;
% a1 = [0.0401 0.2865 0.0508 -0.1680 0.0071 0.0136].*0.99;
% b1 = [-0.0530 0.8497 -0.1236 -0.5506 0.1712 0.0689].*0.99;
% q02 = 0.0417;
% a2 = [-0.3885 -0.0754 -0.1998 -0.1812 0.0637 -0.2131].*0.99;
% b2 = [0.1180 -0.0583 0.2368 -0.0449 0.2776 0.0300].*0.99;
% q03 = -0.0367;
% a3 = [-0.0693 -0.1555 0.3507 0.1921 -0.2333 0.2924].*0.99;
% b3 = [-0.2159 0.0032 0.0134 -0.0092 -0.4652 -0.0492].*0.99;
% large limit
% q01= -0.2565;
% a1= [-0.0432 0.0876 -0.0296 -0.1547 0.1460 0.1630];
% b1= [-0.1513 0.7446 -0.1736 -0.2973 -0.0745 -0.1692];
% q02= -0.3202;
% a2= [-0.6264 -0.0882 -0.4016 -0.0938 0.0178 -0.1788];
% b2= [0.0435 0.0206 0.0459 -0.1101 0.2631 -0.0097];
% q03= -0.3656;
% a3= [-0.0475 -0.2135 0.4993 0.0994 -0.2364 0.0822];
% b3= [-0.4528 -0.0563 0.1351 0.1387 -0.3307 -0.0984];


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
options = optimoptions(@fmincon,'Display','iter','Algorithm','active-set','MaxFunctionEvaluations',20000, 'MaxIterations', 60);
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

