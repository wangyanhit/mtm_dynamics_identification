function [tr] = optimal_exciting_traj(h, n_H, w_f)

tr.n_H = n_H;
tr.w_f = w_f;
tr.period = 2*pi/tr.w_f;
tr.num = 20;
tr.dof = size(h,1);
tr.max_q = 0.5;
tr.min_q = -0.5;
tr.max_dq = 2;
tr.min_dq = -2;
tr.max_ddq = 10;
tr.min_ddq = -10;


q01 = 0;
a1 = ones(1,tr.n_H)*0.2;
b1 = ones(1,tr.n_H)*0.2;
q02 = 0;
a2 = ones(1,tr.n_H)*0.2;
b2 = ones(1,tr.n_H)*0.2;
q03 = 0;
a3 = ones(1,tr.n_H)*0.2;
b3 = ones(1,tr.n_H)*0.2;


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
options = optimoptions(@fmincon,'Display','iter','Algorithm','active-set','MaxFunctionEvaluations',20000);
[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(objfun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)


%cond_traj(x0, h, tr)

end

