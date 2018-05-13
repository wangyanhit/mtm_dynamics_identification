clc;clear;
A = [];
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];
x0 = [0, 0];
objfun = @(x)fun(x,1);
nonlc = @(x)nonlcon(x,1);
[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(objfun,x0,A,b,Aeq,beq,lb,ub,nonlc)

