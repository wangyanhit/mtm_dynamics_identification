function [c, ceq] = nonlcon_traj(v, h, tr)
syms q1 q2 q3 q4 q5 q6 q7 real;
syms dq1 dq2 dq3 dq4 dq5 dq6 dq7 real;
syms ddq1 ddq2 ddq3 ddq4 ddq5 ddq6 ddq7 real;

f_num = tr.n_H*2 + 1;
qi = [v(1); v(f_num+1); v(f_num*2+1)];
dqi = zeros(3,1);
ddqi = zeros(3,1);
H = [];
size(v);
interval = tr.period/tr.num;
vt = 0:interval:tr.period;

c = [];
for t = vt
    for j = 1:tr.dof
        for i =1:tr.n_H
            qi(j) = qi(j)+v((j-1)*f_num+1+i)/(tr.w_f*i)*sin(tr.w_f*i*t)-v((j-1)*f_num+tr.n_H+i)/(tr.w_f*i)*cos(tr.w_f*i*t);
            dqi(j) = dqi(j)+v((j-1)*f_num+1+i)*cos(tr.w_f*i*t)+v((j-1)*f_num+tr.n_H+i)*sin(tr.w_f*i*t);
            ddqi(j) = ddqi(j)-v((j-1)*f_num+1+i)*(tr.w_f*i)*sin(tr.w_f*i*t)+v((j-1)*f_num+tr.n_H+i)*(tr.w_f*i)*cos(tr.w_f*i*t);    
        end
    end
    c = [c; qi-tr.max_q; dqi-tr.max_dq; ddqi-tr.max_ddq; -qi+tr.min_q; -dqi+tr.min_dq; -ddqi+tr.min_ddq;];
end
c = double(c);
ceq =[];
end

