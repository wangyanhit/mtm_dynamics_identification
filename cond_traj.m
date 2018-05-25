function [cond_num] = cond_traj(v, tr)
syms q1 q2 q3 real;
syms dq1 dq2 dq3 real;
syms ddq1 ddq2 ddq3 real;

% q01 = v(1);
% a1 = v(1+1:1+tr.n_H);
% b1 = v(2+tr.n_H:1+2*tr.n_H);
% b = 1+2*tr.n_H;
% q02 = v(1+b);
% a2 = v(1+1+b:1+tr.n_H+b);
% b2 = v(2+tr.n_H+b:1+2*tr.n_H+b);
% b = (1+2*tr.n_H)*2;
% q02 = v(1+b);
% a2 = v(1+1+b:1+tr.n_H+b);
% b2 = v(2+tr.n_H+b:1+2*tr.n_H+b);
% q03 = v(1);
% a3 = v(1+1:1+tr.n_H);
% b3 = v(2+tr.n_H:1+2*tr.n_H);

f_num = tr.n_H*2 + 1;

%H = [];

size(v);
interval = tr.period/tr.num;
vt = 0:interval:tr.period;
t_num = length(vt);
H = zeros(tr.dof*t_num, tr.base_num);

for k = 1:t_num
    t = vt(k);
    qi = [v(1); v(f_num+1); v(f_num*2+1)];
    dqi = zeros(3,1);
    ddqi = zeros(3,1);
    for j = 1:tr.dof
        for i =1:tr.n_H
            qi(j) = qi(j)+v((j-1)*f_num+1+i)/(tr.w_f*i)*sin(tr.w_f*i*t)...
                -v((j-1)*f_num+tr.n_H+1+i)/(tr.w_f*i)*cos(tr.w_f*i*t);
            dqi(j) = dqi(j)+v((j-1)*f_num+1+i)*cos(tr.w_f*i*t)...
                +v((j-1)*f_num+tr.n_H+1+i)*sin(tr.w_f*i*t);
            ddqi(j) = ddqi(j)-v((j-1)*f_num+1+i)*(tr.w_f*i)*sin(tr.w_f*i*t)...
                +v((j-1)*f_num+tr.n_H+1+i)*(tr.w_f*i)*cos(tr.w_f*i*t);    
        end
    end
%     H = [H;subs(h, {q1, q2, q3, dq1, dq2, dq3, ddq1, ddq2, ddq3},...
%         {qi(1), qi(2), qi(3), dqi(1), dqi(2), dqi(3), ddqi(1), ddqi(2), ddqi(3)})];
%     H = [H; h_b_func(qi(1), qi(2), qi(3), dqi(1), dqi(2), dqi(3), ddqi(1), ddqi(2), ddqi(3))];
    H(3*(k-1)+1:3*(k-1)+3,:) = double(subs(tr.h, {q1, q2, q3, dq1, dq2, dq3, ddq1, ddq2, ddq3}, {qi(1), qi(2), qi(3), dqi(1), dqi(2), dqi(3), ddqi(1), ddqi(2), ddqi(3)}));% h_b_func(qi(1), qi(2), qi(3), dqi(1), dqi(2), dqi(3), ddqi(1), ddqi(2), ddqi(3));
end
%vpa(H,2)
%size(H)
%H = double(H);
cond_num = cond(H);

end

