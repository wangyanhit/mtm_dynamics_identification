function [A, b] = A_traj(tr)

f_num = tr.n_H*2 + 1;
l_v = tr.dof*f_num;
interval = tr.period/tr.num;
vt = 0:interval:tr.period;

A = [];
b = [];

for t = vt
    qi = zeros(3,l_v);
    dqi = zeros(3,l_v);
    ddqi = zeros(3,l_v);
    
    q_tmp = zeros(1,f_num);
    q_tmp(1,1) = 1;
    for i = 1:tr.n_H
        q_tmp(1,i+1) = 1/(tr.w_f*i)*sin(tr.w_f*i*t);
        q_tmp(1,i+1+tr.n_H) = -1/(tr.w_f*i)*cos(tr.w_f*i*t);
    end
    qi(1,1:f_num) = q_tmp;
    qi(2,1+f_num:2*f_num) = q_tmp;
    qi(3,1+2*f_num:3*f_num) = q_tmp;

    q_tmp(1,1) = 0;
    for i = 1:tr.n_H
        q_tmp(1,i+1) = cos(tr.w_f*i*t);
        q_tmp(1,i+1+tr.n_H) = sin(tr.w_f*i*t);
    end
    dqi(1,1:f_num) = q_tmp;
    dqi(2,1+f_num:2*f_num) = q_tmp;
    dqi(3,1+2*f_num:3*f_num) = q_tmp;
    
    q_tmp(1,1) = 0;
    for i = 1:tr.n_H
        q_tmp(1,i+1) = -tr.w_f*i*sin(tr.w_f*i*t);
        q_tmp(1,i+1+tr.n_H) = tr.w_f*i*cos(tr.w_f*i*t);
    end
    ddqi(1,1:f_num) = q_tmp;
    ddqi(2,1+f_num:2*f_num) = q_tmp;
    ddqi(3,1+2*f_num:3*f_num) = q_tmp;
    
%     A = [A; qi; dqi; ddqi; -qi; -dqi; -ddqi];
%     b_tmp1 = [tr.max_q1; tr.max_q2; tr.max_q3; tr.max_dq; tr.max_dq; tr.max_dq; tr.max_ddq; tr.max_ddq; tr.max_ddq;];
%     b_tmp2 = [-tr.min_q1; -tr.min_q2; -tr.min_q3; tr.max_dq; tr.max_dq; tr.max_dq; tr.max_ddq; tr.max_ddq; tr.max_ddq;];
%     b = [b; b_tmp1; b_tmp2];

    % without acceleration limitation
    A = [A; qi; dqi; -qi; -dqi];
    b_tmp1 = [tr.max_q1; tr.max_q2; tr.max_q3; tr.max_dq; tr.max_dq; tr.max_dq;];
    b_tmp2 = [-tr.min_q1; -tr.min_q2; -tr.min_q3; tr.max_dq; tr.max_dq; tr.max_dq;];
    b = [b; b_tmp1; b_tmp2];
end
end

