function [W, b] = generate_regression_mat(q1_data, q2_data, q3_data, h)
syms q1 q2 q3 q4 q5 q6 q7 real;
syms dq1 dq2 dq3 dq4 dq5 dq6 dq7 real;
syms ddq1 ddq2 ddq3 ddq4 ddq5 ddq6 ddq7 real;

l = size(q1_data, 1);
w = size(h,2);
W = zeros(l*3, w);
b = zeros(l*3, 1);
for i = 1:l
    W(3*(i-1)+1:3*(i-1)+3, :) = double(subs(h, {q1, q2, q3, dq1, dq2, dq3, ddq1, ddq2, ddq3},...
        {q1_data(i,1), q2_data(i,1), q3_data(i,1),...
        q1_data(i,2), q2_data(i,2), q3_data(i,2),...
        q1_data(i,3), q2_data(i,3), q3_data(i,3)}));
    b(3*(i-1)+1:3*(i-1)+3, :) = [q1_data(i,4); q2_data(i,4); q3_data(i,4)];
end
end

