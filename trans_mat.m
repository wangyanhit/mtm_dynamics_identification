function [T] = trans_mat(v)
x = v(1);
y = v(2);
z = v(3);
T = [[1 0 0 x]; [0 1 0 y]; [0 0 1 z]; [0 0 0 1]];
end

