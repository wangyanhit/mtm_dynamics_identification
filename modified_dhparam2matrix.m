function [T] = modified_dhparam2matrix(param)
a = param(1);
alpha = param(2);
d = param(3);
theta = param(4);

T = [[cos(theta), -sin(theta), 0, a];
    [sin(theta)*cos(alpha), cos(theta)*cos(alpha), -sin(alpha), -d*sin(alpha)];
    [sin(theta)*sin(alpha), cos(theta)*sin(alpha), cos(alpha), cos(alpha)*d];
    [0, 0, 0, 1]];
T = simplify(T);
end

