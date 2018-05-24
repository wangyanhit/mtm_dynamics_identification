function [T] = dhparam2matrix(param)
a = param(1);
alpha = param(2);
d = param(3);
theta = param(4);

T = [[cos(theta), -sin(theta)*cos(alpha), sin(theta)*sin(alpha), a*cos(theta)];
    [sin(theta), cos(theta)*cos(alpha), -cos(theta)*sin(alpha), a*sin(theta)];
    [0, sin(alpha), cos(alpha), d];
    [0, 0, 0, 1]];
T = simplify(T);
end

