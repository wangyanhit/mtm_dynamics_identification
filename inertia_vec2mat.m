function [I] = inertia_vec2mat(vec)
Ixx = vec(1);
Ixy = vec(2);
Ixz = vec(3);
Iyy = vec(4);
Iyz = vec(5);
Izz = vec(6);

I = [[Ixx Ixy Ixz]; [Ixy Iyy Iyz]; [Ixz Iyz Izz]];
end

