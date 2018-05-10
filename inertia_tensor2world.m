function [Iworld] = inertia_tensor2world(T, I)
R = T(1:3, 1:3);
Iworld = R*I*R.';
end

