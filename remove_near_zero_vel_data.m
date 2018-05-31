function [q1_out, q2_out] = remove_near_zero_vel_data(q1, q2, th)
mask1 = q1(:,2) > th | q1(:,2) < -th;
mask2 = q2(:,2) > th | q2(:,2) < -th;

mask = mask1 & mask2;
q1_out = q1(mask, :);
q2_out = q2(mask, :);
end

