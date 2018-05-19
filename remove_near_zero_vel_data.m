function [q1_out, q2_out, q3_out] = remove_near_zero_vel_data(q1, q2, q3, th)
mask1 = q1(:,2) > th | q1(:,2) < -th;
mask2 = q2(:,2) > th | q2(:,2) < -th;
mask3 = q3(:,2) > th | q3(:,2) < -th;

mask = mask1 & mask2 & mask3;
q1_out = q1(mask, :);
q2_out = q2(mask, :);
q3_out = q3(mask, :);
end

