function [q1_out, q2_out, q3_out] = remove_abnormal_acc_data(q1, q2, q3, th, max_dv)
mask1 = q1(:,3) < th & q1(:,3) > -th;
mask2 = q2(:,3) < th & q2(:,3) > -th;
mask3 = q3(:,3) < th & q3(:,3) > -th;

for i = 1:(size(mask1,1)-1)
    if q1(i+1,2) > q1(i,2) + max_dv | q1(i+1,2) < q1(i,2) - max_dv |...
            q2(i+1,2) > q2(i,2) + max_dv | q2(i+1,2) < q2(i,2) - max_dv |...
            q3(i+1,2) > q3(i,2) + max_dv | q3(i+1,2) < q3(i,2) - max_dv
        mask1(i,1) = 0;
    end
            
end

mask = mask1 & mask2 & mask3;
q1_out = q1(mask, :);
q2_out = q2(mask, :);
q3_out = q3(mask, :);
end

