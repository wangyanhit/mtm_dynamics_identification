function [delta_L, r, I] = inertia_bary2std(Lxx, Lxy, Lxz, Lyy, Lyz, Lzz, lx, ly, lz, m)
delta_L = [Lxx Lxy Lxz Lyy Lyz Lzz lx ly lz m].';
l = [lx ly lz].';
r = l./m;
Lv = [Lxx, Lxy, Lxz, Lyy, Lyz, Lzz];
L = inertia_vec2mat(Lv);
I = L - m*VecToso3(r).'*VecToso3(r);
end

