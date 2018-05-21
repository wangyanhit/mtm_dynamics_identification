n = 10;
m = 2;
Wn = 0.2;

[b,a] = maxflat(n,m,Wn);
fvtool(b,a)