function [q, dq, ddq] = fourier_func(q0, a, b, w_f, t)
q = q0;
dq = 0;
ddq = 0;
n_H = length(a);
for i =1:n_H
    q = q + a(i)/(w_f*i)*sin(w_f*i*t) - b(i)/(w_f*i)*cos(w_f*i*t);
    dq = dq + a(i)*cos(w_f*i*t) + b(i)*sin(w_f*i*t);
    ddq = ddq - a(i)*(w_f*i)*sin(w_f*i*t) + b(i)*(w_f*i)*cos(w_f*i*t); 
end
end

