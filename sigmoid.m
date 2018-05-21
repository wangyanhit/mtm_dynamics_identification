function [y] = sigmoid(x, alpha, deadzone)

y= 2/(1+exp(-alpha*x)) - 1;
end

