function w = windowMS(x, alpha)
  w = exp(-alpha.*x.*x) + exp(-alpha.*(x+2).*(x+2)) + exp(-alpha.*(x-2).*(x-2)) ...
     - (2*exp(-alpha)+exp(-9*alpha));
end