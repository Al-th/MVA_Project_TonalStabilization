function [a1,a2,a3] = getAffinity(p1,p2,sigma1,sigma2,sigma3)
a1 = exp(-(abs(p1(1)-p2(1))^2)/(2*(sigma1^2)));
a2 = exp(-(abs(p1(2)-p2(2))^2)/(2*(sigma2^2)));
a3 = exp(-(abs(p1(3)-p2(3))^2)/(2*(sigma3^2)));
end