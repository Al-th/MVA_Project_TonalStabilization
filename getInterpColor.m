function [c1,c2,c3] = getInterpColor(xi,xj,filteredFrame2,A_init,R,sigma1,sigma2,sigma3)
[height,width,~] = size(filteredFrame2);
num1 = 0;
denum1 = 0;

num2 = 0;
denum2 = 0;

num3 = 0;
denum3 = 0;
for i = 1:height
    for j = 1:width
    [w1,w2,w3] = getAffinity(filteredFrame2(xi,xj,:),filteredFrame2(i,j,:),sigma1,sigma2,sigma3);
    num1 = num1 + w1*A_init(i,j,1);
    denum1 = denum1 + w1*R(i,j);
    
    num2 = num2 + w2*A_init(i,j,2);
    denum2 = denum2 + w2*R(i,j);
    
    num3 = num3 + w3*A_init(i,j,3);
    denum3 = denum3 + w3*R(i,j);
    end
end

c1 = num1/denum1;
c2 = num2/denum2;
c3 = num3/denum3;

end