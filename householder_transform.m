% 第87页 计算Householder变换
function [v,beta,H]=householder_transform(x)
% clear;clc;
% x=rand(10,1);
n=length(x);
eta=norm(x,inf);
x=x/eta;
% sigma=x(2:n)'*x(2:n);
sigma=0;
for i=2:n
    sigma=sigma+x(i)^2;
end
v=x;
xx=x;
if sigma==0
    beta=0;
else
    alpha=sqrt(x(1)^2+sigma);
    if x(1)<=0
        v(1)=x(1)-alpha;
    else
        v(1)=-sigma/(x(1)+alpha);
    end
    beta=2*v(1)^2/(sigma+v(1)^2);
    v=v/v(1);
end
H=eye(n)-beta*v*v';
% e=H*x