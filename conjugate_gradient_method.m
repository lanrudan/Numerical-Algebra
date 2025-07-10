% 第146页 解对称正定方程组：共轭梯度法
function [k,x]=conjugate_gradient_method(A,b,error)
% clear;clc;
% n=10;
% A=randn(n);
% A=A*A';
[n,n]=size(A);
x0=ones(n,1);
% b=rand(n,1);
x0=ones(n,1);
r0=b-A*x0;k=0;
x=x0;r=r0;
while norm(r,2)>error
    k=k+1;
    if k==1
        p=r;
        alpha=r'*r/(p'*(A*p));
        x=x+alpha*p;
        rafter=r-alpha*A*p;      
    else
        beta=rafter'*rafter/(r'*r);
        p=rafter+beta*p;
        alpha=rafter'*rafter/(p'*(A*p));
        x=x+alpha*p;
        r=rafter;
        rafter=rafter-alpha*A*p;
    end
end
% norm(A*x-b,2)
