% 第141页 解对称正定方程组最速下降法
function x=steepest_descent_method(A,b,error)
% clear;clc;
% n=100;
n=length(A);
% A=randn(n);
% A=A*A';
% b=rand(n,1);
x0=ones(n,1);
r0=b-A*x0;k=0;
x=x0;r=r0;
while norm(r,2)>error 
    k=k+1;
    alpha=r'*r/(r'*(A*r));
    x=x+alpha*r;
    r=b-A*x;
end
% norm(A*x-b,2);
