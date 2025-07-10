% 第150页 解对称正定方程组：实用共轭梯度法
function [x,k]=conjugate_gradient_solver_practical(A,b,error,kmax)
% clear;clc;
% n=10;kmax=100;
% A=randn(n);
% A=A*A';
% b=rand(n,1);
n=length(A);
x0=ones(n,1);
x=x0;
k=0;r=b-A*x;rho=r'*r;
while (sqrt(rho)>error*norm(b,2))&&(k<kmax)
    k=k+1;
    if k==1
        p=r;
    else
        beta=rho/rho_;
        p=r+beta*p;
    end
    w=A*p;alpha=rho/(p'*w);
    x=x+alpha*p;
    r=r-alpha*w;rho_=rho;rho=r'*r;
end
% norm(A*x-b,2);