% 第153页 解对称正定方程组：预优共轭梯度法
function [x,k]=preconditioned_conjugate_gradient(A,b,error,kmax)
% clear;clc;
% n=10;kmax=500;epsi=1e-5;
% A=randn(n);
% A=A*A';
% b=rand(n,1);
n=length(A);
x0=ones(n,1);
x=x0;
k=0;r=b-A*x;
M=diag(diag(A));
while (sqrt(r'*r)>error*norm(b,2))&&(k<kmax)
    z=M\r;
    k=k+1;
    if k==1
        p=z;rho=r'*z;
    else
        rho_=rho;rho=r'*z;
        beta=rho/rho_;
        p=z+beta*p;
    end
    w=A*p;alpha=rho/(p'*w);
    x=x+alpha*p;
    r=r-alpha*w;
end
% norm(A*x-b,2)