% 207页 计算三对角分解：Householder变换法
function [Q,T]=symmetric_tridiagonalization_householder(A)
% clear;clc;
% A=rand(10);
% A=A'*A;
% A0=A;
n=length(A);
Q=eye(n);
for k=1:n-2
    [v,beta,H]=householder_transform(A(k+1:n,k));
    Hk=[eye(k),zeros(k,n-k);
        zeros(n-k,k),H];
    u=beta*A(k+1:n,k+1:n)*v;
    w=u-(beta*u'*v/2)*v;
    A(k+1,k)=norm(A(k+1:n,k),2);
    A(k,k+1)=A(k+1,k);
    A(k+1:n,k+1:n)=A(k+1:n,k+1:n)-v*w'-w*v';
    Q=Q*Hk;
end
% T=Q'*A0*Q;
T=diag(diag(A))+diag(diag(A,-1),-1)+diag(diag(A,1),1);