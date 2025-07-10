% 第95页 计算QR分解：Householder方法
% 此程序已直接求解方程Ax=b
function [x,R]=qr_decomposition_householder(A,b)
% clear;clc;
[m,n]=size(A);
% m=500;n=300;
% A=rand(m,n);
% b=rand(m,1);
A0=A;
Q=eye(m);
for j=1:n
    if j<m
        [v,beta,H]=householder_transform(A(j:m,j));
        A(j:m,j:n)=A(j:m,j:n)-beta*v*v'*A(j:m,j:n);%R储存在A
        d(j)=beta;
        V(j:m-1,j)=v(2:m-j+1);
        H=[[eye(j-1),zeros(j-1,m-j+1)];
            [zeros(m-j+1,j-1),H]];
        Q=Q*H;
    end
end
R=triu(A(1:n,:));
Q1=Q(:,1:n);
c1=Q1'*b;
% norm(A0-Q1*R,2);
[error1,x]=back_sub(R,c1);
% error=norm(A0*x-Q1*R*x,2)