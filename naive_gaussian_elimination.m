% 第18页 Gauss消去法
function [error,x]=naive_gaussian_elimination(A,b)
% n=500;
% % % A=[2 1 -1;-3 -1 2;-2 1 2];
% % % A=[2 1 2;5 -1 1;1 -3 -4];
% A=rand(n);
[n,n]=size(A);
A0=A;
L=eye(n);
for k=1:n-1
    L(k+1:n,k)=A(k+1:n,k)/A(k,k);
    A(k+1:n,k+1:n)=A(k+1:n,k+1:n)-L(k+1:n,k)*A(k,k+1:n);
end
%下面对A对角线以下的元素赋0
for k=1:n-1
    A(k+1:n,k)=zeros(n-k,1);
end
% error=norm(A0-L*A,2);%L为下三角矩阵，A变为上三角矩阵
% b=rand(n,1);
[error1,y]=fwd_sub(L,b);
[error2,x]=back_sub(A,y);
error=norm(A0*x-b,2);
% a=L*y-b;
% b=A0*x-y;
