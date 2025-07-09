% 第30页 计算Cholesky分解：平方根法
function [error,x]=cholesky_solve(A,b)
% clc
% clear
% n=800;
% A=randn(n);
% A=A'*A;
A0=A;
% b=rand(n,1);
[n,n]=size(A);
L=zeros(n);
for k=1:n
    L(k,k)=sqrt(A(k,k));
    L(k+1:n,k)=A(k+1:n,k)/L(k,k);
        for j=k+1:n
            A(j:n,j)=A(j:n,j)-L(j:n,k)*L(j,k);
        end
end
% error3=norm(A0-L*L',2);
[error1,y]=OKqiandaifa(L,b);
[error2,x]=OKhuidaifa(L',y);
error=norm(A0*x-b,2);

