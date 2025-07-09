% 第12页 前代法解下三角方程组
function [error,x]=fwd_sub(A,b)
% n=30;
% A=[1.2 0 0;3.2 -1.5 0;3 -2.7 4.1];
% b=[1.32;4.72;9.63];
% A=rand(n);
% A=tril(A);
n=length(A);
% b=rand(n,1);
bb=b;%储存bb=b,因为后面b改变了
x=b;
for j=1:n-1
    x(j)=b(j)/A(j,j);
    b(j+1:n)=b(j+1:n)-x(j)*A(j+1:n,j);
end
x(n)=b(n)/A(n,n);
error=norm(A*x-bb,2);
    