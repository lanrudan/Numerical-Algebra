% 第13页 回代法解上三角方程
function [error,x]=back_sub(A,b)
% n=50;
% % % A=[4 -1 2 3;0 -2 7 -4;0 0 6 5;0 0 0 3];
% % % b=[20;-7;4;6];
% A=rand(n);
% b=rand(n,1);
% A=triu(A);
n=length(A);
bb=b;
x=b;
for j=n:-1:2
    x(j)=b(j)/A(j,j);
    b(1:j-1)=b(1:j-1)-x(j)*A(1:j-1,j);
end
x(1)=b(1)/A(1,1);
error=norm(A*x-bb,2);
    