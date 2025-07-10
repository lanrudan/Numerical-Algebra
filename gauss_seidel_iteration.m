% G-S迭代法
function [i,x]=gauss_seidel_iteration(error,A,b,max)
% clear;clc;
% n=100;
% n=size(A);
% for m=n
%     A=zeros(m,m);
%     for m=1:m
%         A(m,m)=20;
%     end
%     for m=2:m
%         A(m,m-1)=-8;
%         A(m-1,m)=-8;
%     end
%     for m=3:m
%         A(m,m-2)=1;
%         A(m-2,m)=1;
%     end
% end
% max=200;%最大迭代次数
% b=rand(n,1);
% error=1e-5;%误差限
[n,n]=size(A);
x0=ones(n,1);%初始迭代向量
D=diag(diag(A));
L=-(tril(A)-D);
U=-(triu(A)-D);
x=x0;
i=1;
while i<=max
    xafter=(D-L)\(U*x+b);
    if norm(xafter-x,2)<error
        break;
    end
    x=xafter;
    i=i+1;
end
x=xafter;
% norm(A*x-b,2);