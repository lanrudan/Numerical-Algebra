% Jacobi迭代法
function [i,x]=jacobi_iteration(error,A,b,max)
% clear;clc;
% n=100;
% n=size(A);
[n,n]=size(A);
x0=ones(n,1);%初始迭代向量
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
% max=100;%最大迭代次数
% error=1e-6;%误差限
% b=rand(n,1);
D=diag(diag(A));
L=-tril(A,-1);
U=-triu(A,1);
x=x0;
i=1;
while i<=max

    xafter=D\((L+U)*x+b);
    if norm(xafter-x,2)<error
        break;
    end
    x=xafter;
    i=i+1;
end
x=xafter;


