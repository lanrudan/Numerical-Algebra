%  第26页 列主元Gauss消去法
function [error,x]=gaussian_elimination_pivoting(A,b)
% n=3;
% A=[2 1 2;5 -1 1;1 -3 -4];
% % % % A=[2 1 -1;-3 -1 2;-2 1 2];
% A=rand(n);
[n,n]=size(A);
A0=A;
L=eye(n);
P=eye(n);
for k=1:n-1
    S=abs(A(k:n,k:n));%取出剩余矩阵
    [p,l]=find(S==max(S));
    A([k,p(1)+k-1],:)=A([p(1)+k-1,k],:);%交换第k行和第p(1)+k-1行(只要其中第一个元素)
    P([k,p(1)+k-1],:)=P([p(1)+k-1,k],:);
    u(k)=p(1);%记录置换矩阵Pk
    if A(k,k)~=0
        if k==1
            L(k+1:n,k)=A(k+1:n,k)/A(k,k);
        else
            L(k+1:n,k)=A(k+1:n,k)/A(k,k);
            L([k,p(1)+k-1],1:k-1)=L([p(1)+k-1,k],1:k-1);
        end
        A(k+1:n,k+1:n)=A(k+1:n,k+1:n)-L(k+1:n,k)*A(k,k+1:n);
    else
        stop
    end
end
for k=1:n-1
    A(k+1:n,k)=zeros(n-k,1);
end 
error3=norm(P*A0-L*A,2);
% b=rand(n,1);
[error1,y]=fwd_sub(L,P*b);
[error2,x]=back_sub(A,y);
error=norm(A0*x-b,2);


