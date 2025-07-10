% 第24页 全主元高斯消去法
function [error,A]=gaussian_elimination_complete_pivoting(A)
% n=50;
% % A=[2 1 2;5 -1 1;1 -3 -4];
% % A=[2 1 -1;-3 -1 2;-2 1 2];
% % A=[2,1,-3,-1;3,1,0,7;-1,2,4,-2;1,0,-1,5]
% A=rand(n);
[n,n]=size(A);
A0=A;
L=eye(n);
P=eye(n);
Q=eye(n);
for k=1:n-1
    S=abs(A(k:n,k:n));%取出剩余矩阵，并取绝对值
    [p,q]=find(S==max(max(S)));%找出绝对值最大值的位置
    A([k,p+k-1],:)=A([p+k-1,k],:);%交换当时矩阵的第p行(p+k-1)和第k行
    A(:,[k,q+k-1])=A(:,[q+k-1,k]);%交换当时矩阵的第q列(q+k-1)和第k列
    P([k,p+k-1],:)=P([p+k-1,k],:);
    Q(:,[k,q+k-1])=Q(:,[q+k-1,k]);
    u(k)=p+k-1;%记录置换矩阵Pk
    v(k)=q+k-1;%记录置换矩阵Qk
    if A(k,k)~=0
        if k==1
        L(k+1:n,k)=A(k+1:n,k)/A(k,k);
        else  
        L(k+1:n,k)=A(k+1:n,k)/A(k,k);
        L([k,p+k-1],1:k-1)=L([p+k-1,k],1:k-1);
        end
        A(k+1:n,k+1:n)=A(k+1:n,k+1:n)-L(k+1:n,k)*A(k,k+1:n);
    else
        stop
    end
end
for k=1:n-1
    A(k+1:n,k)=zeros(n-k,1);
end
% P*A0*Q
% L*A
% P*A0*Q-L*A
error=norm(P*A0*Q-L*A,2);


    