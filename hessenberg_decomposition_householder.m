% 第181页 计算上Hessenberg分解：Householder变换法
function [Q,A]=hessenberg_decomposition_householder(A)
% clear;clc;
% n=100;
% A=rand(n);A0=A;
[n,n]=size(A);
Q=eye(n);
for k=1:n-2
    [v,beta,~]=householder_transform(A(k+1:n,k));
    A(k+1:n,k:n)=(eye(n-k)-beta*v*v')*A(k+1:n,k:n);
    A(1:n,k+1:n)=A(1:n,k+1:n)*(eye(n-k)-beta*v*v');
    Htilde=(eye(n-k)-beta*v*v');
    H=[eye(k),zeros(k,n-k);
        zeros(n-k,k),Htilde];
    Q=Q*H;
end
for k=1:n-2
    A(k+2:end,k)=0;
end
% norm(A0*Q-Q*A,2);

