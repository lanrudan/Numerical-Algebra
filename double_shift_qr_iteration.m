% 第193页 双重歩位移的QR迭代
function [H,Q1]=double_shift_qr_iteration(A)
% clear;clc;
% n=10;A=rand(n);
n=length(A);
[Q,H]=hessenberg_decomposition_householder(A);
H0=H;
% norm(Q'*A*Q-H,2)
m=n-1;
s=H(m,m)+H(n,n);
t=H(m,m)*H(n,n)-H(m,n)*H(n,m);
x=H(1,1)*H(1,1)+H(1,2)*H(2,1)-s*H(1,1)+t;
y=H(2,1)*(H(1,1)+H(2,2)-s);
z=H(2,1)*H(3,2);
Q1=eye(n);Q2=eye(n);
for k=0:n-3
    [v,beta]=householder_transform([x,y,z]');
    q=max(1,k);
    I=eye(3);
    H(k+1:k+3,q:n)=(I-beta*v*v')*H(k+1:k+3,q:n);
    r=min(k+4,n);
    H(1:r,k+1:k+3)=H(1:r,k+1:k+3)*(I-beta*v*v');
    x=H(k+2,k+1);y=H(k+3,k+1);
    if k<n-3
        z=H(k+4,k+1);
    end
    Htilde=(eye(3)-beta*v*v');
    P=[eye(k,k),      zeros(k,3),        zeros(k,n-k-3);
       zeros(3,k),    Htilde,zeros(3,n-k-3);
       zeros(n-k-3,k),zeros(n-k-3,3),    eye(n-k-3,n-k-3)];
    Q1=Q1*P;
%     Q2=P*Q2;
end
[v,beta]=householder_transform([x,y]');
I=eye(2);
H(n-1:n,n-2:n)=(I-beta*v*v')*(H(n-1:n,n-2:n));
H(1:n,n-1:n)=H(1:n,n-1:n)*(I-beta*v*v');
for k=1:n-2
    H(k+2:end,k)=0;
end
Htilde=(eye(2)-beta*v*v');
P=[eye(n-2,n-2),zeros(n-2,2);
    zeros(2,n-2),Htilde];
Q1=Q1*P;
% Q2=P*Q2;
% norm(Q2*H0*Q1-H,2);
