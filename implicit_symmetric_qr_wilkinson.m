% 209 带Wilkinson位移的隐式对称QR迭代
function T=implicit_symmetric_qr_wilkinson(A)
% clear;clc;
% A=rand(100);
% A=A'*A;
n=length(A);
[~,T]=symmetric_tridiagonalization_householder(A);
% T0=T;
d=(T(n-1,n-1)-T(n,n))/2;
mu=T(n,n)-T(n,n-1)^2/(d+sign(d)*sqrt(d^2+T(n,n-1)^2));
x=T(1,1)-mu;z=T(2,1);
for k=1:n-1
    [c,s]=givens_rotation(x,z);
    Gk=eye(n);
    Gk([k,k+1],[k,k+1])=[c,s;
                         -s,c];
    T=Gk*T*Gk';
    if k<n-1
        x=T(k+1,k);z=T(k+2,k);
    end
end
T=diag(diag(T))+diag(diag(T,-1),-1)+diag(diag(T,1),1);