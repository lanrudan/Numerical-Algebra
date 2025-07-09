% 第31页 计算LDLT分解：改进的平方根法
function [error,x]=ldlt_decomposition(A,b)
% clc
% clear
% n=500;
% A=rand(n);
% A=A'*A;
A0=A;
[n,n]=size(A);
% b=rand(n,1);
D=zeros(n);
L=eye(n);
D(1,1)=A(1,1);
v=zeros(n,1);
for j=1:n
    for i=1:j-1
        v(i)=L(j,i)*D(i,i);
    end 
    D(j,j)=A(j,j)-L(j,1:j-1)*v(1:j-1);  
    L(j+1:n,j)=(A(j+1:n,j)-L(j+1:n,1:j-1)*v(1:j-1))./D(j,j);
end
error3=norm(L*D*L'-A0,2);
[error1,y]=fwd_sub(L,b);
[error2,x]=back_sub(D*L',y);
error=norm(A0*x-b);
