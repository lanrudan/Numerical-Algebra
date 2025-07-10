% 71页 估计矩阵的1范数：优化法
% 输出A的转置的逆的1范数，也即A的逆的无穷范数，便于求条件数
function v=hager_norm_estimator(A)
% clear;clc;
% n=5;
% A=hilb(n);
n=length(A);
x=1./n*ones(n,1);
k=1;
while k==1
    [error1,w]=gaussian_elimination_pivoting(A',x);
    v=sign(w);
    [error2,z]=gaussian_elimination_pivoting(A,v);
    if max(abs(z))<=z'*x
        v=norm(w,1);
        k=0;
    else
        j=find(abs(z)==max(abs(z)));
        x=[zeros(1,j(1)-1),1,zeros(1,n-j(1))]';
        k=1;
    end
end

