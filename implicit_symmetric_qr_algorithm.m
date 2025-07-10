% 209 隐式对称QR算法
function T=implicit_symmetric_qr_algorithm(A)
% clear;clc;
% A=rand(15);
% A=A'*A;
n=length(A);
[Q,T]=symmetric_tridiagonalization_householder(A);
eig(A);
u=10^-8;
while (1)
    d=abs(diag(H));%提取H的对角线元素绝对值
    d2=abs(diag(H,-1));%提取H的次对角线元素绝对值
    d3=(d(1:n-1)+d(2:n))*u;%计算用于比较的元素
    pl=find((d2-d3)<=0);
    for i=1:length(pl)
        H(pl(i)+1,pl(i))=0;
    end
    n=length(H);
    m=0;i=n;
    while (2)
        %考虑特殊情况
        if i==2
            m=m+2;
            break;
        end
        if i==1
            m=m+1;
            break;
        end
        %查找连续次对角线元素不为0
        %如果次对角线元素为0，利用continue命令跳出此次循环，进入下次循环
        if H(i,i-1)==0
            m=m+1;
            i=i-1;
            continue;
        end
        %只有当前一个次对角线元素不为0时，才可能进行这一步
        if H(i-1,i-2)==0
            m=m+2;
            i=i-2;
            continue;
        end
%         if H(i,i-1)==0 && H(i-1,i-2)==0
%             m=m+2;i=i-2;
%             break;
%         else
%             i=i-1;m=m+1;
%         end
%         continue;
        break;
    end
    %下面确定l的值
    %count 记录H22的大小
    %不确定l的值后面H22进行[H22,~]=doubleQR(H22)后会变成NAN l=n-m-count+1
    if m==n
        return;
    end
    count=2;
    for i=n-m-1:-1:2
        if H(i,i-1)==0
            break;
        else
            count=count+1;
        end
    end
    if m==n
        break;
    end
    H22=H(n-m-count+1:n-m,n-m-count+1:n-m);
    [H22,~]=double_shift_qr_iteration(H22);
    H(n-m-count+1:n-m,n-m-count+1:n-m)=H22;
end