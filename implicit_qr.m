% 194 隐式QR
function H=implicit_qr(A)
%A=rand(6);
n=length(A);
[~,H]=hessenberg_decomposition_householder(A);
% eig(A);
u=10^-6;
while (1)
    d=abs(diag(H));
    d2=abs(diag(H,-1));
    d3=(d(1:n-1)+d(2:n))*u;
    pl=find((d2-d3)<=0);
    for i=1:length(pl)
        H(pl(i)+1,pl(i))=0;
    end
    n=length(H);
    i=n;
    while (2)

        if i==2
            break;
        end
        if i==1
            break;
        end
       
        if H(i,i-1)==0
            i=i-1;
            continue;
        end
       
        if H(i-1,i-2)==0
            i=i-2;
            continue;
        end
        break;
    end
  
    if i==1
        break;
    end
    count=2;
    for  k=i-1:-1:2
        if H(k,k-1)==0
            break;
        else
            count=count+1;
        end
    end
    H22=H(i-count+1:i,i-count+1:i);
    [H22,~]=double_shift_qr_iteration(H22);
    H(i-count+1:i,i-count+1:i)=H22;
end