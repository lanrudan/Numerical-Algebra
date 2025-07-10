% 幂法
function [maxlambda,max_eig]=power_method(A,ep,N)
% ep=1e-4;N=100;
n=length(A);
u=ones(n,1);
k=0;
m1=0;
while k<=N
   y=A*u;
   ab=abs(y);
   [i,~]=find(ab==max(ab));
   u=y/y(i(1,1),1);
   if norm(y(i(1,1),1)-m1,2)<ep
        break;
   end          
   m1=y(i(1,1),1);
   k=k+1;
end
maxlambda=y(i(1,1),1);
max_eig=y;
% norm(A*max_eig-maxlambda*max_eig,2);

