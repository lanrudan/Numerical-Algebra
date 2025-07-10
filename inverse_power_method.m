% 169 反幂法
function [minlambda,min_eig]=inverse_power_method(A,ep,N)
n=length(A);
z=ones(n,1);
k=0;
% ep=1e-4;N=100;
m1=0;
while k<=N
   y=gaussian_elimination_pivoting(A,z);
   ab=abs(y);
   [i,~]=find(ab==max(ab));
   z=y/y(i(1,1),1);
   if norm(y(i(1,1),1)-m1,2)<ep
        break;
   end          
   m1=y(i(1,1),1);
   k=k+1;
end
minlambda=y(i(1,1),1);
min_eig=y;
% norm(A*max_eig-maxlambda*max_eig,2);