% 第90页 Givens变换
function [c,s]=givens_rotation(a,b)
% clear;clc;
% x=rand(10,1);
% a=x(5);b=x(7);
if b==0
    c=1;
    s=0;
else
    if abs(b)>abs(a)
        tau=a/b;
        s=1/sqrt(1+tau^2);
        c=s*tau;
    else 
        tau=b/a;
        c=1/sqrt(1+tau^2);
        s=c*tau;
    end
end
y1=c*a+s*b;
y2=-s*a+c*b;
a=y1;b=y2;