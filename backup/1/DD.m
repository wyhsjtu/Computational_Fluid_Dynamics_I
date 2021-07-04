function A = DD(n,h)
% DD(n,h)
%
% One-dimensional finite-difference derivative matrix
% of size n times n for second derivative:
%
% This function belongs to SG2212.m

A = zeros(n);  
A(1,1)=-1;
A(1,2)=1;
A(n,n-1)=1;
A(n,n)=-1;

for i=2:n-1
    A(i,i-1)=1;
    A(i,i)=-2;
    A(i,i+1)=1;
end
A=A/h^2;