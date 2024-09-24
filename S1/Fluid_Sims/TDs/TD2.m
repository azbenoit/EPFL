clear all;close all;clc;
% Vars
L = 1;%m
k = 400;
Ta = 300;
Tb = 320;
S = 5000;
n = 21;
dx = L/(n-1);

x = linspace(0,L,n);
% T = zeros(1,n);
A = zeros(n,n);


b = S*dx*ones(n,1);
b(1) = Ta;
b(n) = Tb;
for i=2:n-1
    A(i,i) = 2*k/dx;
    A(i,i+1) = -k/dx;
    A(i, i-1) = -k/dx;
end
A(1,1) = 1;
A(n,n) = 1;

T = A\b;
T = reshape(T,1,n);

figure(1)
plot(x,T, 'b-x', 'DisplayName','T');