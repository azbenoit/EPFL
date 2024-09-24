clear all;close all;clc;
% Vars
L = 1;%m
k = 400;
Ta = 20;
Tb = 40;
Sc = 4;
Sl = -5;
n = 21;
dx = L/(n-1);
Dx = dx;
tol = 1e-3;
residual = 10;
residual_N = 10;
resids_P = [];
resids_N = [];

x = linspace(0,L,n);
T_guess = ones(1,n);
T_guess_N = ones(1,n);
A = zeros(n,n);

b = Sc*Dx*ones(n,1);
b(1) = Ta;
b(n) = Tb;
for i=2:n-1
    A(i,i) = 2*k/dx - Sl*Dx*T_guess(i)^2;
    A(i,i+1) = -k/dx;
    A(i, i-1) = -k/dx;
end
A(1,1) = 1;
A(n,n) = 1;

% Piccard
while residual >= tol
    for i=2:n-1
        A(i,i) = 2*k/dx - Sl*Dx*T_guess(i)^2;
    end
    T = A\b;
    T = reshape(T,1,n);
    residual = abs(T-T_guess);
    T_guess = T;
    % resids_P(end+1) = residual;
end

% Newton
% Sc = 4+10*T_guess_N;
% Sl = -15;
% while residual_N >= tol
%     b = Sc*Dx*ones(n,1);
%     b(1) = Ta;
%     b(n) = Tb;
%     for i=2:n-1
%         A(i,i) = 2*k/dx - Sl*Dx*T_guess_N(i)^2;
%     end
%     A(1,1) = 1;
%     A(n,n) = 1;
%     TN = A\b;
%     TN = reshape(TN,1,n);
%     residual = abs(TN-TN_guess);
%     T_guess_N = TN;
%     % resids_N(end+1) = residual_N;
% end

figure(1)
plot(x,T, 'b-x', 'DisplayName','T');