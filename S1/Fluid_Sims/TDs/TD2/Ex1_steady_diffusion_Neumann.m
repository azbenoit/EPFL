clear all
close all
%clc

set(0,'DefaultLineLineWidth',2);
set(0,'DefaultAxesFontSize',18);
set(0,'DefaultLineMarkerSize',10);

%%%%%%%%%%%%% Physical parameters
L = 1;
Ta = 300;
k = 400;
Sc = 5000;
Sl = 0;

qb = -5000;

%%%%%%%%%%%%% Numerical parameters
n = 21;

%%%%%%%%%%%%%% Grid generation
x0=linspace(0,L,n);

dx=L/(n-1);
Dx=dx;

%%%%%%%%%%%%%% Creating the matrix
A = zeros(n,n);
b = zeros(n,1);

for i=2:n-1 
    A(i,i-1) = -k/dx;
    A(i,i+1) = -k/dx;
    A(i,i)   = 2*k/dx -Sl*Dx;
    b(i)     = Sc*Dx;
end

% Boundary conditions
A(1,1) = 1;
b(1)   = Ta;

A(n,n)   = k/dx;
A(n,n-1) = -k/dx;
b(n)     = Sc*Dx/2 -qb;

%%%%%%%%%%%%%% Numerical solution  
%T=inv(A)*b;
T = A\b;

%%%%%%%%%%%%%% Solution plot
figure('color','w')
plot(x0,T,'bo')
hold on

%%%%%%%%%%%%%% Theoretical solution
if (Sl==0)
    alpha = (Sc*L-qb)/k;
    T_theo = -Sc/(2*k)*x0.^2 + alpha*x0 + Ta;
else    
end

plot(x0,T_theo,'r-')

%%%%%%%%%%%%%% Finalize the figure
grid on, box on,
xlabel('x [m]'), ylabel('T [K]')
title(['Neumann T_a=' num2str(Ta) ', q_b=' num2str(qb)...
    ', k=' num2str(k) ', S_c=' num2str(Sc) ', S_l=' num2str(Sl)...
    ', n=' num2str(n)])
legend('Numerical','Theoretical', 'location','best')
saveas(gcf,['Neumann-Ta' num2str(Ta)  '_k' num2str(k)...
    '_Sc' num2str(Sc) '_Sl' num2str(Sl) '_qb' num2str(qb) '_n_' num2str(n) '.png'])
