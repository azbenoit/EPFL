clear all
close all
%clc

set(0,'DefaultLineLineWidth',2);
set(0,'DefaultAxesFontSize',18);
set(0,'DefaultLineMarkerSize',10);

%%%%%%%%%%%%% Physical parameters
L = 1;
Ta = 300;
Tb = 320;
k = 400;
Sc = 5000;   Sl = -100;

%%%%%%%%%%%%% Numerical parameters
n = 9;

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
A(n,n) = 1;
b(n)   = Tb;

%%%%%%%%%%%%%% Numerical solution  
%T=inv(A)*b;
T = A\b;

%%%%%%%%%%%%%% Solution plot
figure('color','w')
plot(x0,T,'bo')
hold on

%%%%%%%%%%%%%% Theoretical solution
if (Sl==0)
    alpha = (Tb-Ta)/L + Sc*L/(2*k);
    T_theo = -Sc/(2*k)*x0.^2 + alpha*x0 + Ta;
else
    mu1 =  sqrt(abs(Sl)/k);
    mu2 = -sqrt(abs(Sl)/k);
    c1  = (Tb-(Sc/Sl+Ta)*exp(mu2*L)+Sc/Sl)/(exp(mu1*L)-exp(mu2*L));
    c2  = Ta+Sc/Sl-c1;
    T_theo = c1*exp(mu1*x0)+c2*exp(mu2*x0)-Sc/Sl;
end

plot(x0,T_theo,'r-')

err = sum(abs( T-T_theo(:) ))/n

%%%%%%%%%%%%%% Finalize the figure
grid on, box on,
xlabel('x [m]'), ylabel('T [K]')
title(['Dirichlet T_a=' num2str(Ta), ', T_b=' num2str(Tb)...
    ', k=' num2str(k) ', S_c=' num2str(Sc) ', S_l=' num2str(Sl)...
    ', n=' num2str(n)])
legend('Numerical','Theoretical', 'location','best')
saveas(gcf,['Dirichlet-Ta' num2str(Ta), '_Tb' num2str(Tb) '_k' num2str(k)...
    '_Sc' num2str(Sc) '_Sl' num2str(Sl) '_n_' num2str(n) '.png'])
