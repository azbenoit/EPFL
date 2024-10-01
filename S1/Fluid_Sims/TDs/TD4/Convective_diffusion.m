clear all
close all
%clc

set(0,'DefaultLineLineWidth',2);
set(0,'DefaultAxesFontSize',18);
set(0,'DefaultLineMarkerSize',10);

%%%%%%%%%%%%% Physical parameters
L = 1;
gamma = 0.1;
rho = 1;
u = 2;
phi_0 = 1;
phi_L = 0;


%%%%%%%%%%%%% Numerical parameters
n = 20;


%%%%%%%%%%%%%% Grid generation
x0=linspace(0,L,n);

dx=L/(n-1);
Dx=dx;

F = rho*u;
D = gamma /dx;

%%%%%%%%%%%%%% Creating the matrix
A = zeros(n,n);
b = zeros(n,1);

for i=2:n-1 
    A(i,i-1) = D+F/2;
    A(i,i+1) = D-F/2;
    A(i,i)   = 2*D ;
    % b(i)     = 0;
end

% Boundary conditions
A(1,1) = 1;
b(1)   = phi_0;
A(n,n) = 1;
b(n)   = phi_L;

%%%%%%%%%%%%%% Numerical solution  
%T=inv(A)*b;
phi = A\b;

%%%%%%%%%%%%%% Solution plot
figure('color','w')
plot(x0,phi,'bo')
hold on

%%%%%%%%%%%%%% Theoretical solution
% if (Sl==0)
%     alpha = (Tb-Ta)/L + Sc*L/(2*k);
%     T_theo = -Sc/(2*k)*x0.^2 + alpha*x0 + Ta;
% else
%     mu1 =  sqrt(abs(Sl)/k);
%     mu2 = -sqrt(abs(Sl)/k);
%     c1  = (Tb-(Sc/Sl+Ta)*exp(mu2*L)+Sc/Sl)/(exp(mu1*L)-exp(mu2*L));
%     c2  = Ta+Sc/Sl-c1;
%     T_theo = c1*exp(mu1*x0)+c2*exp(mu2*x0)-Sc/Sl;
% end
Pe = (rho*u*L)/gamma;
phi_theo = phi_0 + (phi_L-phi_0)*(exp(Pe*x0/L)-1)/(exp(Pe)-1);

plot(x0,phi_theo,'r-')

err = sum(abs( phi-phi_theo(:) ))/n

%%%%%%%%%%%%%% Finalize the figure
grid on, box on,
xlabel('x [m]'), ylabel('phi')
title(['Dirichlet T_a=' num2str(phi_0), ', T_b=' num2str(phi_L)...
    ', gamma=' num2str(gamma) ', n=' num2str(n)])
legend('Numerical','Theoretical', 'location','best')
% saveas(gcf,['Dirichlet-Ta' num2str(Ta), '_Tb' num2str(Tb) '_k' num2str(k)...
    % '_Sc' num2str(Sc) '_Sl' num2str(Sl) '_n_' num2str(n) '.png'])
