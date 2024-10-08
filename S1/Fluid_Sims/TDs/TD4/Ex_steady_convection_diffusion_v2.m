clear all
% close all
clc

set(0,'DefaultLineLineWidth',2);
fs=18;   set(0,'DefaultAxesFontSize',fs);
set(0,'DefaultLineMarkerSize',10);

%%%%%%%%%%%%% Physical parameters
L = 1;
Phi_a = 1;
Phi_b = 0;
u = 10;
rho = 1;
gamma = 0.1;

%%%%%%%%%%%%% Numerical parameters
n = 41;

%%%%%%%%%%%%%% Grid generation
dx = L/(n-1);
x0 = 0:dx:L;

%%%%%%%%%%%%%%
F = rho*u;
D = gamma/dx;

%%%%%%%%%%%%%% Theoretical solution
x_theo = linspace(0,L,10001);
Phi_theo = Phi_a + (Phi_b-Phi_a)*(exp(rho*u*x_theo/gamma)-1)/(exp(rho*u*L/gamma)-1);

figure('color','w'), hold on, grid on, box on
xlabel('x [m]'), ylabel('\phi','fontsize',fs+2)
plot(x_theo, Phi_theo, 'k-')
title(['L=1 m, \phi_0=1, \phi_L=0, \rho=1 kg/m^3, \Gamma=0.1 kg/(m.s), u=' num2str(u) ' m/s'])
set(gca,'ytick',0:0.2:2)
drawnow



%%%%%%%%%%%%%% 1. CD
A = zeros(n,n);
b = zeros(n,1);
for i=2:n-1
    A(i,i-1) = -(D+F/2);
    A(i,i+1) = -(D-F/2);
    A(i,i)   = -(A(i,i-1)+A(i,i+1));
end

%--- Boundary conditions in CVs 1 and n
A(1,1) = 1;
b(1)   = Phi_a;
A(n,n) = 1;
b(n)   = Phi_b;

%--- Solution
Phi = A\b;
plot(x0, Phi, 'ro--','MarkerFaceColor','r'), drawnow

x_theo=0:L/(n-1):L;
Phi_theo=Phi_a+(Phi_b-Phi_a)*(exp(rho*u*x_theo/gamma)-1)/(exp(rho*u*L/gamma)-1);
errCD = mean(abs(Phi(1:end-1)-Phi_theo(1:end-1)')./Phi_theo(1:end-1)')



%%%%%%%%%%%%% 2. UD
A = zeros(n,n);   
b = zeros(n,1);
for i=2:n-1    
    A(i,i-1) = -(D+F);
    A(i,i+1) = -(D);
    A(i,i)   = -(A(i,i-1)+A(i,i+1));
end

%--- Boundary conditions in CVs 1 and n
A(1,1) = 1;
b(1)   = Phi_a;
A(n,n) = 1;
b(n)   = Phi_b;

%--- Solution
Phi = A\b;
plot(x0, Phi, 'd--','color',[0 0.5 0],'MarkerFaceColor',[0 0.5 0]), drawnow

errUD = mean(abs(Phi(1:end-1)-Phi_theo(1:end-1)')./Phi_theo(1:end-1)')



%%%%%%%%%%%%% 3. QUICK + dc
tol = 1e-12;
maxIter = 1000;

A = zeros(n,n);
b = zeros(n,1);
Phi_new = ones(n,1);
nIter = 1;

dPhi_relative = 1;
while dPhi_relative>tol && nIter<maxIter
       
    Phi = Phi_new;    
    
    for i=3:n-1     
        A(i,i-1) = -(D+F);
        A(i,i+1) = -(D);
        A(i,i)   = -(A(i,i-1)+A(i,i+1));
        b(i)     = 1/8*F * (5*Phi(i)-3*Phi(i+1)-Phi(i-1)-Phi(i-2)); % deferred correction
    end
    
    %--- Boundary conditions in CVs 1 and n
    A(1,1) = 1;
    b(1)   = Phi_a;
    A(n,n) = 1;
    b(n)   = Phi_b;
    
    %--- Boundary condition in CV 2
    A(2,1) = -(D+F/2);      % West
    A(2,3) = -(D);          % East
    A(2,2) = 2*D + F - F/2; % P
    b(2) = F/8 *(Phi(1) +2*Phi(2) -3*Phi(3)); % deferred correction
        
    %--- Solution
    Phi_new = A\b;
    dPhi_relative = norm(abs(Phi_new-Phi))/norm(abs(Phi));    
    nIter = nIter+1;
end
Phi = Phi_new;
plot(x0, Phi, 'bs--','MarkerFaceColor','b'), drawnow

errQUICK = mean(abs(Phi(1:end-1)-Phi_theo(1:end-1)')./Phi_theo(1:end-1)')



%%%%%%%%%%%%% 4. TVD + dc (flux limiter: Van Leer psi=(r+|r|)/(1+r))
tol = 1e-12;
maxIter = 1000;

A = zeros(n,n);
b = zeros(n,1);
Phi_new = ones(n,1);
nIter = 1;

dPhi_relative = 1;
while dPhi_relative>tol && nIter<maxIter
       
    Phi = Phi_new;    
    
    for i=3:n-1        
        A(i,i-1) = -(D+F);
        A(i,i+1) = -(D);
        A(i,i)   = -(A(i,i-1)+A(i,i+1));
        
        %--- Define r (ratio of gradients) with test to avoid dividing by 0
        if abs(Phi(i+1)-Phi(i))<1e-12
            if abs(Phi(i)-Phi(i-1))>1e-12
                r_e = 1e30 * sign(Phi(i)-Phi(i-1)); 
            else
                r_e = 0;
            end
        else
            r_e = (Phi(i)-Phi(i-1)) / (Phi(i+1)-Phi(i));
        end
        if abs(Phi(i)-Phi(i-1))<1e-12
            if abs(Phi(i-1)-Phi(i-2))>1e-12
                r_w = 1e30 * sign(Phi(i-1)-Phi(i-2));
            else
                r_w = 0;
            end
        else
            r_w = (Phi(i-1)-Phi(i-2)) / (Phi(i)-Phi(i-1));
        end
        
        %--- Define psi
        psi_e = (r_e+abs(r_e)) / (1+r_e);
        psi_w = (r_w+abs(r_w)) / (1+r_w);
        b(i)  = -F/2 * (psi_e*(Phi(i+1)-Phi(i))-psi_w*(Phi(i)-Phi(i-1))); % deferred correction
    end
    
    %--- Boundary conditions in CVs 1 and n
    A(1,1) = 1;
    b(1)   = Phi_a;
    A(n,n) = 1;
    b(n)   = Phi_b;
    
    %--- Boundary condition in CV 2
    A(2,1) = -(D+F);           % West
    A(2,3) = -(D);             % East
    A(2,2) = -(A(2,1)+A(2,3)); % P
    if abs(Phi(3)-Phi(2))<1e-12
        if abs(Phi(2)-Phi(1))>1e-12
            r_e2 = 1e30 * sign(Phi(2)-Phi(1));
        else
            r_e2 = 0;
        end
    else
        r_e2 = (Phi(2)-Phi(1))/(Phi(3)-Phi(2));
    end
    r_w2 = 1;
    psi_e2 = (r_e2+abs(r_e2))/(1+r_e2);
    psi_w2 = (r_w2+abs(r_w2))/(1+r_w2);       
    b(2)   = -F/2 * (psi_e2*(Phi(3)-Phi(2)) -psi_w2*(Phi(2)-Phi(1))); % deferred correction
        
    %--- Solution
    Phi_new = A\b;
    dPhi_relative = norm(abs(Phi_new-Phi))/norm(abs(Phi));
    nIter = nIter+1;
end
Phi = Phi_new;
plot(x0, Phi, 'mv--','MarkerFaceColor','m')

errTVD = mean(abs(Phi(1:end-1)-Phi_theo(1:end-1)')./Phi_theo(1:end-1)')

legend('Exact solution','CD','UD','QUICK','TVD (Van Leer)', 'location','best')
