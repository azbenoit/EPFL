function Ex_nonlinearity

clear all
close all

set(0,'DefaultLineLineWidth',2);
set(0,'DefaultAxesFontSize',14);
set(0,'DefaultLineMarkerSize',10);

%%%%%%%%%%%%% Physical parameters
L = 1;
Ta = 20;
Tb = 40;
k = 400;

%%%%%%%%%%%%% Numerical parameters
n = 50;

%lin_method = 'Picard';   % choose Picard or Newton
 lin_method = 'Newton';   % choose Picard or Newton

tol = 1e-12;   % Tolerance
omega = 1;     % Relaxation coefficient

 T_guess = 30 * ones(n,1);   % initial guess
%T_guess =  0 * ones(n,1);   % initial guess

%%%%%%%%%%%%%% Grid generation
x0 = linspace(0,L,n);
x0 = x0(:);   % Make sure it's a column vector

dx=L/(n-1);
Dx=dx;

%%%%%%%%%%%%%%
fig_T  = figure('color','w'); hold on, box on, grid on, 
ylim([5 45]), xlabel('x [m]'), ylabel('T [K]'), title('iteration 0')
plot(x0,T_guess,'o-'), 

dT_rel = 1;
dT_rel_vec = [];
iter = 0;
while dT_rel>tol   % Loop until increment varies less than tolerance
    iter = iter + 1;
        
    %%%%%%%%%%%%%% Creating the matrix
    A = zeros(n,n);
    b = zeros(n,1);
            
    %--- Linearize the source term S = 4-5*T^3
    switch lin_method
        case 'Picard',   Sc = 4*ones(n,1);        Sl =  -5*T_guess.^2;
        case 'Newton',   Sc = 4 +10*T_guess.^3;   Sl = -15*T_guess.^2;
    end
    
    for i=2:n-1
        A(i,i-1) = -k/dx;
        A(i,i+1) = -k/dx;
        A(i,i)   = 2*k/dx -Sl(i)*Dx;
        b(i)     = Sc(i)*Dx;
    end
    
    %--- Boundary conditions
    A(1,1) = 1;
    b(1)   = Ta;
    A(n,n) = 1;
    b(n)   = Tb;
    
    %%%%%%%%%%%%%% Numerical solution
    T = A\b;
    
    dT_rel = sum(abs( T-T_guess )) / sum(abs( T ));
    dT_rel_vec(iter) = dT_rel;

    figure(fig_T), hold off, plot(x0,T,'o-'), hold on, grid on,
    xlabel('x [m]'), ylabel('T [K]'), title(['iteration ' num2str(iter)])
    drawnow 
    
    %--- Update the current solution, possibly with relaxation
    T_guess = omega*T + (1-omega)*T_guess;
end

res = sum(abs(A*T - b)) / sum(abs( diag(A).*T )),   % Final normalized residual

figure('color','w'); hold on, box on, grid on, 
xlabel('iter'), ylabel('relative ||dT||')
plot(1:iter, dT_rel_vec,'sk-'), set(gca,'yscale','log')
ylim([1e-15 1])

%%%%%%%%%%%%%% Compare with Matlab's bvp4c
% Rewrite the equation kT''+S(T)=0 as a 1st-order system 
% of 2 variables T1=T and T2=T':
% T1' = T2 
% T2' = T1'' = -S(T)/k
function dTdx = two_ode(x,T)   
    dTdx = [T(2);
            -(4-5*T(1)^3)/k];  
end
function res = two_bc(Tleft,Tright)
    res = [Tleft(1)-Ta; 
           Tright(1)-Tb];
end
solinit = bvpinit(x0, [30, 0]);
sol = bvp4c(@two_ode, @two_bc, solinit);
T_bvp4c = deval(sol,x0);

figure(fig_T), 
title(['T_a=' num2str(Ta), ', T_b=' num2str(Tb) ', k=' num2str(k) ', n=' num2str(n)])
plot(x0, T_bvp4c(1,:), 'r--')
legend('This FVM code','Matlab bvp4c', 'location','best')

diff_FVM_bvp4c = sum(abs( T-T_bvp4c(1,:)' ))/n,

end