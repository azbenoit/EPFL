clear all
close all

set(0,'DefaultLineLineWidth',2);
set(0,'DefaultAxesFontSize',14);
set(0,'DefaultLineMarkerSize',10);

%%%%%%%%%%%%% Physical parameters
L = 1;
Ta = 300;
Tb = 320;
lambda = 400;
Sc = 5000;   Sp = -100;

%%%%%%%%%%%%% Numerical parameters
nvec = [10 20 50 100 500 1000 5000 1e4];
nn = length(nvec);

for in=1:nn
    disp('---------------'),
    n = nvec(in),
    
    %%%%%%%%%%%%%% Grid generation
    x0 = linspace(0,L,n);
    x0 = x0(:); % Make sure it's a column vector
    
    dx=L/(n-1);
    Dx=dx;
    
    %%%%%%%%%%%%%% Creating the matrix
    A = zeros(n,n);
    b = zeros(n,1);
    
    for i=2:n-1
        A(i,i-1) = -lambda/dx;
        A(i,i+1) = -lambda/dx;
        A(i,i)   = 2*lambda/dx-Sp*Dx;
        b(i)     = Sc*Dx;
    end
    
    % Boundary conditions
    A(1,1) = 1;
    b(1)   = Ta;
    A(n,n) = 1;
    b(n)   = Tb;
    
    %%%%%%%%%%%%%% Numerical solution
    disp('inv(A)*b'), tic, T1 = inv(A)*b;       t=toc, t1(in) = t;
    disp('A\b'),      tic, T2 = A\b;            t=toc, t2(in) = t;
    disp('TDMA'),     tic, T3 = my_TDMA(  A,b); t=toc, t3(in) = t;
    
    %%%%%%%%%%%%%% Solution plot
    figure('color','w'), hold on, grid on, box on,
    plot(x0,T1,'bo')
    plot(x0,T2,'kx')
    plot(x0,T3,'ms')
    
    %%%%%%%%%%%%%% Theoretical solution
    mu1 =  sqrt(abs(Sp)/lambda);
    mu2 = -sqrt(abs(Sp)/lambda);
    c1  = (Tb-(Sc/Sp+Ta)*exp(mu2*L)+Sc/Sp)/(exp(mu1*L)-exp(mu2*L));
    c2  = Ta+Sc/Sp-c1;
    Ttheo = c1*exp(mu1*x0)+c2*exp(mu2*x0)-Sc/Sp;
    
    plot(x0,Ttheo,'r-')
    
    %%%%%%%%%%%%%% Finalize the figure    
    xlabel('x [m]'), ylabel('T [K]')
    title(['x_1=0, x_2=' num2str(L) ', T_a=' num2str(Ta), ', T_b=' num2str(Tb)...
        ', \lambda=' num2str(lambda) ', S_C=' num2str(Sc) ', S_P=' num2str(Sp)...
        ', N=' num2str(n)])
    legend('inv(A)*b','A\b','TDMA','Theoretical', 'location','best')
    
    %%%%%%%%%%%%%% Error and residual
    err1(in) = sum(abs(T1 - Ttheo))/n;
    err2(in) = sum(abs(T2 - Ttheo))/n;
    err3(in) = sum(abs(T3 - Ttheo))/n;
    
    %err1(in) = norm(T1 - Ttheo, 1)/n;
    %err2(in) = norm(T2 - Ttheo, 1)/n;
    %err3(in) = norm(T3 - Ttheo, 1)/n;
    
    %err1(in) = norm(T1 - Ttheo, 2)/n;
    %err2(in) = norm(T2 - Ttheo, 2)/n;
    %err3(in) = norm(T3 - Ttheo, 2)/n;
    
    res1(in) = sum(abs(A*T1 - b)) / sum(abs( diag(A).*T1 ));
    res2(in) = sum(abs(A*T2 - b)) / sum(abs( diag(A).*T2 ));
    res3(in) = sum(abs(A*T3 - b)) / sum(abs( diag(A).*T3 )); 
end

figure, hold on, box on, grid on,
xlabel('n'), ylabel('calculation time')
plot(nvec, t1, 'bo')
plot(nvec, t2, 'kx')
plot(nvec, t3, 'ms')    
set(gca,'xscale','log'), set(gca,'yscale','log')
legend('inv(A)*b','A\b','TDMA', 'location','best')

figure, hold on, box on, grid on,
xlabel('n'), ylabel('error')
plot(nvec, err1, 'bo')
plot(nvec, err2, 'kx')
plot(nvec, err3, 'ms')
set(gca,'xscale','log'), set(gca,'yscale','log')
legend('inv(A)*b','A\b','TDMA', 'location','best')

figure, hold on, box on, grid on,
xlabel('n'), ylabel('scaled residual')
plot(nvec, res1, 'bo')
plot(nvec, res2, 'kx')
plot(nvec, res3, 'ms')
set(gca,'xscale','log'), set(gca,'yscale','log')
legend('inv(A)*b','A\b','TDMA', 'location','best')
