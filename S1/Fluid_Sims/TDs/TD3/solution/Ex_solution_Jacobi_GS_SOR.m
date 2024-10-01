clear all
close all

set(0,'DefaultLineLineWidth',2);
set(0,'DefaultAxesFontSize',14);
set(0,'DefaultLineMarkerSize',10);

%%%%%%%%%%%%% Physical parameters
L = 1;
Ta = 300;
Tb = 320;
k = 400;
Sc = 5000;   Sl = 0 %-100;

%%%%%%%%%%%%% Numerical parameters
nvec = [10 20 50 100 200];
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
    disp('inv(A)*b'), tic,  T1            = inv(A)*b;                          t=toc, t1(in) = t;
    disp('A\b'),      tic,  T2            = A\b;                               t=toc, t2(in) = t;
    disp('Jacobi'),   tic, [T4,iter4(in)] = my_Jacobi(A,b, Ta*ones(n,1));      t=toc, t4(in) = t;
    disp('GS'),       tic, [T5,iter5(in)] = my_GS(    A,b, Ta*ones(n,1));      t=toc, t5(in) = t;
    disp('SOR'),      tic, [T6,iter6(in)] = my_SOR(   A,b, Ta*ones(n,1), 1.8); t=toc, t6(in) = t;
    
    %%%%%%%%%%%%%% Solution plot
    figure('color','w'), hold on, grid on, box on,
    plot(x0,T1,'bo')
    plot(x0,T2,'kx')
    plot(x0,T4,'r+')
    plot(x0,T5,'gd')
    plot(x0,T6,'m*')
    
    %%%%%%%%%%%%%% Theoretical solution
    mu1 =  sqrt(abs(Sl)/k);
    mu2 = -sqrt(abs(Sl)/k);
    c1  = (Tb-(Sc/Sl+Ta)*exp(mu2*L)+Sc/Sl)/(exp(mu1*L)-exp(mu2*L));
    c2  = Ta+Sc/Sl-c1;
    Ttheo = c1*exp(mu1*x0)+c2*exp(mu2*x0)-Sc/Sl;
    
    plot(x0,Ttheo,'r-')
    
    %%%%%%%%%%%%%% Finalize the figure    
    xlabel('x [m]'), ylabel('T [K]')
    title(['T_a=' num2str(Ta), ', T_b=' num2str(Tb)...
        ', k=' num2str(k) ', S_c=' num2str(Sc) ', S_l=' num2str(Sl)...
        ', n=' num2str(n)])
    legend('inv(A)*b','A\b','Jacobi','GS','SOR','Theoretical', 'location','best')
    
    %%%%%%%%%%%%%% Error and residual
    err1(in) = sum(abs(T1 - Ttheo))/n;
    err2(in) = sum(abs(T2 - Ttheo))/n;
    err4(in) = sum(abs(T4 - Ttheo))/n;
    err5(in) = sum(abs(T5 - Ttheo))/n;
    err6(in) = sum(abs(T6 - Ttheo))/n;
        
    res1(in) = sum(abs(A*T1 - b)) / sum(abs( diag(A).*T1 ));
    res2(in) = sum(abs(A*T2 - b)) / sum(abs( diag(A).*T2 ));
    res4(in) = sum(abs(A*T4 - b)) / sum(abs( diag(A).*T4 ));         
    res5(in) = sum(abs(A*T5 - b)) / sum(abs( diag(A).*T5 ));         
    res6(in) = sum(abs(A*T6 - b)) / sum(abs( diag(A).*T6 ));         
end

figure, hold on, box on, grid on,
xlabel('n'), ylabel('calculation time')
plot(nvec, t1, 'bo')
plot(nvec, t2, 'kx')
plot(nvec, t4, 'r+') 
plot(nvec, t5, 'gd')
plot(nvec, t6, 'm*')
set(gca,'xscale','log'), set(gca,'yscale','log')
legend('inv(A)*b','A\b','Jacobi','GS','SOR', 'location','best')

figure, hold on, box on, grid on,
xlabel('n'), ylabel('number of iterations')
plot(nvec, iter4, 'r+') 
plot(nvec, iter5, 'gd')
plot(nvec, iter6, 'm*')
set(gca,'xscale','log'), set(gca,'yscale','log')
legend('Jacobi','GS','SOR', 'location','best')

figure, hold on, box on, grid on,
xlabel('n'), ylabel('error')
plot(nvec, err1, 'bo')
plot(nvec, err2, 'kx')
plot(nvec, err4, 'r+')
plot(nvec, err5, 'gd')
plot(nvec, err6, 'm*')
set(gca,'xscale','log'), set(gca,'yscale','log')
legend('inv(A)*b','A\b','Jacobi','GS','SOR', 'location','best')

figure, hold on, box on, grid on,
xlabel('n'), ylabel('scaled residual')
plot(nvec, res1, 'bo')
plot(nvec, res2, 'kx')
plot(nvec, res4, 'r+')
plot(nvec, res5, 'gd')
plot(nvec, res6, 'm*')
set(gca,'xscale','log'), set(gca,'yscale','log')
legend('inv(A)*b','A\b','Jacobi','GS','SOR', 'location','best')
