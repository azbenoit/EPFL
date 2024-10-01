function Ex_solution_multigrid_simple_full

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
n = 17; 
tol = 1e-10;
T_guess = zeros(n,1);

%%%%%%%%%%%%%% Grid generation
x0 = linspace(0,L,n);
x0 = x0(:);   % Make sure it's a column vector

dx=L/(n-1);
Dx=dx;

%%%%%%%%%%%%%% Creating the original matrix and RHS
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

%%%%%%%%%%%%%% Interpolation operators
n1 = n,
n2 = (n1+1)/2,
n3 = (n2+1)/2,

TT1 = zeros(n,n2);
for i=1:2:n,   TT1(i,(i+1)/2) = 1; end
for i=2:2:n-1, TT1(i,i/2:i/2+1) = 1/2; end
TT1,

TT2 = zeros(n2,n3);
for i=1:2:n2,   TT2(i,(i+1)/2) = 1; end
for i=2:2:n2-1, TT2(i,i/2:i/2+1) = 1/2; end
TT2,

%%%%%%%%%%%%%% Restriction operators
RR1 = transpose(TT1)/2;
RR1(1,1:2) = [1,0];   
RR1(end,end-1:end) = [0,1],

RR2 = transpose(TT2)/2;
RR2(1,1:2) = [1,0];   
RR2(end,end-1:end) = [0,1],

%%%%%%%%%%%%%% System matrices 
A1 = A,
A2 = RR1*A1*TT1,
A3 = RR2*A2*TT2,


%%%%%%%%%%%%%% Main loop for V-cycles
figure; 
res = 1;
res_vec = [];
iter      = 0;
iter_fine = 0;
T = T_guess;

while res>tol
    iter = iter + 1,

    %%%%%%%%%%%%%% Start on mesh 1
    for j=1:2                          % Solve A1*T=b, with only a few GS iterations
        iter_fine = iter_fine + 1;
        T = GS_ONE(A1, b, T);
        res = sum(abs(A*T - b)) / sum(abs( diag(A).*T ));
        res_vec = [res_vec,res];
    end
    plot(x0,T,'o-'), drawnow, pause(0.2)
    r1 = b - A1*T;                     % Residual vector for A1*T=b
    
    
    %%%%%%%%%%%%%% Go to mesh 2
    r2 = RR1*r1;                       % Interpolate residual r1 onto mesh 2
    e2 = zeros(n2,1);                  % Initial guess
    
    for j=1:2                          % Solve A2*e2=r2 with only a few GS iterations
        e2 = GS_ONE(A2, r2, e2);
    end
    s2 = r2 - A2*e2;                   % Residual vector for A2*e2=r2
    
    
    %%%%%%%%%%%%%% Go to mesh 3
    s3 = RR2*s2;                       % Interpolate residual s2 onto mesh 3
    f3 = GS(A3, s3, zeros(n3,1));      % Solve A3*f3=s3 with GS until convergence
    
    
    %%%%%%%%%%%%%% Come back to mesh 2
    f2 = TT2*f3;                       % Interpolate error f3 onto mesh 2
    e2 = e2+f2;                        % Update e2
    
    for j=1:2                          % Solve A2*e2=r2 with only a few GS iterations
        e2 = GS_ONE(A2, r2, e2);
    end
    
    
    %%%%%%%%%%%%%% Come back to mesh 1
    e1 = TT1*e2;                       % Interpolate error e2 onto mesh 1
    T = T+e1;                          % Update T
    plot(x0,T,'o-'), drawnow, pause(0.2)
    
    for j=1:2                          % Solve A1*T=b, with only a few GS iterations
        iter_fine = iter_fine + 1;
        T = GS_ONE(A1, b, T);
        res = sum(abs(A*T - b)) / sum(abs( diag(A).*T ));
        res_vec = [res_vec,res];
    end    
    plot(x0,T,'o-'), drawnow, pause(0.2)
end

xlabel('x [m]'), ylabel('T [K]')    

figure, hold on, 
plot(1:iter_fine, res_vec, 'm'), 
set(gca,'yscale','log')
box on, 
xlabel('Iterations'), ylabel('Normalized residuals'),

end

%--------------------------------------------------------------------
function x = GS(M, rhs, x_guess)
% Solve the linear system M*x=rhs until convergence, 
% using Gauss-Seidel, and with initial guess x_guess.
   n = length(rhs);
   tol = 1e-10;
   res = 1;
   x = x_guess;   % initial guess value
   iter = 0;
   while res>tol
       iter = iter + 1;
       for i=1:n
           x(i)= (rhs(i) - M(i,:)*x + M(i,i)*x(i)) / M(i,i);
       end
       res = sum(abs(M*x - rhs)) / sum(abs( diag(M).*x ));
   end
end   
%--------------------------------------------------------------------
function x = GS_ONE(M, rhs, x_guess)
% Perform ony ONE Gauss-Seidel iteration for the 
% linear system M*x=rhs, with initial guess x_guess.
   n = length(rhs);
   x = x_guess;   % initial guess value
   for i=1:n
       x(i)= (rhs(i) - M(i,:)*x + M(i,i)*x(i)) / M(i,i);
   end
end
%--------------------------------------------------------------------
