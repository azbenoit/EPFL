function [x,iter] = my_GS(A,b, x_guess)

n = length(b);

tol = 1e-10;
res = 1;

x = x_guess; % initial guess value

iter = 0;
while res>tol
    iter = iter + 1;
    
    for i=1:n
        % [b(i) - (sum of A(i,j)*x(i) for j=1...n and j different from i)]/A(i,i).
        % In order to compute the sum efficiently (without a for loop), 
        % rewrite it as sum of A(i,j)*x(i) for j=1...n  (which can be computed
        % as a single matrix-vector product) minus the (j=i) term.
        x(i)= (b(i) - A(i,:)*x + A(i,i)*x(i)) / A(i,i);
    end
    
    res = sum(abs(A*x - b)) / sum(abs( diag(A).*x ));
end
