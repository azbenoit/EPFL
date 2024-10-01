A = load('thermomech_TC.mat'); A = A.Problem; A = A.A;
n = size(A,1);
n_it = 35;
Bfun = @(x) A*x;
f = @(Y) logm(Y);
Afun = @(Y) matmat(n,Bfun,Y,f,n_it);
tr = -546786.561681857; %true solution







function t = product_trace(A,B)
%Computes trace(A*B) without computing all entries of A*B
    t = sum(sum(A.*B',2));
end


function y = matvec(n,Afun,x,f,n_it)

    [T,U] = lanczos(n,Afun,n_it,x);
    F = f(T);
    y = norm(x)*U*F(:,1);

end


function C = matmat(n,Afun,B,f,n_it)
%Computes C = f(A)*B by n_it iterations of Lanczos

%m = size(A,1);
    m = size(B,2);
    C = zeros(n,m);
    
    for col = 1:m
    
        C(:,col) = matvec(n,Afun,B(:,col),f,n_it);
    
    end

end