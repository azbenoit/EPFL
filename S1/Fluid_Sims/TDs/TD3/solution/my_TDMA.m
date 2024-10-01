function [x] = my_TDMA(A,b)

n = length(b);
x = zeros(n,1);
P = zeros(n,1);
Q = zeros(n,1);
P(1) = -A(1,2)/A(1,1);
Q(1) = b(1)/A(1,1);
for i=2:n-1,
    tmp = A(i,i) + A(i,i-1)*P(i-1);
    P(i) =  -A(i,i+1) / tmp;
    Q(i) = (-A(i,i-1)*Q(i-1) + b(i)) / tmp;
end
Q(n) = (-A(n,n-1)*Q(n-1) + b(n)) / (A(n,n) + A(n,n-1)*P(n-1));

x(n) = Q(n);

for i=(n-1):-1:1
    x(i) = P(i)*x(i+1) + Q(i);
end
