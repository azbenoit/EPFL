
T = 10000;
computed = [];
expected = [];
ns = 5:50;
for n = 5:50
    sigma = 0.05;
    expected = [expected, sqrt(2*n./pi) * sigma];
    cont = 0;
    for i = 1:T
        x = randn(n, 1);
        x = x / norm(x);
        if (abs(x(1)) < sigma)
            cont = cont + 1/T;
        end
    end
    computed = [computed, cont];
end
subplot(2, 1, 1)
plot(ns, computed, 'Linewidth', 2);
hold on
plot(ns, expected, '--', 'Linewidth', 2);
legend('computed', 'expected');
xlabel('n')
ylabel('failure probability')


clear;
T = 10000;
computed = [];
expected = [];
expected = [];
n = 20;
sigmas = 0.01:0.005:0.1;
for sigma = sigmas
    expected = [expected, sqrt(2*n./pi) * sigma];
    cont = 0;
    for i = 1:T
        x = randn(n, 1);
        x = x / norm(x);
        if (abs(x(1)) < sigma)
            cont = cont + 1/T;
        end
    end
    computed = [computed, cont];
end
subplot(2, 1, 2)
plot(sigmas, computed, 'Linewidth', 2);
hold on
plot(sigmas, expected, '--', 'Linewidth', 2);
legend('computed', 'expected')
ylabel('failure probability')