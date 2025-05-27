clear; clc;

f = @(x, y) y - x.^2 + 1;
y_exact = @(x) -0.5 * exp(x) + x.^2 + 2.*x + 1;

h_values = [1/8, 1/16, 1/32];
errors_AB4 = zeros(size(h_values));
errors_AM3 = zeros(size(h_values));

for idx = 1 : length(h_values)
    h = h_values(idx);
    N = 2 / h;
    x = 0 : h : 2;

    y_AB4 = zeros(1, N+1);
    y_AM3 = zeros(1, N+1);

    y_AB4(1) = 0.5;
    y_AM3(1) = 0.5;

    for n = 1 : 3
        k1 = f(x(n),           y_AB4(n));
        k2 = f(x(n) + h/2,     y_AB4(n) + h*k1/2);
        k3 = f(x(n) + h/2,     y_AB4(n) + h*k2/2);
        k4 = f(x(n) + h,       y_AB4(n) + h*k3);
        y_AB4(n+1) = y_AB4(n) + h * (k1 + 2*k2 + 2*k3 + k4) / 6;

        k1 = f(x(n),           y_AM3(n));
        k2 = f(x(n) + h/2,     y_AM3(n) + h*k1/2);
        k3 = f(x(n) + h/2,     y_AM3(n) + h*k2/2);
        k4 = f(x(n) + h,       y_AM3(n) + h*k3);
        y_AM3(n+1) = y_AM3(n) + h * (k1 + 2*k2 + 2*k3 + k4) / 6;
    end

    for n = 3 : N
        R = y_AM3(n) ...
          + (h/24) * (  19 * f(x(n),   y_AM3(n)) ...
                      -  5 * f(x(n-1), y_AM3(n-1)) ...
                      +      f(x(n-2), y_AM3(n-2)) ...
                      +  9 * ( - x(n+1)^2 + 1 ) );

        denom_AM3 = 1 - (9 * h / 24);
        y_AM3(n+1) = R / denom_AM3;

        if n >= 4
            y_AB4(n+1) = y_AB4(n) ...
                + (h/24) * (  55 * f(x(n),   y_AB4(n)) ...
                            - 59 * f(x(n-1), y_AB4(n-1)) ...
                            + 37 * f(x(n-2), y_AB4(n-2)) ...
                            -  9 * f(x(n-3), y_AB4(n-3)) );
        end
    end

    err_AB4 = max( abs( y_AB4 - y_exact(x) ) );
    err_AM3 = max( abs( y_AM3 - y_exact(x) ) );
    y_AB4
    y_AM3
    errors_AB4(idx) = err_AB4;
    errors_AM3(idx) = err_AM3;
end

orders_AB4 = zeros(1, length(h_values)-1);
orders_AM3 = zeros(1, length(h_values)-1);

for idx = 1 : length(h_values)-1
    orders_AB4(idx) = log( errors_AB4(idx) / errors_AB4(idx+1) ) / log(2);
    orders_AM3(idx) = log( errors_AM3(idx) / errors_AM3(idx+1) ) / log(2);
end

fprintf('=== 步长 h       AB4 最大误差      AM3 最大误差 ===\n');
for idx = 1:length(h_values)
    fprintf('h = %-6.4f    %.4e    %.4e\n', ...
        h_values(idx), errors_AB4(idx), errors_AM3(idx));
end
fprintf('\n');
fprintf('=== 估计收敛阶 (相邻 h 之间) ===\n');
for idx = 1:length(orders_AB4)
    fprintf('h=%.4f->%.4f, AB4 阶 ≈ %.4f,   AM3 阶 ≈ %.4f\n', ...
        h_values(idx), h_values(idx+1), orders_AB4(idx), orders_AM3(idx));
end
