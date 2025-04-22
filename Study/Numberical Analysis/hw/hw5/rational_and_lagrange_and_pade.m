function rational_and_lagrange_and_pade(N,n,m)
    f = @(x) exp(x);
    syms x;

    x_uniform = linspace(-1, 1, N+1);
    y_uniform = f(x_uniform);
    x_cheb = cos((2*(0:N)+1)*pi/(2*(N+1)));
    y_cheb = f(x_cheb);
    % 有理插值 (n,m)
    xdata = linspace(-1, 1, N+1);
    ydata = f(xdata);
    [p, q] = rational_interp(xdata, ydata, n, m);
    fprintf('有理插值 均匀 R(x) = P(x)/Q(x):\n');
    p_poly = poly2sym(p, x);
    q_poly = poly2sym(q, x);
    disp(vpa(p_poly / q_poly,6));
    xdata = x_cheb;
    ydata = y_cheb;
    [p, q] = rational_interp(xdata, ydata, n, m);
    fprintf('有理插值 Chebyshev R(x) = P(x)/Q(x):\n');
    p_poly = poly2sym(p, x);
    q_poly = poly2sym(q, x);
    disp(vpa(p_poly / q_poly,6));

    % Lagarange 插值
    L_uniform = lagrange_poly(x_uniform, y_uniform);
    L_cheb = lagrange_poly(x_cheb,y_cheb);
    fprintf('\n拉格朗日插值 - 均匀节点:\n');
    disp(vpa(L_uniform,6));
    fprintf('\n拉格朗日插值 - Chebyshev 节点:\n');
    disp(vpa(L_cheb,6));

    % Pade 逼近 (n,m)
    [numPade, denPade] = pade_approx_exp(n, m);
    fprintf('\nPade 逼近 (n=%d, m=%d):\n', n, m);
    disp('分子：');
    disp(numPade);
    disp('分母：');
    disp(denPade);
end

%% Rational interpolation: R(x) = P(x)/Q(x)
function [p, q] = rational_interp(x, y, n, m)
    N = length(x) - 1;
    A = zeros(N+1, n + m + 1);
    for i = 1:N+1
        xi = x(i);
        row = [xi.^(0:n), -y(i)*xi.^(1:m)];
        A(i,:) = row;
    end
    coeff = A\y(:);  % 求解线性方程组
    p = coeff(1:n+1);
    q = [1; coeff(n+2:end)];
end

function L = lagrange_poly(X, y)
    syms x;
    n = length(X);
    L_expr = sym(0);
    for i = 1:n
        Li = sym(1);
        for j = 1:n
            if i ~= j
                Li = Li * (x - X(j)) / (X(i) - X(j));
            end
        end
        L_expr = L_expr + y(i) * Li;
    end
    L_expr = expand(L_expr);
    coeffs_sym = sym2poly(L_expr);
    coeffs_double = double(coeffs_sym);
    L = poly2sym(coeffs_double, x);
end


%% Pade approximation at x=0 for e^x
function [P, Q] = pade_approx_exp(n, m)
    N = n + m;
    syms x;
    f_taylor = taylor(exp(x), x, 'Order', N+1);
    c = flip(sym2poly(f_taylor)); 
    A = zeros(m, m);
    b = zeros(m,1);
    for i = 1:m
        for j = 1:m
            if n+i-j+1 > 0 && n+i-j+1 <= length(c)
                A(i,j) = c(n+i-j+1);
            end
        end
        b(i) = -c(n+i+1);
    end

    den_coef = [1; A\b];
    Q = poly2sym(flip(den_coef'), x);

    num_coef = zeros(n+1, 1);
    for k = 0:n
        s = 0;
        for j = 0:min(k,m)
            s = s + den_coef(j+1)*c(k-j+1);
        end
        num_coef(k+1) = s;
    end
    P = poly2sym(flip(num_coef'), x);
end
