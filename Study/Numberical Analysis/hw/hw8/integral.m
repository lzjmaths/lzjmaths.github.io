%==================================================================
%   修正：四步显式 Adams–Bashforth (AB4) 与 三步隐式 Adams–Moulton (AM3)
%   用于 y' = y - x^2 + 1,  x∈[0,2],  y(0) = 1/2
%==================================================================

clear; clc;
% 清空工作区变量，清除命令窗口

%% 1. 定义方程右端函数 f(x,y) 以及精确解函数 y_exact(x)
f = @(x, y) y - x.^2 + 1;
% f：匿名函数，表示 y' = y - x^2 + 1

y_exact = @(x) -0.5 * exp(x) + x.^2 + 2.*x + 1;
% y_exact：精确解的匿名函数 y(x) = -0.5 e^x + x^2 + 2x + 1

%% 2. 设置步长数组以及预分配误差存储空间
h_values = [1/8, 1/16, 1/32];
% h_values：待测试的三个步长

errors_AB4 = zeros(size(h_values));
errors_AM3 = zeros(size(h_values));
% errors_AB4, errors_AM3：用于存储不同步长下 AB4 和 AM3 的最大误差

%% 3. 对每个步长分别计算
for idx = 1 : length(h_values)
    h = h_values(idx);       % 当前步长 h
    N = 2 / h;               % 区间 [0,2] 上总步数 N
    x = 0 : h : 2;           % 网格节点向量，共 N+1 个点（索引从 1 到 N+1）

    % 初始化数值解向量
    y_AB4 = zeros(1, N+1);
    y_AM3 = zeros(1, N+1);

    % 设置初值 y(0) = 1/2
    y_AB4(1) = 0.5;
    y_AM3(1) = 0.5;

    %% 3.1 用四阶 Runge–Kutta 方法 RK4 生成前 3 个点
    % 为了为 AB4 和 AM3 提供启动所需的 y(1), y(2), y(3), y(4)
    for n = 1 : 3
        % ———— RK4 更新 y_AB4 ————
        k1 = f(x(n),           y_AB4(n));
        k2 = f(x(n) + h/2,     y_AB4(n) + h*k1/2);
        k3 = f(x(n) + h/2,     y_AB4(n) + h*k2/2);
        k4 = f(x(n) + h,       y_AB4(n) + h*k3);
        y_AB4(n+1) = y_AB4(n) + h * (k1 + 2*k2 + 2*k3 + k4) / 6;
        % 四阶 Runge–Kutta 增量：h*(k1 + 2*k2 + 2*k3 + k4)/6

        % ———— RK4 更新 y_AM3 ————
        k1 = f(x(n),           y_AM3(n));
        k2 = f(x(n) + h/2,     y_AM3(n) + h*k1/2);
        k3 = f(x(n) + h/2,     y_AM3(n) + h*k2/2);
        k4 = f(x(n) + h,       y_AM3(n) + h*k3);
        y_AM3(n+1) = y_AM3(n) + h * (k1 + 2*k2 + 2*k3 + k4) / 6;
    end
    % 到这里，y_AB4(1:4) 和 y_AM3(1:4) 都已经用 RK4 计算完成
    %% 3.2 从 n = 3 开始，用 AM3 更新；从 n = 4 开始，用 AB4 更新
    for n = 3 : N
        %% 3.2.1 三步隐式 Adams–Moulton (AM3) 更新
        % AM3 公式：y_{n+1} = y_n + h/24*( 9 f_{n+1} + 19 f_n - 5 f_{n-1} + f_{n-2} )
        % 对于 f(x_{n+1}, y_{n+1})，这里因为 f(x,y) = y - x^2 +1 关于 y 是线性的，
        % 可将 y_{n+1} 显式解出。下面先把含 f_{n+1} 的常数项和已知项合并为 R，再除以系数求 y_{n+1}。
        R = y_AM3(n) ...
          + (h/24) * (  19 * f(x(n),   y_AM3(n)) ...
                      -  5 * f(x(n-1), y_AM3(n-1)) ...
                      +      f(x(n-2), y_AM3(n-2)) ...
                      +  9 * ( - x(n+1)^2 + 1 ) );
        % 其中 9*(-x(n+1)^2+1) 来自于 9 * f(x_{n+1}, 0) 部分（先把 y_{n+1} 分离出去）

        denom_AM3 = 1 - (9 * h / 24);
        y_AM3(n+1) = R / denom_AM3;
        % 直接解析地解出 y_{n+1}

        %% 3.2.2 四步显式 Adams–Bashforth (AB4) 更新（仅当 n ≥ 4）
        if n >= 4
            % AB4 公式：y_{n+1} = y_n + h/24*( 55 f_n - 59 f_{n-1} + 37 f_{n-2} - 9 f_{n-3} )
            y_AB4(n+1) = y_AB4(n) ...
                + (h/24) * (  55 * f(x(n),   y_AB4(n)) ...
                            - 59 * f(x(n-1), y_AB4(n-1)) ...
                            + 37 * f(x(n-2), y_AB4(n-2)) ...
                            -  9 * f(x(n-3), y_AB4(n-3)) );
        end
        % 注意：当 n=3 时，n-3 = 0，会造成索引 x(0)、y(0) 错误，所以必须从 n=4 开始做 AB4。
    end

    %% 3.3 计算数值解与精确解在各节点处的最大范数误差
    err_AB4 = max( abs( y_AB4 - y_exact(x) ) );
    err_AM3 = max( abs( y_AM3 - y_exact(x) ) );
    y_AB4
    y_AM3
    errors_AB4(idx) = err_AB4;
    errors_AM3(idx) = err_AM3;
end

%% 4. 计算相邻步长误差比并估计收敛阶
orders_AB4 = zeros(1, length(h_values)-1);
orders_AM3 = zeros(1, length(h_values)-1);

for idx = 1 : length(h_values)-1
    orders_AB4(idx) = log( errors_AB4(idx) / errors_AB4(idx+1) ) / log(2);
    orders_AM3(idx) = log( errors_AM3(idx) / errors_AM3(idx+1) ) / log(2);
end

%% 5. 在命令窗口输出结果
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
