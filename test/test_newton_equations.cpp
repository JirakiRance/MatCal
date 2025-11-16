#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include"..\src\Iteration.hpp"

#include<windows.h>

using NewtonForEquations=MatCal::Algorithm::Iteration::NewtonForEquations;

// 假设你的 NewtonForEquations 类和其他矩阵库已经包含

void test_newton_solver() {
    std::cout << "=== 非线性方程组求解器测试 ===" << std::endl;
    
    // 测试案例1：简单的二次方程组
    {
        std::cout << "\n--- 测试案例1: 圆和直线交点 ---" << std::endl;
        
        // 方程组：
        // f1(x,y) = x² + y² - 4 = 0  (半径为2的圆)
        // f2(x,y) = x + y - 2 = 0    (直线)
        // 预期解：有两个交点 (0,2) 和 (2,0)
        
        std::vector<NewtonForEquations::Function> funcs = {
            [](const std::vector<double>& x) { return x[0]*x[0] + x[1]*x[1] - 4.0; },
            [](const std::vector<double>& x) { return x[0] + x[1] - 2.0; }
        };
        
        // 测试第一个解附近
        std::vector<double> initial_guess1 = {0.5, 1.5};
        auto result1 = NewtonForEquations::solve(2, funcs, initial_guess1);
        
        std::cout << "初始猜测: (" << initial_guess1[0] << ", " << initial_guess1[1] << ")" << std::endl;
        std::cout << "找到的解: (" << result1.root[0] << ", " << result1.root[1] << ")" << std::endl;
        std::cout << "迭代次数: " << result1.iterations << std::endl;
        std::cout << "最终误差: " << result1.error << std::endl;
        std::cout << "是否收敛: " << (result1.converged ? "是" : "否") << std::endl;
        std::cout << "消息: " << result1.message << std::endl;
        
        // 验证解的正确性
        double f1 = funcs[0](result1.root);
        double f2 = funcs[1](result1.root);
        std::cout << "验证 - f1(x) = " << f1 << ", f2(x) = " << f2 << std::endl;
        
        // 测试第二个解附近
        std::vector<double> initial_guess2 = {1.5, 0.5};
        auto result2 = NewtonForEquations::solve(2, funcs, initial_guess2);
        
        std::cout << "\n初始猜测: (" << initial_guess2[0] << ", " << initial_guess2[1] << ")" << std::endl;
        std::cout << "找到的解: (" << result2.root[0] << ", " << result2.root[1] << ")" << std::endl;
        std::cout << "迭代次数: " << result2.iterations << std::endl;
        std::cout << "最终误差: " << result2.error << std::endl;
        std::cout << "是否收敛: " << (result2.converged ? "是" : "否") << std::endl;
    }
    
    // 测试案例2：包含指数和三角函数的复杂方程组
    {
        std::cout << "\n--- 测试案例2: 复杂非线性方程组 ---" << std::endl;
        
        // 方程组：
        // f1(x,y) = e^x + sin(y) - 2 = 0
        // f2(x,y) = x^2 + cos(y) - 1 = 0
        // 预期解在 (0, π/2) 附近
        
        std::vector<NewtonForEquations::Function> funcs = {
            [](const std::vector<double>& x) { return std::exp(x[0]) + std::sin(x[1]) - 2.0; },
            [](const std::vector<double>& x) { return x[0]*x[0] + std::cos(x[1]) - 1.0; }
        };
        
        std::vector<double> initial_guess = {0.5, 1.0};
        auto result = NewtonForEquations::solve(2, funcs, initial_guess);
        
        std::cout << "初始猜测: (" << initial_guess[0] << ", " << initial_guess[1] << ")" << std::endl;
        std::cout << "找到的解: (" << result.root[0] << ", " << result.root[1] << ")" << std::endl;
        std::cout << "迭代次数: " << result.iterations << std::endl;
        std::cout << "最终误差: " << result.error << std::endl;
        std::cout << "是否收敛: " << (result.converged ? "是" : "否") << std::endl;
        std::cout << "消息: " << result.message << std::endl;
        
        // 验证解的正确性
        double f1 = funcs[0](result.root);
        double f2 = funcs[1](result.root);
        std::cout << "验证 - f1(x) = " << f1 << ", f2(x) = " << f2 << std::endl;
    }
    
    // 测试案例3：三维方程组
    {
        std::cout << "\n--- 测试案例3: 三维方程组 ---" << std::endl;
        
        // 方程组：
        // f1(x,y,z) = x + y + z - 6 = 0
        // f2(x,y,z) = x*y + y*z + z*x - 11 = 0  
        // f3(x,y,z) = x*y*z - 6 = 0
        // 预期解： (1,2,3) 及其排列
        
        std::vector<NewtonForEquations::Function> funcs = {
            [](const std::vector<double>& x) { return x[0] + x[1] + x[2] - 6.0; },
            [](const std::vector<double>& x) { return x[0]*x[1] + x[1]*x[2] + x[2]*x[0] - 11.0; },
            [](const std::vector<double>& x) { return x[0]*x[1]*x[2] - 6.0; }
        };
        
        std::vector<double> initial_guess = {1.5, 2.5, 2.0};
        auto result = NewtonForEquations::solve(3, funcs, initial_guess);
        
        std::cout << "初始猜测: (" << initial_guess[0] << ", " << initial_guess[1] << ", " << initial_guess[2] << ")" << std::endl;
        std::cout << "找到的解: (" << result.root[0] << ", " << result.root[1] << ", " << result.root[2] << ")" << std::endl;
        std::cout << "迭代次数: " << result.iterations << std::endl;
        std::cout << "最终误差: " << result.error << std::endl;
        std::cout << "是否收敛: " << (result.converged ? "是" : "否") << std::endl;
        std::cout << "消息: " << result.message << std::endl;
        
        // 验证解的正确性
        double f1 = funcs[0](result.root);
        double f2 = funcs[1](result.root);
        double f3 = funcs[2](result.root);
        std::cout << "验证 - f1(x) = " << f1 << ", f2(x) = " << f2 << ", f3(x) = " << f3 << std::endl;
    }
    
    // 测试案例4：收敛性测试（困难的初始值）
    {
        std::cout << "\n--- 测试案例4: 收敛性测试 ---" << std::endl;
        
        // 同样的圆和直线方程，但使用较差的初始猜测
        std::vector<NewtonForEquations::Function> funcs = {
            [](const std::vector<double>& x) { return x[0]*x[0] + x[1]*x[1] - 4.0; },
            [](const std::vector<double>& x) { return x[0] + x[1] - 2.0; }
        };
        
        std::vector<double> poor_guess = {10.0, 10.0};  // 很差的初始猜测
        auto result = NewtonForEquations::solve(2, funcs, poor_guess, 50);  // 增加迭代次数
        
        std::cout << "初始猜测: (" << poor_guess[0] << ", " << poor_guess[1] << ")" << std::endl;
        std::cout << "找到的解: (" << result.root[0] << ", " << result.root[1] << ")" << std::endl;
        std::cout << "迭代次数: " << result.iterations << std::endl;
        std::cout << "最终误差: " << result.error << std::endl;
        std::cout << "是否收敛: " << (result.converged ? "是" : "否") << std::endl;
        std::cout << "消息: " << result.message << std::endl;
    }
    
    std::cout << "\n=== 测试完成 ===" << std::endl;
}

int main() {

    SetConsoleOutputCP(CP_UTF8);

    try {
        test_newton_solver();
    } catch (const std::exception& e) {
        std::cerr << "测试过程中发生错误: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}