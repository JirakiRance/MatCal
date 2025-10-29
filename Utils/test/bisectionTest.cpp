#include "..\QinJiuShao.hpp"
#include "..\Algorithms\Jiraki_Algorithm.hpp"
#include <iostream>
#include <cmath>
#include <iomanip>

void testSimplePolynomial() {
    std::cout << "=== 多项式求根测试 ===" << std::endl;
    std::cout << std::fixed << std::setprecision(10);
    
    // 测试：求解 2x³ - 4x² - 3x + 6 = 0
    // 这个多项式在区间 [1, 2] 有一个根
    std::cout << "\n求解多项式: 2x³ - 4x² - 3x + 6 = 0" << std::endl;
    
    // 创建多项式
    MatCal::Utils::QinJiuShao poly({{3, 2.0}, {2, -4.0}, {1, -3.0}, {0, 6.0}});
    
    std::cout << "多项式表达式: ";
    poly.show();
    
    // 转换为函数对象
    auto polyFunc = poly.toFunction();
    
    // 使用二分法求根
    std::cout << "\n使用二分法在区间 [1.0, 2.0] 求根..." << std::endl;
    
    try {
        auto result = MatCal::Algorithm::Bisection::solveDetailed(polyFunc, 1.0, 2.0, 1e-12, 100);
        
        std::cout << "求根结果:" << std::endl;
        std::cout << "  根: " << result.root << std::endl;
        std::cout << "  迭代次数: " << result.iterations << std::endl;
        std::cout << "  函数值: f(" << result.root << ") = " << result.error << std::endl;
        std::cout << "  收敛状态: " << (result.converged ? "成功收敛" : "未收敛") << std::endl;
        
        // 验证结果
        std::cout << "\n验证:" << std::endl;
        std::cout << "  f(1.0) = " << polyFunc(1.0) << std::endl;
        std::cout << "  f(2.0) = " << polyFunc(2.0) << std::endl;
        std::cout << "  f(" << result.root << ") = " << polyFunc(result.root) << std::endl;
        
        // 手动验证多项式计算
        std::cout << "\n手动验证多项式计算:" << std::endl;
        double x = result.root;
        double manual = 2.0*x*x*x - 4.0*x*x - 3.0*x + 6.0;
        std::cout << "  手动计算: " << manual << std::endl;
        std::cout << "  秦九韶计算: " << polyFunc(x) << std::endl;
        std::cout << "  差异: " << std::abs(manual - polyFunc(x)) << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "错误: " << e.what() << std::endl;
    }
}

void testComplexPolynomial() {
    std::cout << "=== 高次多项式测试 ===" << std::endl;
    std::cout << std::fixed << std::setprecision(12);
    
    // 12项复杂多项式：x¹¹ - 2x¹⁰ + 3x⁹ - 4x⁸ + 5x⁷ - 6x⁶ + 7x⁵ - 8x⁴ + 9x³ - 10x² + 11x - 12
    std::cout << "\n求解多项式: x¹¹ - 2x¹⁰ + 3x⁹ - 4x⁸ + 5x⁷ - 6x⁶ + 7x⁵ - 8x⁴ + 9x³ - 10x² + 11x - 12 = 0" << std::endl;
    
    MatCal::Utils::QinJiuShao poly({
        {11, 1.0},   // x¹¹
        {10, -2.0},  // -2x¹⁰
        {9, 3.0},    // 3x⁹
        {8, -4.0},   // -4x⁸
        {7, 5.0},    // 5x⁷
        {6, -6.0},   // -6x⁶
        {5, 7.0},    // 7x⁵
        {4, -8.0},   // -8x⁴
        {3, 9.0},    // 9x³
        {2, -10.0},  // -10x²
        {1, 11.0},   // 11x
        {0, -12.0}   // -12
    });
    
    std::cout << "多项式表达式: ";
    poly.show();
    
    std::cout << "\n多项式信息:" << std::endl;
    std::cout << "  最高次数: " << poly.getHighestDegree() << std::endl;
    std::cout << "  项数: " << poly.size() << std::endl;
    
    auto polyFunc = poly.toFunction();
    
    // 在区间 [1.0, 2.0] 寻找根
    std::cout << "\n使用二分法在区间 [1.0, 2.0] 求根..." << std::endl;
    
    try {
        auto result = MatCal::Algorithm::Bisection::solveDetailed(polyFunc, 1.0, 2.0, 1e-10, 200);
        
        std::cout << "求根结果:" << std::endl;
        std::cout << "  根: " << result.root << std::endl;
        std::cout << "  迭代次数: " << result.iterations << std::endl;
        std::cout << "  最终误差: |f(root)| = " << result.error << std::endl;
        std::cout << "  收敛状态: " << (result.converged ? "成功收敛" : "未收敛") << std::endl;
        
        // 验证区间端点值
        std::cout << "\n区间端点验证:" << std::endl;
        std::cout << "  f(1.0) = " << polyFunc(1.0) << std::endl;
        std::cout << "  f(1.5) = " << polyFunc(1.5) << std::endl;
        std::cout << "  f(2.0) = " << polyFunc(2.0) << std::endl;
        std::cout << "  f(" << result.root << ") = " << polyFunc(result.root) << std::endl;
        
        // 验证多项式计算正确性
        std::cout << "\n多项式计算验证:" << std::endl;
        double x = result.root;
        double manual = std::pow(x, 11) - 2*std::pow(x, 10) + 3*std::pow(x, 9) - 4*std::pow(x, 8) + 
                       5*std::pow(x, 7) - 6*std::pow(x, 6) + 7*std::pow(x, 5) - 8*std::pow(x, 4) + 
                       9*std::pow(x, 3) - 10*std::pow(x, 2) + 11*x - 12;
        std::cout << "  手动计算: " << manual << std::endl;
        std::cout << "  秦九韶计算: " << polyFunc(x) << std::endl;
        std::cout << "  计算差异: " << std::abs(manual - polyFunc(x)) << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "错误: " << e.what() << std::endl;
    }
}

int main() {
    testComplexPolynomial();
    return 0;
}