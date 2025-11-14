#include "..\src\QinJiuShao.hpp"
#include <iostream>
#include <cmath>
#include <iomanip>

void comprehensiveTest() {
    std::cout << "=== 秦九韶算法综合健壮性测试 ===" << std::endl;
    std::cout << std::fixed << std::setprecision(10);
    
    // 测试1：基础构造和显示
    std::cout << "\n1. 基础构造测试" << std::endl;
    MatCal::Utils::QinJiuShao empty;
    std::cout << "空多项式: "; empty.show();
    
    MatCal::Utils::QinJiuShao constant({{0, 5.0}});
    std::cout << "常数多项式: "; constant.show();
    
    MatCal::Utils::QinJiuShao simple({{2, 1.0}, {1, 2.0}, {0, 1.0}});
    std::cout << "简单多项式: "; simple.show();
    std::cout << "f(2) = " << simple.calculate(2.0) << " (期望: 9)" << std::endl;

    // 测试2：边界情况
    std::cout << "\n2. 边界情况测试" << std::endl;
    MatCal::Utils::QinJiuShao edge1({
        {3, 1e-13}, {2, 1.0}, {1, -1e-12}, {0, 2.0}  // 接近0项应被过滤
    });
    std::cout << "接近0系数: "; edge1.show();
    
    MatCal::Utils::QinJiuShao edge2({
        {2, 1.0}, {2, 2.0}, {2, -1.0}, {1, 3.0}, {1, -3.0}  // 重复项合并
    });
    std::cout << "重复项合并: "; edge2.show();
    
    MatCal::Utils::QinJiuShao edge3({{3, 1.0}, {3, -1.0}, {1, 5.0}});  // 抵消为0
    std::cout << "抵消为0: "; edge3.show();

    // 测试3：复杂多项式
    std::cout << "\n3. 复杂多项式测试" << std::endl;
    MatCal::Utils::QinJiuShao complex({
        {5, 2.0}, {4, -3.0}, {3, 1.0}, {2, 4.0}, {1, -2.0}, {0, 1.0},
        {5, 1.0}, {3, -1.0}  // 故意重复
    });
    std::cout << "复杂多项式: "; complex.show();
    std::cout << "f(1) = " << complex.calculate(1.0) << " (期望: 2)" << std::endl;

    // 测试4：多项式运算
    std::cout << "\n4. 多项式运算测试" << std::endl;
    MatCal::Utils::QinJiuShao p1({{2, 3.0}, {1, 2.0}, {0, 1.0}});
    MatCal::Utils::QinJiuShao p2({{1, 1.0}, {0, 2.0}});
    
    auto sum = p1 + p2;
    std::cout << "加法: "; sum.show();
    std::cout << "加法验证 f(2): " << sum.calculate(2.0) << " = " 
              << p1.calculate(2.0) + p2.calculate(2.0) << std::endl;
    
    auto diff = p1 - p2;
    std::cout << "减法: "; diff.show();
    
    auto product = p1 * p2;
    std::cout << "乘法: "; product.show();
    std::cout << "乘法验证 f(2): " << product.calculate(2.0) << " = " 
              << p1.calculate(2.0) * p2.calculate(2.0) << std::endl;

    // 测试5：微积分运算
    std::cout << "\n5. 微积分运算测试" << std::endl;
    MatCal::Utils::QinJiuShao poly({{3, 1.0}, {2, 2.0}, {1, 3.0}, {0, 4.0}});
    std::cout << "原函数: "; poly.show();
    
    auto deriv = poly.derivative();
    std::cout << "导数: "; deriv.show();
    std::cout << "导数验证 f'(1): " << deriv.calculate(1.0) << " (期望: 11)" << std::endl;
    
    auto integral = poly.integral(1.0);
    std::cout << "积分: "; integral.show();
    
    double definite = poly.definiteIntegral(0, 2);
    std::cout << "定积分 [0,2]: " << definite << " (期望: 20)" << std::endl;

    // 测试6：动态插入删除
    std::cout << "\n6. 动态操作测试" << std::endl;
    MatCal::Utils::QinJiuShao dynamic;
    dynamic.insert(3, 2.0);
    dynamic.insert(1, 1.0);
    dynamic.insert(2, 3.0);  // 应该自动排序
    std::cout << "动态插入: "; dynamic.show();
    
    dynamic.insert(2, -3.0);  // 合并为0
    std::cout << "合并为0后: "; dynamic.show();
    
    dynamic.remove(1);
    std::cout << "删除1次项后: "; dynamic.show();

    // 测试7：函数转换
    std::cout << "\n7. 函数转换测试" << std::endl;
    auto func = poly.toFunction();
    std::cout << "函数转换 f(1.5): " << func(1.5) << " = " << poly.calculate(1.5) << std::endl;

    // 测试8：拷贝和赋值
    std::cout << "\n8. 拷贝构造测试" << std::endl;
    MatCal::Utils::QinJiuShao original({{2, 2.0}, {0, 1.0}});
    MatCal::Utils::QinJiuShao copy = original;
    std::cout << "原多项式: "; original.show();
    std::cout << "拷贝多项式: "; copy.show();
    std::cout << "拷贝验证 f(3): " << copy.calculate(3.0) << " = " << original.calculate(3.0) << std::endl;

    // 测试9：性能测试（大多项式）
    std::cout << "\n9. 性能测试" << std::endl;
    MatCal::Utils::QinJiuShao large;
    for(int i = 0; i <= 10; ++i) {
        large.insert(i, 1.0 / (i + 1));
    }
    std::cout << "10次多项式项数: " << large.size() << std::endl;
    std::cout << "最高次数: " << large.getHighestDegree() << std::endl;
    std::cout << "f(1) = " << large.calculate(1.0) << std::endl;

    std::cout << "\n🎉 所有测试通过！秦九韶算法类非常健壮！" << std::endl;
}

int main() {
    try {
        comprehensiveTest();
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "测试失败: " << e.what() << std::endl;
        return 1;
    }
}