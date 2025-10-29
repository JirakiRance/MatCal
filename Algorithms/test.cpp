#include <iostream>
#include"Jiraki_Algorithm.hpp"
#include"..\Utils\QinJiuShao.hpp"

using namespace MatCal::Utils;
using namespace MatCal::Algorithm;

int main() {
    // 示例：f(x) = x^3 - x - 1
    QinJiuShao poly({
        {3,1},
        {1,-1},
        {0,-1}
    });

    auto f = poly.toFunction();
    auto df = poly.derivative().toFunction();
    poly.derivative().show();

    try {
        // 1. Newton 法
        auto newton_result = Newton::solve(f, df, 1.5);
        std::cout << "Newton root: " << newton_result.root 
                  << " (iterations: " << newton_result.iterations << ")\n";

        // 2. 下山 Newton 法
        auto downhill_result = Newton::solve_downhill(f, df, 1.5);
        std::cout << "Downhill Newton root: " << downhill_result.root << "\n";

        // 3. 割线法（双点）
        auto secant_result = Secant::solve_two_point(f, 1.0, 2.0);
        std::cout << "Secant root: " << secant_result.root << "\n";

        // 4. Picard 迭代法（示例 varphi(x)）
        auto varphi = [](double x) { return std::cbrt(2*x*x - 3*x + 5); };
        auto picard_result = Picard::solveDetailed(varphi, 1.0);
        std::cout << "Picard root: " << picard_result.root << "\n";

    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << "\n";
    }

    return 0;
}
