#include <iostream>
#include "..\src\Matrix.hpp"   // 你的头文件
#include<windows.h>
using namespace MatCal::Utils;

int main(){

    SetConsoleOutputCP(CP_UTF8);
    /*---------- 1. 单右端 ----------*/
    // 三对角线（n=4）
    std::vector<double> low  = { 1.0, 1.0, 1.0 };
    std::vector<double> diag = { 4.0, 4.0, 4.0, 4.0 };
    std::vector<double> up   = { 1.0, 1.0, 1.0 };

    TridiagonalMatrix T(low, diag, up);
    std::cout << "T ="; T.show();

    // 右端向量 b = [5 6 7 8]^T
    Matrix b({{5},{6},{7},{8}});
    auto x = T.solve(b);
    std::cout << "x ="; x->show();

    // 残差 ||T*x - b||_inf
    auto r = T.multiply(*x)->add(*b.scalarMultiply(-1.0));
    double res = 0.0;
    for(int i = 0; i < 4; ++i) res = std::max(res, std::abs(r->get(i,0)));
    std::cout << "单右端残差 ||Tx-b||_inf = " << res << "\n\n";

    /*---------- 2. 多右端 ----------*/
    Matrix B({
        { 5,  10},
        { 6,  12},
        { 7,  14},
        { 8,  16}
    });
    auto X = T.solve(B);
    std::cout << "X (2 右端) ="; X->show();

    // 残差矩阵
    auto R = T.multiply(*X)->add(*B.scalarMultiply(-1.0));
    double res2 = 0.0;
    for(int i = 0; i < 4; ++i)
        for(int j = 0; j < 2; ++j)
            res2 = std::max(res2, std::abs(R->get(i,j)));
    std::cout << "多右端残差 ||TX-B||_inf = " << res2 << '\n';
    return 0;
}