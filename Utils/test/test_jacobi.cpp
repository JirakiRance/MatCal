#include "..\Matrix.hpp"
#include <iostream>

using namespace MatCal::Utils;
using namespace MatCal::Algorithm::Matrix;


void test_Gauss_Seidel() {
    std::cout << "=== Gauss-Seidel迭代法测试 ===" << std::endl;
    
    // 测试案例1: 3x3矩阵
    Matrix A1 ( {
        {4, -1, 0},
        {-1, 4, -1},
        {0, -1, 4}
    });
    Matrix b1 ( {
        {1},
        {2},
        {3}
    });
    
    std::cout << "系数矩阵 A:" << std::endl;
    A1.show();
    std::cout << "右端向量 b:" << std::endl;
    b1.show();
    
    try {
        auto result1 = Gauss_Seidel(A1, b1, 1e-6, 100);
        std::cout << "Gauss-Seidel迭代结果:" << std::endl;
        std::cout << "是否收敛: " << (result1.converged ? "是" : "否") << std::endl;
        std::cout << "迭代次数: " << result1.iterations << std::endl;
        std::cout << "最终误差: " << result1.error << std::endl;
        std::cout << "解向量 x:" << std::endl;
        result1.root.show();
        
        // 验证解的正确性
        std::cout << "验证 Ax = b:" << std::endl;
        auto Ax = A1.multiply(result1.root);
        Ax->show();
        
    } catch (const std::exception& e) {
        std::cout << "错误: " << e.what() << std::endl;
    }
    
    // 与Jacobi比较收敛速度
    std::cout << "\n=== 与Jacobi方法比较 ===" << std::endl;
    auto jacobi_result = Jacobi(A1, b1, 1e-6, 100);
    auto gauss_result = Gauss_Seidel(A1, b1, 1e-6, 100);
    
    std::cout << "Jacobi迭代次数: " << jacobi_result.iterations << std::endl;
    std::cout << "Gauss-Seidel迭代次数: " << gauss_result.iterations << std::endl;
}

void test_Jacobi() {
    std::cout << "=== Jacobi迭代法测试 ===" << std::endl;
    
    // 测试案例1: 3x3矩阵，有唯一解
    std::cout << "\n测试案例1: 3x3矩阵" << std::endl;
    Matrix A1 ( {
        {4, -1, 0},
        {-1, 4, -1},
        {0, -1, 4}
    });
    Matrix b1  ({
        {1},
        {2},
        {3}
    });
    
    std::cout << "系数矩阵 A:" << std::endl;
    A1.show();
    std::cout << "右端向量 b:" << std::endl;
    b1.show();
    
    try {
        auto result1 = Jacobi(A1, b1, 1e-6, 100);
        std::cout << "Jacobi迭代结果:" << std::endl;
        std::cout << "是否收敛: " << (result1.converged ? "是" : "否") << std::endl;
        std::cout << "迭代次数: " << result1.iterations << std::endl;
        std::cout << "最终误差: " << result1.error << std::endl;
        std::cout << "解向量 x:" << std::endl;
        result1.root.show();
        
        // 验证解的正确性
        std::cout << "验证 Ax = b:" << std::endl;
        auto Ax = A1.multiply(result1.root);
        Ax->show();
        
    } catch (const std::exception& e) {
        std::cout << "错误: " << e.what() << std::endl;
    }
    
    // 测试案例2: 对角占优矩阵
    std::cout << "\n测试案例2: 对角占优矩阵" << std::endl;
    Matrix A2 ({
        {5, 1, 0},
        {1, 5, 1},
        {0, 1, 5}
    });
    Matrix b2 ( {
        {6},
        {7},
        {8}
    });
    
    std::cout << "系数矩阵 A:" << std::endl;
    A2.show();
    std::cout << "右端向量 b:" << std::endl;
    b2.show();
    
    try {
        auto result2 = Jacobi(A2, b2, 1e-6, 100);
        std::cout << "Jacobi迭代结果:" << std::endl;
        std::cout << "是否收敛: " << (result2.converged ? "是" : "否") << std::endl;
        std::cout << "迭代次数: " << result2.iterations << std::endl;
        std::cout << "最终误差: " << result2.error << std::endl;
        std::cout << "解向量 x:" << std::endl;
        result2.root.show();
        
    } catch (const std::exception& e) {
        std::cout << "错误: " << e.what() << std::endl;
    }
    
    // 测试案例3: 多组右端项
    std::cout << "\n测试案例3: 多组右端项" << std::endl;
    Matrix A3 ( {
        {3, 1},
        {1, 3}
    });
    Matrix b3 ( {
        {1, 2},
        {3, 4}
    });
    
    std::cout << "系数矩阵 A:" << std::endl;
    A3.show();
    std::cout << "右端矩阵 b (2组):" << std::endl;
    b3.show();
    
    try {
        auto result3 = Jacobi(A3, b3, 1e-6, 100);
        std::cout << "Jacobi迭代结果:" << std::endl;
        std::cout << "是否收敛: " << (result3.converged ? "是" : "否") << std::endl;
        std::cout << "迭代次数: " << result3.iterations << std::endl;
        std::cout << "最终误差: " << result3.error << std::endl;
        std::cout << "解矩阵 x (2组解):" << std::endl;
        result3.root.show();
        
    } catch (const std::exception& e) {
        std::cout << "错误: " << e.what() << std::endl;
    }
    
    // 测试案例4: 不收敛的情况（非对角占优）
    std::cout << "\n测试案例4: 非对角占优矩阵" << std::endl;
    Matrix A4 ( {
        {1, 2, 3},
        {4, 5, 6},
        {7, 8, 9}
    });
    Matrix b4 ( {
        {1},
        {2},
        {3}
    });
    
    std::cout << "系数矩阵 A:" << std::endl;
    A4.show();
    std::cout << "右端向量 b:" << std::endl;
    b4.show();
    
    try {
        auto result4 = Jacobi(A4, b4, 1e-6, 50);  // 减少迭代次数
        std::cout << "Jacobi迭代结果:" << std::endl;
        std::cout << "是否收敛: " << (result4.converged ? "是" : "否") << std::endl;
        std::cout << "迭代次数: " << result4.iterations << std::endl;
        std::cout << "最终误差: " << result4.error << std::endl;
        if (!result4.converged) {
            std::cout << "注意: 矩阵非对角占优，Jacobi方法可能不收敛" << std::endl;
        }
        std::cout << "解向量 x:" << std::endl;
        result4.root.show();
        
    } catch (const std::exception& e) {
        std::cout << "错误: " << e.what() << std::endl;
    }
    
    // 测试案例5: 使用稀疏矩阵
    std::cout << "\n测试案例5: 稀疏矩阵" << std::endl;
    std::vector<std::tuple<int, int, double>> triplets = {
        {0, 0, 4.0}, {0, 1, -1.0},
        {1, 0, -1.0}, {1, 1, 4.0}, {1, 2, -1.0},
        {2, 1, -1.0}, {2, 2, 4.0}
    };
    SparseMatrix A5(3, 3, triplets);
    Matrix b5 ( {
        {1},
        {2},
        {3}
    });
    
    std::cout << "稀疏系数矩阵 A:" << std::endl;
    A5.show();
    std::cout << "右端向量 b:" << std::endl;
    b5.show();
    
    try {
        auto result5 = Jacobi(A5, b5, 1e-6, 100);
        std::cout << "Jacobi迭代结果:" << std::endl;
        std::cout << "是否收敛: " << (result5.converged ? "是" : "否") << std::endl;
        std::cout << "迭代次数: " << result5.iterations << std::endl;
        std::cout << "最终误差: " << result5.error << std::endl;
        std::cout << "解向量 x:" << std::endl;
        result5.root.show();
        
    } catch (const std::exception& e) {
        std::cout << "错误: " << e.what() << std::endl;
    }
}

void test_MatrixOperations() {
    std::cout << "\n=== 矩阵运算测试 ===" << std::endl;
    
    // 测试矩阵乘法
    Matrix A ( {
        {1, 2},
        {3, 4}
    });
    Matrix B ( {
        {5, 6},
        {7, 8}
    });
    
    std::cout << "矩阵 A:" << std::endl;
    A.show();
    std::cout << "矩阵 B:" << std::endl;
    B.show();
    
    auto C = A.multiply(B);
    std::cout << "A * B:" << std::endl;
    C->show();
    
    // 测试矩阵加法
    auto D = A.add(B);
    std::cout << "A + B:" << std::endl;
    D->show();
    
    // 测试标量乘法
    auto E = A.scalarMultiply(2.5);
    std::cout << "2.5 * A:" << std::endl;
    E->show();
}

void test1104(){
    Matrix A({
        {2,6},
        {2,1e-4}
    });

    Matrix b1={
        {8},
        {8.00001}
    };

    Matrix b2={
        {8},
        {8.00002}
    };

    auto lu=LU_Decompose(A);
    auto ret_ce=solve_columnElimination(A,b1);
    auto ret_gs=Gauss_Seidel(A,b1);
    //cE
    std::cout<<"\nret_ce is:\n"<<ret_ce->toString()<<"\n\n";
    //GS
    std::cout<<"\nret_gs:\nconverged:"<<ret_gs.converged<<"\nerror:"<<ret_gs.error<<"\niterations:"<<ret_gs.iterations<<"\n\nret is:"<<ret_gs.root.toString()<<"\n\n";
    //LU
    std::cout<<"\n\n x1:\n"<<lu.solve(b1)->toString()<<"\n";
    std::cout<<"\n\n x2:\n"<<lu.solve(b2)->toString()<<"\n";
}

void test_isolation() {
    std::cout << "=== 隔离测试 ===" << std::endl;
    
    // 测试1: 只运行CE
    {
        Matrix A({{2,6},{2,1e-4}});
        Matrix b1({{8},{8.00001}});
        auto ret_ce = solve_columnElimination(A, b1);
        std::cout << "单独CE解: " << ret_ce->toString() << std::endl;
    }
    
    // 测试2: 只运行LU  
    {
        Matrix A({{2,6},{2,1e-4}});
        Matrix b1({{8},{8.00001}});
        auto lu = LU_Decompose(A);
        auto x_lu = lu.solve(b1);
        std::cout << "单独LU解: " << x_lu->toString() << std::endl;
    }
    
    // 测试3: 只运行Jacobi
    {
        Matrix A({{2,6},{2,1e-4}});
        Matrix b1({{8},{8.00001}});
        auto ret_jacobi = Jacobi(A, b1);
        std::cout << "单独Jacobi解: " << ret_jacobi.root.toString() << std::endl;
        std::cout << "收敛: " << ret_jacobi.converged << std::endl;
    }
}

void test_with_fresh_objects() {
    std::cout << "=== 使用全新对象测试 ===" << std::endl;
    
    // 场景1: CE开启
    {
        Matrix A_fresh({{2,6},{2,1e-4}});
        Matrix b1_fresh({{8},{8.00001}});

        A_fresh.show();
        b1_fresh.show();
        
        auto ret_ce = solve_columnElimination(A_fresh, b1_fresh);
        std::cout << "CE解: " << ret_ce->toString() << std::endl;

        A_fresh.show();
        b1_fresh.show();
        
        auto ret_jacobi = Jacobi(A_fresh, b1_fresh);
        std::cout << "Jacobi收敛: " << ret_jacobi.converged << std::endl;
        std::cout << "Jacobi解: " << ret_jacobi.root.toString() << std::endl;

        A_fresh.show();
        b1_fresh.show();
        
        auto lu = LU_Decompose(A_fresh);
        auto x_lu = lu.solve(b1_fresh);
        std::cout << "LU解: " << x_lu->toString() << std::endl;

        A_fresh.show();
        b1_fresh.show();
    }
    
    std::cout << "---" << std::endl;
    
    // 场景2: CE关闭（注释掉）
    {
        Matrix A_fresh({{2,6},{2,1e-4}});
        Matrix b1_fresh({{8},{8.00001}});
        
        // auto ret_ce = solve_columnElimination(A_fresh, b1_fresh); // 注释掉
        
        auto ret_jacobi = Jacobi(A_fresh, b1_fresh);
        std::cout << "Jacobi收敛: " << ret_jacobi.converged << std::endl;
        std::cout << "Jacobi解: " << ret_jacobi.root.toString() << std::endl;
        
        auto lu = LU_Decompose(A_fresh);
        auto x_lu = lu.solve(b1_fresh);
        std::cout << "LU解: " << x_lu->toString() << std::endl;
    }
}

int main() {
    std::cout << "矩阵计算库测试程序" << std::endl;
    std::cout << "==================" << std::endl;
    
    try {
        //test_MatrixOperations();
        //test_Jacobi();
        //test_Gauss_Seidel();
        test1104();
        //test_isolation();
        //test_with_fresh_objects();
        
        std::cout << "\n=== 所有测试完成 ===" << std::endl;
        
    } catch (const std::exception& e) {
        std::cout << "测试过程中发生错误: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}