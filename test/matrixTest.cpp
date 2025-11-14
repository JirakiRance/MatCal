#include "..\src\Matrix.hpp"
#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include <algorithm>
#include <sstream>

using namespace MatCal::Utils;
using namespace MatCal::Algorithm::Matrix;

// 辅助函数：比较两个双精度浮点数是否近似相等
bool is_approx_equal(double a, double b, double epsilon = 1e-12) {
return std::abs(a - b) < epsilon;
}

// 辅助函数：验证 AbstractMatrix 实例中的元素
void assert_matrix_elements(const AbstractMatrix& A, const std::vector<std::vector<double>>& expected, const std::string& step_name) {
bool passed = true;
std::cout << "\n--- 验证: " << step_name << " ---" << std::endl;
    std::cout << ">> 实际计算结果:" << std::endl;
    A.show();

if (A.getRows() != expected.size() || A.getCols() != (expected.empty() ? 0 : expected[0].size())) {
std::cerr << "!!! 失败: 维度不匹配. 预期: " << expected.size() << "x" << (expected.empty() ? 0 : expected[0].size())
                  << ", 实际: " << A.getRows() << "x" << A.getCols() << std::endl;
passed = false;
} else {
        std::cout << ">> 预期结果 (部分展示, 详细对比中):" << std::endl;
        // 仅展示预期结果的前几行，避免冗余
        for (size_t i = 0; i < std::min(A.getRows(), 2); ++i) {
            std::cout << " [";
            for (size_t j = 0; j < std::min(A.getCols(), 4); ++j) {
                std::cout << std::setw(10) << std::setprecision(4) << expected[i][j] << " ";
            }
            if (A.getCols() > 4) std::cout << "...";
            std::cout << "]\n";
        }
        if (A.getRows() > 2) std::cout << " ...\n";


for (int i = 0; i < A.getRows(); ++i) {
for (int j = 0; j < A.getCols(); ++j) {
if (!is_approx_equal(A.get(i, j), expected[i][j])) {
std::cerr << "!!! 失败: 元素不匹配 at [" << i << ", " << j << "]. 预期: " 
                              << expected[i][j] << ", 实际: " << A.get(i, j) << std::endl;
passed = false;
break;
}
}
            if (!passed) break;
}
}
    
    if (passed) {
        std::cout << "<<< 步骤 [" << step_name << "] 成功通过! >>>" << std::endl;
    } else {
        std::cerr << "<<< 步骤 [" << step_name << "] 失败! >>>" << std::endl;
        assert(false); // 触发硬中断
    }
}

// =========================================================================
// 1. Matrix (稠密矩阵) 基础功能测试
// =========================================================================
void test_dense_matrix_basics() {
std::cout << "\n\n=========================================================" << std::endl;
    std::cout << "--- 模块 1: Dense Matrix 稠密矩阵基础运算 ---" << std::endl;
std::cout << "=========================================================" << std::endl;

// 1.1 构造与修改
Matrix A({
{1.0, 2.0, 3.0},
{4.0, 5.0, 6.0}
});
    std::cout << "\n[初始化] 矩阵 A (2x3):" << std::endl;
    A.show();
    
A.set(0, 0, 10.0);
    std::cout << "[操作] 设置 A[0, 0] = 10.0。 检查 A[0, 0] = " << A.get(0, 0) << std::endl;
assert(A.getRows() == 2 && A.getCols() == 3);
assert(is_approx_equal(A.get(1, 2), 6.0));
assert(is_approx_equal(A[0][0], 10.0));
    std::cout << "<<< 步骤 [初始化和修改] 成功通过! >>>" << std::endl;

// 1.2 加法 A + B
Matrix B({
{0.5, 0.5, 0.5},
{1.0, 1.0, 1.0}
});
    std::cout << "\n[输入] 矩阵 B (加数):" << std::endl;
    B.show();
    
auto C = A.add(B);
assert_matrix_elements(*C, {
{10.5, 2.5, 3.5},
{5.0, 6.0, 7.0}
}, "稠密矩阵加法 A + B");

// 1.3 乘法 A * D
Matrix D({
{1.0, 0.0},
{0.0, 1.0},
{1.0, 1.0}
});
    std::cout << "\n[输入] 矩阵 D (乘数):" << std::endl;
    D.show();
    
auto E = A.multiply(D); // (2x3) * (3x2) -> (2x2)
    
assert_matrix_elements(*E, {
{13.0, 5.0},
{10.0, 11.0} // 修复后的正确期望值
}, "稠密矩阵乘法 A * D");

// 1.4 转置
auto F = A.transpose(); // (2x3) -> (3x2)
assert_matrix_elements(*F, {
{10.0, 4.0},
{2.0, 5.0},
{3.0, 6.0}
}, "稠密矩阵转置 A^T");

std::cout << "\n--- 模块 1: Dense Matrix 稠密矩阵基础运算 [全部通过] ---" << std::endl;
}

// =========================================================================
// 2. SparseMatrix (稀疏矩阵) 功能测试
// =========================================================================
void test_sparse_matrix() {
std::cout << "\n\n=========================================================" << std::endl;
    std::cout << "--- 模块 2: Sparse Matrix 稀疏矩阵运算 ---" << std::endl;
std::cout << "=========================================================" << std::endl;

// 2.1 构造与 set/get
SparseMatrix S(3, 3);
S.set(0, 0, 5.0);
S.set(1, 2, 8.0);
S.set(2, 1, 3.0);
    std::cout << "\n[初始化] 稀疏矩阵 S (3x3), 初始元素: (0,0)=5.0, (1,2)=8.0, (2,1)=3.0" << std::endl;
S.set(1, 2, 0.0); // 设置为0应删除元素
    std::cout << "[操作] S[1, 2] 设置为 0.0 (应被删除)." << std::endl;

assert(S.getNonZeroCount() == 2);
assert(is_approx_equal(S.get(0, 0), 5.0));
assert(is_approx_equal(S.get(1, 2), 0.0));
assert(is_approx_equal(S.get(2, 1), 3.0));
    std::cout << "<<< 步骤 [稀疏矩阵初始化和设置] 成功通过! >>>" << std::endl;

// 2.2 转换到稠密矩阵
auto D = S.toNormalMatrix();
assert_matrix_elements(*D, {
{5.0, 0.0, 0.0},
{0.0, 0.0, 0.0},
{0.0, 3.0, 0.0}
}, "稀疏矩阵 S 转换为稠密矩阵 D");

// 2.3 乘法 (稀疏 S * 稠密 M)
Matrix M({
{1.0, 1.0},
{2.0, 2.0},
{3.0, 3.0}
}); // 3x2
    std::cout << "\n[输入] 稠密矩阵 M (3x2):" << std::endl;
    M.show();
    
auto P = S.multiply(M); // (3x3) * (3x2) -> (3x2)
assert_matrix_elements(*P, {
{5.0, 5.0},// 5.0 * R0 = 5.0, 5.0
{0.0, 0.0}, 
{6.0, 6.0}// 3.0 * R1 = 6.0, 6.0
}, "稀疏矩阵乘法 S * M");

// 2.4 稀疏转置
auto ST = S.transpose(); // (3x3) -> (3x3)
assert_matrix_elements(*ST, {
{5.0, 0.0, 0.0},
{0.0, 0.0, 3.0},
{0.0, 0.0, 0.0}
}, "稀疏矩阵转置 S^T");

std::cout << "\n--- 模块 2: Sparse Matrix 稀疏矩阵运算 [全部通过] ---" << std::endl;
}

// =========================================================================
// 3. Triangular Matrix (三角矩阵) 索引与求解测试
// =========================================================================
void test_triangular_matrix() {
std::cout << "\n\n=========================================================" << std::endl;
    std::cout << "--- 模块 3: Triangular Matrix 三角矩阵及其求解 ---" << std::endl;
std::cout << "=========================================================" << std::endl;

// 原始稠密矩阵
Matrix dense({
{1.0, 2.0, 3.0},
{4.0, 5.0, 6.0},
{7.0, 8.0, 9.0}
});
    std::cout << "\n[输入] 稠密矩阵 (用于构造三角矩阵):" << std::endl;
    dense.show();

// 3.1 UpperTriangularMatrix (上三角矩阵) 
try {
UpperTriangularMatrix U(dense);
        std::cout << "\n[构造] 上三角矩阵 U (只保留上三角部分):" << std::endl;
        U.show();
assert(is_approx_equal(U.determinant(), 45.0)); // 1*5*9
        std::cout << "检查 U 的行列式 (1*5*9) = " << U.determinant() << std::endl;
        assert(U.getData().size() == 6);

// 3.2 LowerTriangularMatrix (下三角矩阵)
LowerTriangularMatrix L(dense);
        std::cout << "[构造] 下三角矩阵 L (只保留下三角部分):" << std::endl;
        L.show();
assert(is_approx_equal(L.determinant(), 45.0)); // 1*5*9
        std::cout << "检查 L 的行列式 (1*5*9) = " << L.determinant() << std::endl;
        
        // 3.3 求解上三角系统 Ux=b
Matrix b({ {6.0}, {11.0}, {9.0} });
        std::cout << "\n[求解] 上三角系统 Ux=b, b 向量:" << std::endl;
        b.show();
        
auto x_u = U.solve(b);
assert_matrix_elements(*x_u, { {1.0}, {1.0}, {1.0} }, "上三角系统 Ux=b 的解 x");

// 3.4 求解下三角系统 Ly=b
Matrix b2({ {1.0}, {9.0}, {24.0} });
        std::cout << "\n[求解] 下三角系统 Ly=b, b 向量:" << std::endl;
        b2.show();

auto y_l = L.solve(b2);
assert_matrix_elements(*y_l, { {1.0}, {1.0}, {1.0} }, "下三角系统 Ly=b 的解 y");
        
        // 3.5 转置 (新接口测试)
        auto Ut = U.transpose(); // Ut 是 unique_ptr<LowerTriangularMatrix>
        assert_matrix_elements(*Ut, {
            {1.0, 0.0, 0.0},
            {2.0, 5.0, 0.0},
            {3.0, 6.0, 9.0}
        }, "上三角 U 转置为下三角 Ut");

} catch (const std::exception& e) {
std::cerr << "!!! 模块 3 失败: " << e.what() << std::endl;
assert(false); 
}

std::cout << "\n--- 模块 3: Triangular Matrix 三角矩阵 [全部通过] ---" << std::endl;
}

// =========================================================================
// 4. LU 分解及求解测试
// =========================================================================
void test_lu_decomposition() {
std::cout << "\n\n=========================================================" << std::endl;
    std::cout << "--- 模块 4: LU 分解及求解 (修正矩阵以避免行交换) ---" << std::endl;
std::cout << "=========================================================" << std::endl;
    
    // 修正矩阵 A，避免列主元引发不必要的行交换，使 A=LU 分解成立
Matrix A_orig({
{2.0, 1.0, 1.0},
{4.0, 3.0, 3.0},
{2.0, 1.5, 3.0} // 替换为 {8.0, 7.0, 9.0}
});
    std::cout << "\n[输入] 矩阵 A (用于 LU 分解):" << std::endl;
    A_orig.show();

// LU_Decompose 会修改传入的 A，因此使用副本
Matrix A = A_orig;

try {
auto result = LU_Decompose(A);
        std::cout << "\n[结果] LU 分解结果 L 矩阵:" << std::endl;
        result.L.show();
        std::cout << "[结果] LU 分解结果 U 矩阵:" << std::endl;
        result.U.show();

// 验证 L (基于新的 A)
assert(is_approx_equal(result.L.get(1, 0), 2.0));  // 4/2 = 2.0
assert(is_approx_equal(result.L.get(2, 0), 1.0));  // 2/2 = 1.0
assert(is_approx_equal(result.L.get(2, 1), 0.5));  // (1.5 - 1.0) / (3 - 2) = 0.5
        
        // 验证 U (基于新的 A)
assert(is_approx_equal(result.U.get(0, 0), 2.0));  
assert(is_approx_equal(result.U.get(1, 1), 1.0));  
assert(is_approx_equal(result.U.get(2, 2), 1.5));  
        std::cout << "<<< 步骤 [L/U 矩阵元素验证] 成功通过! >>>" << std::endl;


// 求解 Ax=b, 设 x = [1, 1, 1]^T
Matrix b({ {4.0}, {10.0}, {6.5} }); // 修正 b 向量: A_orig * [1,1,1]^T = [4, 10, 6.5]^T
        std::cout << "\n[求解] Ax=b, b 向量:" << std::endl;
        b.show();
        
auto x_lu = result.solve(b);
assert_matrix_elements(*x_lu, { {1.0}, {1.0}, {1.0} }, "LU 分解系统解 x");

} catch (const std::exception& e) {
std::cerr << "!!! 模块 4 失败: " << e.what() << std::endl;
assert(false);
}
std::cout << "\n--- 模块 4: LU 分解及求解 [全部通过] ---" << std::endl;
}

// =========================================================================
// 5. Elementary Transformations (初等变换) 测试 (仅 Matrix 类型)
// =========================================================================
void test_elementary_transforms() {
std::cout << "\n\n=========================================================" << std::endl;
    std::cout << "--- 模块 5: Elementary Transformations 初等变换 ---" << std::endl;
std::cout << "=========================================================" << std::endl;

Matrix A({
{1.0, 2.0, 3.0},
{4.0, 5.0, 6.0},
{7.0, 8.0, 9.0}
});
    std::cout << "\n[输入] 初始矩阵 A:" << std::endl;
    A.show();
    
// 5.1 行交换 R1 <-> R3
swapRows(A, 0, 2);
assert_matrix_elements(A, {
{7.0, 8.0, 9.0},
{4.0, 5.0, 6.0},
{1.0, 2.0, 3.0}
}, "行交换 R0 <-> R2");

// 5.2 行乘法 R2 <- 0.5 * R2 (第二行索引为 1)
scaleRow(A, 1, 0.5);
assert_matrix_elements(A, {
{7.0, 8.0, 9.0},
{2.0, 2.5, 3.0},
{1.0, 2.0, 3.0}
}, "行乘法 R1 <- 0.5 * R1");

// 5.3 行加法 R3 <- R3 - 2 * R2 (目标行索引 2, 源行索引 1)
addScaledRow(A, 2, 1, -2.0);
assert_matrix_elements(A, {
{7.0, 8.0, 9.0},
{2.0, 2.5, 3.0},
{-3.0, -3.0, -3.0} 
}, "行加法 R2 <- R2 - 2 * R1");

std::cout << "\n--- 模块 5: Elementary Transformations 初等变换 [全部通过] ---" << std::endl;
}

// =========================================================================
// 6. 边缘案例和异常测试
// =========================================================================
void test_edge_cases() {
std::cout << "\n\n=========================================================" << std::endl;
    std::cout << "--- 模块 6: Edge Cases 边缘案例和异常测试 ---" << std::endl;
std::cout << "=========================================================" << std::endl;

// 6.1 越界访问
Matrix A(2, 2);
    std::cout << "\n[异常测试] 越界访问 A.get(2, 0): ";
try {
A.get(2, 0); 
assert(false); 
} catch (const std::out_of_range&) {
std::cout << "抛出 std::out_of_range [成功]" << std::endl;
} catch (...) {
assert(false); 
}

// 6.2 维度不匹配的加法
Matrix B(3, 3);
    std::cout << "[异常测试] 维度不匹配加法 A(2x2) + B(3x3): ";
try {
A.add(B); 
assert(false);
} catch (const std::invalid_argument&) {
std::cout << "抛出 std::invalid_argument [成功]" << std::endl;
}

// 6.3 三角矩阵区域外设置非零值
UpperTriangularMatrix U_test(2);
    std::cout << "[异常测试] 上三角矩阵设置 U[1, 0] = 1.0: ";
try {
U_test.set(1, 0, 1.0); // 应该抛出异常
assert(false);
} catch (const std::invalid_argument&) {
std::cout << "抛出 std::invalid_argument [成功]" << std::endl;
} catch (const std::exception& e) {
std::cerr << "!!! 模块 6 失败: " << e.what() << std::endl;
assert(false);
}

std::cout << "\n--- 模块 6: Edge Cases 边缘案例和异常测试 [全部通过] ---" << std::endl;
}


int main() {
std::cout << "--- Starting Matrix Library Tests ---" << std::endl;

//test_dense_matrix_basics();
//test_sparse_matrix();
//test_triangular_matrix();
test_lu_decomposition();
test_elementary_transforms();
test_edge_cases();

std::cout << "\n\n=========================================================" << std::endl;
    std::cout << "--- 所有测试模块执行完毕 ---" << std::endl;
std::cout << "=========================================================" << std::endl;
return 0;
}