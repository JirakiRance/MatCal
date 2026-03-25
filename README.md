# MatCal C++ 库完整详细 API 文档

本手册提供了 `MatCal` 库所有命名空间、类和函数的极度详细说明，适用于二次开发与深入查阅。

---

## 目录
[toc]

---

## 1. 命名空间：MatCal::Utils (核心数据结构)

### 1.1 `QinJiuShaoNode` (秦九韶结点)
用于构建秦九韶多项式的基本单元。
* **成员变量**:
  * `double a`: 系数。
  * `int n`: 次数 (要求 $\ge 0$)。
  * `static constexpr double ZERO_THRESHOLD = 1e-12`: 零阈值，绝对值小于此值的系数被视为 0。
* **方法**:
  * `bool isFake() const`: 判断该项系数是否小于 `ZERO_THRESHOLD`（是否为无效项）。
  * `std::string toString() const`: 格式化输出为字符串（如 `3.5X^2`）。
  * 包含完整的比较运算符重载 (`>`, `<`, `==`, `!=`, `>=`, `<=`)，主要按**次数**排序。

### 1.2 `QinJiuShao` (秦九韶多项式)
利用秦九韶算法（Horner's Method）实现的高效多项式类。自动维护内部按次数降序排列，且自动合并同类项。
* **构造函数**:
  * `QinJiuShao(std::vector<QinJiuShaoNode> par = {})`: 基础构造，会自动调用内部 `cleanup()` 整理。
  * `QinJiuShao(std::initializer_list<std::pair<int, double>> pairs)`: 方便的代码内联构造，如 `{{2, 3.0}, {0, 1.0}}` ($3x^2+1$)。
  * `QinJiuShao(int A, std::vector<double>& zeros)`: 从最高次项系数 `A` 和一系列零点 (根) 构造展开多项式。
* **增删改查**:
  * `void insert(int n, double a)` / `insert(QinJiuShaoNode node)`: 插入新项，若存在同次项则累加系数。使用二分查找，复杂度 $O(\log N)$。
  * `bool remove(int n)`: 删除指定次数 `n` 的项。
  * `double calculate(double x) const`: 极速求值，复杂度 $O(N)$。
  * `std::function<double(double)> toFunction() const`: 将该多项式打包为标准的 `std::function`，方便向下传递。
* **多项式运算**:
  * `operator+`, `operator-`: 多项式加减法。
  * `operator*(const QinJiuShao& other)`: 多项式乘法，复杂度 $O(N^2)$。
  * `operator*(double scalar)`, `operator/(double scalar)`: 标量乘除。
* **微积分**:
  * `QinJiuShao derivative() const`: 对原多项式求导，返回新多项式。
  * `QinJiuShao integral(double constant = 0) const`: 求不定积分，`constant` 为常数项 $C$。
  * `double definiteIntegral(double a, double b) const`: 计算区间 $[a, b]$ 的定积分。
  * `double product(QinJiuShao& other, double a, double b)`: 计算两多项式在 $[a, b]$ 上的内积。

### 1.3 矩阵抽象基类 `AbstractMatrix`
所有矩阵类型的父类，提供多态接口。
* **核心纯虚函数** (子类必须实现):
  * `void resize(int newRows, int newCols)`
  * `double get(int row, int col) const`
  * `void set(int row, int col, double value)`
  * `std::unique_ptr<AbstractMatrix> toNormalMatrix() const`: 统一转换为稠密矩阵指针。
  * `add`, `multiply`, `scalarMultiply`: 基础四则运算。

### 1.4 `Matrix` (稠密普通矩阵)
使用 `std::vector<std::vector<double>>` 存储数据的标准稠密矩阵。
* **特有方法**:
  * `static Matrix identity(int size)`: 生成单位矩阵。
  * `static Matrix zeros(int rows, int cols)` / `ones(...)`: 生成全0或全1矩阵。
  * `std::unique_ptr<AbstractMatrix> transpose() const`: 转置。
  * `std::vector<double>& operator[](int row)`: 支持 `A[i][j]` 形式的快速引用访问。
  * `std::string toString() const`: 返回标准的 JSON 格式字符串。

### 1.5 `SparseMatrix` (稀疏矩阵)
使用 `(row, col, value)` 三元组列表存储，节省大量内存空间。
* **特有方法**:
  * `SparseMatrix(int rows, int cols, const std::vector<std::tuple<int,int,double>>& triplets)`: 从三元组构造。
  * `double getSparsity() const`: 计算稀疏度（$1 - 非零元素数 / 总容量$）。
  * `static SparseMatrix random(int rows, int cols, double density)`: 生成指定非零密度的随机稀疏矩阵。
* **注意**: 其 `multiply` 方法已针对稀疏结构进行了内部双重循环优化。

### 1.6 `UpperTriangularMatrix` & `LowerTriangularMatrix`
上下三角矩阵，通过一维 `std::vector<double>` 压缩存储。尝试在非存储区域写入非零值会抛出异常。
* **特有求解方法**:
  * `std::unique_ptr<AbstractMatrix> solve(const AbstractMatrix& b) const`: 
    * 对上三角：利用**回代法 (Back Substitution)** 求解 $Ux=b$。
    * 对下三角：利用**前代法 (Forward Substitution)** 求解 $Lx=b$。
    * 支持 $b$ 为多列（即一次求解多组右端向量）。
* **特有计算**:
  * `double determinant() const`: 行列式，即主对角线元素之积。
  * `bool isSingular() const`: 检查主对角线是否有 0 元素。

### 1.7 `TridiagonalMatrix` (三对角矩阵)
仅存储主对角线 `d`、下对角线 `l` 和上对角线 `u`。
* **特有方法**:
  * `std::unique_ptr<AbstractMatrix> solve(const AbstractMatrix& B) const`: 
    * 使用**追赶法 (Thomas Algorithm)** 求解 $Tx = B$。复杂度仅为 $O(N)$，彻底避免高斯消元 $O(N^3)$ 的耗时。

---

## 2. 命名空间：MatCal::Algorithm::Matrix (矩阵高级算法)

### 2.1 矩阵初等变换与范数
* `void swapRows(AbstractMatrix& A, int r1, int r2)`: 行交换。若非普通 `Matrix` 会触发自动转换。
* `void scaleRow(...)` / `addScaledRow(...)`: 行缩放与行叠加变换。
* `void applyRowSwaps(AbstractMatrix& M, const std::vector<std::pair<int,int>>& swaps, bool reverse)`: 按顺序应用或逆向还原行交换记录。
* `double norm_one(...)`, `norm_infinite(...)`, `norm_Frobenius(...)`: 分别计算列和范数 (1-范数)、行和范数 ($\infty$-范数) 和 Frobenius 范数。

### 2.2 线性方程组求解
* **`solve_columnElimination(AbstractMatrix& A, AbstractMatrix& b)`**
  * **机制**: 列主元高斯消去法。通过寻找当前列绝对值最大的元素作为主元进行行交换，有效减少计算机浮点截断误差。
  * **返回**: 包含解矩阵 $x$ 的智能指针。
* **`LU_Decompose(AbstractMatrix& A, double eps)`**
  * **机制**: Doolittle LU 分解 ($A = LU$，其中 $L$ 的对角线为 1)。
  * **返回**: `LUresult` 结构体，内含 `L` 和 `U`。可调用 `result.solve(x)` 极速求解。

### 2.3 迭代求解法
用于求解 $Ax=b$。共同返回值：`Iteration_Result` 结构体（含 `root` 解、`iterations` 次数、`error` 最终误差、`converged` 是否收敛布尔值）。
* **`Jacobi(A, b, epsilon, max_iterations)`**: 雅可比迭代法。
* **`Gauss_Seidel(A, b, epsilon, max_iterations)`**: 高斯-赛德尔迭代法（利用已更新的值计算后续值）。
* **`SOR(A, b, omega, epsilon, max_iterations)`**: 逐次超松弛迭代法。`omega` 需在 $(0, 2)$ 范围内。

### 2.4 特征值求解
* **`PowerMethod(A, eps, max_iter)`**: 规范化幂法。求矩阵按模最大的特征值 $\lambda_1$ 及对应的主特征向量。返回 `PowerMethod_Result`。
* **`PowerMethod_reverse(A, near_num, eps, max_iter)`**: 反幂法。求矩阵距离 `near_num` 最近的特征值。内部结合 LU 分解求解线性系统。

---

## 3. 命名空间：MatCal::Algorithm::Basics (基础数值算法)

### 3.1 `Derivative` (数值求导)
* `double dy_dx(Func_y, x, eps)`: 一元函数前向差分：$(f(x+h)-f(x))/h$。
* `double dy_dx_center(Func_y, x, eps)`: 一元函数中心差分：$(f(x+h/2)-f(x-h/2))/h$，精度更高。
* `double pF_px(Func_F, xs, i, eps)`: 多元函数偏导数（针对输入向量的第 `i` 个分量）。
* `double dF_dx(...)`: 多元函数全导数（各偏导数之和）。

### 3.2 `Least_Square` (最小二乘多项式拟合)
* **参数说明**:
  * `degree`: 期望拟合的多项式最高次数。
  * `x`, `y`: 数据点集。
  * `weights`: 可选。对应数据点的权重数组。
  * `selects`: 可选。布尔数组（长度为 `degree+1`），为 `true` 表示强制在模型中包含该次项，`false` 表示剔除（如只拟合偶数次项）。
* **返回**: `Result_least_square`。包含：系数 `coee` 数组、法方程矩阵 `A` 和 `b`、秦九韶多项式对象 `poly`。内部尝试优先用 Gauss-Seidel 解法方程，失败则回退至列主元消去。

### 3.3 `OrthogonalPolynomials` (正交多项式)
* `Chebyshev(int n, bool second)`: 利用三项递推公式极速生成第 `n` 阶第一类/第二类切比雪夫多项式。
* `ChebyshevZeros(int n)`: 返回第 `n` 阶第一类切比雪夫多项式在 $[-1, 1]$ 内的所有解析根。
* `Legendre(int n)`: 利用递推关系生成勒让德多项式。

### 3.4 `NumericalIntegration` (数值积分)
* `Instant`: 最基础的黎曼和矩形积分。
* `NewtonCotes`: 牛顿-柯特斯公式积分，参数 `n` (1~7) 决定插值多项式的阶数（如 n=1 为梯形，n=2 为辛普森）。
* `CompositeNewtonCotes`: 复合牛顿-柯特斯公式。通过 `segments` 将大区间均分，在子区间应用 NewtonCotes，精度极高。
* `Romberg`: **龙贝格积分**。通过理查森外推法 (Richardson Extrapolation) 将梯形公式 (T) 逐步加速为辛普森 (S)、柯特斯 (C) 和龙贝格 (R) 序列。返回最终值及推演矩阵 (T-表)。
*(注：上述三类积分方法均提供了针对 `std::function` 和针对离散数据表 `vector<pair>` 的两套重载。)*

### 3.5 `ODE` (常微分方程)
求解形如 $Y' = F(x, Y)$ 的初值问题。
* `SimpleEuler` / `Euler`: 显式欧拉法与改进欧拉法（预测-校正）。
* `RungeKutta_44`: 经典的 4 阶 4 段龙格-库塔法 (RK4)。
  * `funcs`: 右端函数集合。
  * `inits`: 初始值集合 (包含自变量 $x_0$ 和 各状态初值)。
  * **返回**: 记录求解步长的完整状态矩阵。

---

## 4. 命名空间：MatCal::Algorithm::Insert (插值算法)

通用设计：所有类均支持使用 `calculate(double X)` 取得插值，并提供 `reconstruct()` 用于重置数据。

### 4.1 `LagrangeInsert` (拉格朗日插值)
* **机制**: 直接根据公式构造基函数 $l_i(x)$，将所有基函数利用秦九韶乘法累加得到最终的 $L_n(x)$ 多项式对象。

### 4.2 `NewtonInsert_Quotient` (牛顿差商插值 - 不等距)
* **机制**: 构建下三角差商表 `lower`。
* **特有功能**: `insertNewTerm(double X, double Y)` 动态增量更新。由于牛顿插值的承继性，新增节点无需重算整个多项式，直接在表末端拓展即可。

### 4.3 `NewtonInsert_Finite` (牛顿差分插值 - 等距)
* **机制**: 专为等距节点优化，构造时仅需提供步长 `h`、起点 `x0` 以及各点的纵坐标 `yData`。内部构建前向差分表。

### 4.4 `Hermite` (埃尔米特插值)
* **机制**: 满足函数值与导数值双重拟合。构造函数接受 `xs`, `ys` 以及 **导数数组 `dy_dxs`**。通过构造带平方因子的特殊基函数，确保平滑过渡，消除拉格朗日插值常见的龙格现象振荡。

### 4.5 `CubicSpline` (三次样条插值)
* **机制**: 默认采用**自然边界条件**（两端二阶导数为0）。构建三对角矩阵方程 $A \cdot M = d$ 求解弯矩，随后代入样条方程。
* **注意**: 不生成统一的多项式对象，而是使用 `findInterval` 二分查找所在区间后，在对应的三次小区间上求值。

### 4.6 `LinearInsert` (分段线性插值)
* **机制**: 相邻点直接连线。如果给定的 `X` 超出插值范围，会自动利用首末段的斜率进行**线性外推**。

---

## 5. 命名空间：MatCal::Algorithm::Iteration (非线性方程求根)

### 5.1 单变量方程根求解
所有类均提供 `solveDetailed`，返回 `Result` 结构体：含 `root` (根)、`iterations` (迭代步数)、`error` (精度)、`converged` (收敛标识)、`series` (每一步的迭代轨迹数组)。

* **`Bisection`** (二分法): 给定初始区间 $[a, b]$，要求端点异号，保证 100% 收敛。
* **`Picard`** (简单迭代法):
  * 要求将方程改写为 $x = \varphi(x)$ 形式。
  * `solve_Aitken`: 使用**埃特金 (Aitken) 加速**技术，对原迭代序列进行外推，大幅加快缓慢收敛的序列。
* **`Newton`** (牛顿-拉夫逊法):
  * `solve`: $x_{k+1} = x_k - f(x_k)/f'(x_k)$。
  * `solve_downhill`: **牛顿下山法**。若新探测点的函数绝对值大于当前点（即偏离根区），内部触发阻尼机制，将搜索步长 $\lambda$ 不断折半，强制保证下降特性。极大增强了对糟糕初值的鲁棒性。
* **`Secant`** (割线法):
  * `solve_two_point`: 给定两个初值，以割线斜率替代切线斜率。
  * `solve_one_point`: 单点割线（基于小步长 $h$ 差分的牛顿法变体）。

### 5.2 多变量方程组求解
* **`NewtonForEquations::solve`**
  * **输入**: 
    * `int n`: 变量个数。
    * `std::vector<Function>& Funcs`: 包含 $n$ 个多元函数的集合（每个接收 `vector<double>` 并返回 `double`）。
    * `std::vector<double> xs`: 初始猜测向量。
  * **机制**:
    1. 在当前点利用 `Derivative::pF_px` 计算出 $n \times n$ 的雅可比矩阵 (Jacobian) $J$。
    2. 计算右端向量 $-F(x)$。
    3. 构建线性系统 $J \cdot \Delta x = -F(x)$。
    4. **高容错回退机制**: 优先尝试利用列主元高斯消去法解此系统；若遇到矩阵奇异或条件数极差导致抛出异常，内部通过 `try-catch` 捕获，并**自动切换**为 Gauss-Seidel 迭代法强行逼近解。
  * **返回**: 包含解向量、迭代次数与详细错误/成功日志 (`message`) 的结构体。