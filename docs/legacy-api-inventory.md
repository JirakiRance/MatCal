# Legacy API Inventory

This inventory describes the current public legacy surface. It is not a promise that every design is good; it is the M0 compatibility baseline.

## Public Headers

- `src/Matrix.hpp`
- `src/QinJiuShao.hpp`
- `src/Insert.hpp`
- `src/Iteration.hpp`
- `src/Basics.hpp`

## Namespaces

- `MatCal::Utils`
- `MatCal::Algorithm::Matrix`
- `MatCal::Algorithm::Insert`
- `MatCal::Algorithm::Iteration`
- `MatCal::Algorithm::Basics`
- `MatCal::Algorithm::Basics::Integrate`

## `MatCal::Utils`

### `QinJiuShaoNode`

- Public data: `double a`, `int n`, `static constexpr double ZERO_THRESHOLD = 1e-12`.
- Constructors: `QinJiuShaoNode(int nn = 0, double aa = 0)`.
- Methods: `isFake() const`, `toString() const`.
- Operators: `>`, `<`, `==`, `!=`, `>=`, `<=`.
- Errors: negative degree prints to stdout and throws `std::invalid_argument`.
- Ownership/lifetime: value type.

### `QinJiuShao`

- Constructors:
  - `QinJiuShao(std::vector<QinJiuShaoNode> par = {})`
  - `QinJiuShao(std::vector<std::pair<int,double>>& vec)`
  - `QinJiuShao(std::initializer_list<std::pair<int,double>> pairs)`
  - `QinJiuShao(int A, std::vector<double>& zeros)`
  - `QinJiuShao(const QinJiuShao& qin)`
- Methods: `insert(int,double)`, `insert(QinJiuShaoNode)`, `remove(int)`, `remove_if(...)`, `calculate(double) const`, `show() const`, `toString() const`, `getParameters() const`, `toFunction() const`, `getHighestDegree() const`, `getParameter(int) const`, `size() const`, `empty() const`, `clear()`.
- Operators: `+`, `-`, polynomial `*`, scalar `*`, scalar `/`, friend scalar-left `*`.
- Calculus: `derivative() const`, `integral(double constant = 0) const`, `definiteIntegral(double,double) const`, `product(QinJiuShao& other,double,double)`.
- Errors: empty highest degree throws; scalar division by near-zero throws.
- Ownership/lifetime: value type. Since M0.1, `toFunction()` captures polynomial coefficients by value; the returned function owns the state it needs and may outlive the polynomial object. Before M0.1 this callable borrowed `this`.

### Matrix Base and Storage Types

`AbstractMatrix`

- Public polymorphic API: `getRows() const`, `getCols() const`, `resize(int,int)`, `get(int,int) const`, `set(int,int,double)`, `show() const`, `toString() const`, `toNormalMatrix() const`, `add(...) const`, `multiply(...) const`, `operator*(const AbstractMatrix&) const`, `scalarMultiply(double) const`, `operator*(double) const`, `isSquare() const`.
- Destructor: yes, virtual pure destructor with definition.
- Ownership/lifetime: returned `std::unique_ptr<AbstractMatrix>` owns the result. Callers usually dynamic-cast to `Matrix`.

`Matrix`

- Constructors: `Matrix(int row = 0,int col = 0)`, `Matrix(const std::vector<std::vector<double>>& input)`, `Matrix(const double** input,int row,int col)`, `Matrix(const Matrix&)`, `Matrix(std::initializer_list<std::initializer_list<double>>)`, `Matrix(AbstractMatrix& A)`.
- Assignment: `operator=(AbstractMatrix&)`, `operator=(const Matrix&)`.
- Methods: overrides for base operations; `transpose() const`, `operator[](int)`, `operator[](int) const`, `getData() const`.
- Static factories: `identity(int)`, `zeros(int,int)`, `ones(int,int)`.
- Errors: negative sizes throw; index errors throw; invalid dimensions throw.
- Ownership: raw `double**` input is copied; MatCal does not take ownership.
- Input modification: arithmetic methods return new matrices. `resize` mutates.

`SparseMatrix`

- Constructors: `SparseMatrix(int rows = 0,int cols = 0)`, `SparseMatrix(const Matrix&)`, `SparseMatrix(int,int,const std::vector<std::tuple<int,int,double>>& triplets)`.
- Public `Element`: `int row`, `int col`, `double value`.
- Methods: base overrides; `transpose() const`, `inner_sort()`, `find(int,int)`, `getNonZeroCount() const`, `getSparsity() const`, `getTriplets() const`, `getElements() const`.
- Static factories: `identity(int)`, `diagonal(const std::vector<double>&)`, `random(int,int,double)`.
- Complexity: `get` and `set` are linear search over triplets; some operations degrade toward dense loops.

`AbstractTriangularMatrix`

- Protected virtual hooks: `calculateIndex`, `isInStoredRegion`, `createEmptyClone`.
- Methods: base overrides; `getData() const`, `determinant() const`, `isSingular() const`.

`UpperTriangularMatrix`

- Constructors: `UpperTriangularMatrix(int size = 0)`, copy/move constructors, `UpperTriangularMatrix(const Matrix&)`.
- Methods: base overrides; `transpose() const`, `solve(const AbstractMatrix& b) const`, `multiplyUpper(const UpperTriangularMatrix&) const`.

`LowerTriangularMatrix`

- Constructors: `LowerTriangularMatrix(int size = 0)`, copy/move constructors, `LowerTriangularMatrix(const Matrix&)`.
- Methods: base overrides; `transpose() const`, `solve(const AbstractMatrix& b) const`, `multiplyLower(const LowerTriangularMatrix&) const`.

`TridiagonalMatrix`

- Constructors: `TridiagonalMatrix(int n = 0)`, `TridiagonalMatrix(const std::vector<double>& lower,const std::vector<double>& diag,const std::vector<double>& upper)`, `explicit TridiagonalMatrix(const Matrix&)`.
- Methods: base overrides; `solve(const AbstractMatrix& B) const`, `lower() const`, `diag() const`, `upper() const`.
- Contract: Thomas/chase method with no pivoting; zero pivots throw. Since M0.1, zero-size matrices are allowed and keep empty bands.

## `MatCal::Algorithm::Matrix`

- Elementary transformations: `swapRows`, `scaleRow`, `addScaledRow`, `swapCols`, `scaleCol`, `addScaledCol`, `applyRowSwaps`.
- Norms: `norm_one`, `norm_infinite`, `norm_Frobenius`.
- Direct solvers: `solve_columnElimination(AbstractMatrix& A, AbstractMatrix& b)`, `columnElimination_Transformation(AbstractMatrix& A)`, `determinant(AbstractMatrix& A)`, `LU_Decompose(AbstractMatrix& A,double eps = 1e-12)`.
- Result types: `LUresult { LowerTriangularMatrix L; UpperTriangularMatrix U; solve(AbstractMatrix& x); }`.
- Iterative solvers: `Iteration_Result`, `Jacobi`, `Gauss_Seidel`, `SOR(AbstractMatrix&,AbstractMatrix&,double omega,double epsilon = 1e-6,int max_iterations = 100)`, legacy `SOR(AbstractMatrix&,AbstractMatrix&,int omega,double epsilon = 1e-6,int max_iterations = 100)`.
- Eigen solvers: `PowerMethod_Result`, `PowerMethod`, `PowerMethod_reverse`.
- Input mutation: signatures use non-const references. Some functions copy internally and do not mutate input; elementary transformations do mutate input.

## `MatCal::Algorithm::Insert`

- `LagrangeInsert`: constructors from degree, vector pairs, initializer-list; `reconstruct`, `calculate`, `getDegree`, `getPoly`.
- `NewtonInsert_Quotient`: constructors, `reconstruct`, `calculate`, `insertNewTerm`, `getPoly`, `getSheet`, `getXs`, `getDegree`.
- `NewtonInsert_Finite`: constructors, `calculate`, `insertNewTerm`, `reSet_X0`, `getPoly`, `getSheet`, `getDegree`, `getH`, `getX0`.
- `Hermite`: constructors, `reconstruct`, `calculate`, `getDegree`, `getPoly`.
- `CubicSpline`: constructors, `reconstruct`, `calculate`, `getXs`, `getYs`, `getM`.
- `LinearInsert`: constructors from pair list or `xs/ys`, `reconstruct`, `calculate`, `getXs`, `getYs`.
- M6 note: all listed interpolation classes now keep these legacy signatures while delegating their core construction/evaluation formulas to `MatCal::Interpolation`.

## `MatCal::Algorithm::Iteration`

- `Bisection`: `solve`, `solveDetailed`, `Result`.
- `Picard`: `solve`, `solveDetailed`, `solve_Aitken`, `Result`.
- `Newton`: `solve`, `solve_downhill`, `Result`.
- `Secant`: `solve_two_point`, `solve_one_point`, `Result`.
- `NewtonForEquations`: `using Function = std::function<double(const std::vector<double>&)>`, `Result`, `solve`.
- M6 note: scalar classes delegate to `MatCal::Roots`; `NewtonForEquations::solve` now delegates to `MatCal::Nonlinear`.

## `MatCal::Algorithm::Basics`

- Constants: `PI`, `E`.
- `solve_Linear_System(Matrix& A, Matrix& b)`.
- `Derivative`: `dy_dx`, `dy_dx_center`, `pF_px`, `dF_dx`.
- `Least_Square`: nested `Result_least_square`, overloads of `solve` for weighted/unweighted and selected terms.
- `OrthogonalPolynomials`: `Chebyshev`, `ChebyshevZeros`, `Legendre`.
- `NumericalIntegration`: `Instant`, `NewtonCotes`, `CompositeNewtonCotes`, `Romberg`, overloads for discrete data, `CotesSheet`.
- `ODE`: `SimpleEuler`, `Euler`, `RungeKutta_44`.
- `Integrate::RK4`: `step`, `step2`; this area is explicitly marked in source as PT-specific legacy support.
- M6 note: `Least_Square::solve` overloads now delegate to `MatCal::LeastSquares`. Multivariable derivative helpers `pF_px` and `dF_dx` remain legacy-only.

## Possible PT/External Dependency Surface

Treat these as likely used until proven otherwise:

- All public headers and namespaces listed above.
- `MatCal::Algorithm::Basics::Integrate::RK4::step` and `step2`, due source comments and commit history.
- `CubicSpline` and `LinearInsert`, due PT mention in history.
- `Matrix`, `QinJiuShao`, `solve_columnElimination`, `Gauss_Seidel`, and interpolation classes, because README presents them as public API.

## Signature Reference

This compact reference preserves return types, parameter types, and default parameters for the most dependency-sensitive public calls.

### Matrix Algorithms

```cpp
void swapRows(AbstractMatrix& A, int r1, int r2);
void scaleRow(AbstractMatrix& A, int r, double scalar);
void addScaledRow(AbstractMatrix& A, int r_target, int r_source, double scalar);
void swapCols(AbstractMatrix& A, int c1, int c2);
void scaleCol(AbstractMatrix& A, int c, double scalar);
void addScaledCol(AbstractMatrix& A, int c_target, int c_source, double scalar);
void applyRowSwaps(AbstractMatrix& M, const std::vector<std::pair<int,int>>& swaps, bool reverse = false);
double norm_one(AbstractMatrix& matrix);
double norm_infinite(AbstractMatrix& matrix);
double norm_Frobenius(AbstractMatrix& matrix);
std::unique_ptr<AbstractMatrix> solve_columnElimination(AbstractMatrix& A, AbstractMatrix& b);
std::pair<std::unique_ptr<AbstractMatrix>, std::vector<std::pair<int,int>>> columnElimination_Transformation(AbstractMatrix& A);
double determinant(AbstractMatrix& A);
LUresult LU_Decompose(AbstractMatrix& A, double eps = 1e-12);
Iteration_Result Jacobi(AbstractMatrix& A, AbstractMatrix& b, double epsilon = 1e-6, int max_iterations = 100);
Iteration_Result Gauss_Seidel(AbstractMatrix& A, AbstractMatrix& b, double epsilon = 1e-6, int max_iterations = 100);
Iteration_Result SOR(AbstractMatrix& A, AbstractMatrix& b, double omega, double epsilon = 1e-6, int max_iterations = 100);
Iteration_Result SOR(AbstractMatrix& A, AbstractMatrix& b, int omega, double epsilon = 1e-6, int max_iterations = 100);
PowerMethod_Result PowerMethod(AbstractMatrix& A, double eps = 1e-6, int max_iter = 1000);
PowerMethod_Result PowerMethod_reverse(AbstractMatrix& A, double near_num = 0, double eps = 1e-6, int max_iter = 1000);
```

### Basics

```cpp
std::pair<Matrix, std::string> solve_Linear_System(Matrix& A, Matrix& b);
double Derivative::dy_dx(Func_y f, double x, double eps = 1e-6);
double Derivative::dy_dx_center(Func_y f, double x, double eps = 1e-6);
double Derivative::pF_px(Func_F f, std::vector<double>& xs, int i, double eps = 1e-6);
double Derivative::dF_dx(Func_F f, std::vector<double>& xs, double eps = 1e-6);
Least_Square::Result_least_square Least_Square::solve(int degree, std::vector<double>& x, std::vector<double>& y, std::vector<double>& weights);
Least_Square::Result_least_square Least_Square::solve(int degree, std::vector<double>& x, std::vector<double>& y);
Least_Square::Result_least_square Least_Square::solve(int degree, std::vector<double>& x, std::vector<double>& y, std::vector<double>& weights, std::vector<bool>& selects);
Least_Square::Result_least_square Least_Square::solve(int degree, std::vector<double>& x, std::vector<double>& y, std::vector<bool>& selects);
QinJiuShao OrthogonalPolynomials::Chebyshev(int n, bool second = false);
std::vector<double> OrthogonalPolynomials::ChebyshevZeros(int n);
QinJiuShao OrthogonalPolynomials::Legendre(int n);
double NumericalIntegration::Instant(std::function<double(double)> f, double a, double b, double eps = 1e-6);
double NumericalIntegration::NewtonCotes(std::function<double(double)> f, double a, double b, int n = 4);
double NumericalIntegration::CompositeNewtonCotes(std::function<double(double)> f, double a, double b, int segments = 10, int n = 4);
std::pair<double, Matrix> NumericalIntegration::Romberg(std::function<double(double)> f, double a, double b, double eps = 1e-6, int maxIterations = 20);
Matrix ODE::SimpleEuler(int n, std::vector<std::function<double(std::vector<double>&)>>& funcs, std::vector<double>& inits, double h = 1e-2, int count = 100);
std::pair<Matrix, Matrix> ODE::Euler(int n, std::vector<std::function<double(std::vector<double>&)>>& funcs, std::vector<double>& inits, double h = 1e-2, int count = 100);
Matrix ODE::RungeKutta_44(int n, std::vector<std::function<double(std::vector<double>&)>>& funcs, std::vector<double>& inits, double h = 1e-2, int count = 100);
void Integrate::RK4::step(const RHS& f, const std::vector<double>& y, double dt, std::vector<double>& y_out);
void Integrate::RK4::step2(const RHS2& f, double theta, double omega, double dt, double& theta_out, double& omega_out);
```

M5 adds independent 0.x targets while preserving every legacy signature above:

```cpp
// MatCal::Calculus
DerivativeResult forward_difference(Function f, double x, double step);
DerivativeResult central_difference(Function f, double x, double step);
IntegrationResult integrate_instant(Function f, double a, double b, double step);
IntegrationResult integrate_newton_cotes(Function f, double a, double b, int order);
IntegrationResult integrate_composite_newton_cotes(Function f, double a, double b, int segments, int order);
IntegrationResult integrate_romberg(Function f, double a, double b, const IntegrationOptions& options = {});

// MatCal::ODE
OdeStepResult euler_step(Rhs rhs, double t, const std::vector<double>& state, double dt, const OdeOptions& options = {});
std::pair<OdeStepResult, std::vector<double>> improved_euler_step(Rhs rhs, double t, const std::vector<double>& state, double dt, const OdeOptions& options = {});
OdeStepResult rk4_step(Rhs rhs, double t, const std::vector<double>& state, double dt, const OdeOptions& options = {});
OdeTrajectoryResult integrate_euler(Rhs rhs, double t0, const std::vector<double>& state0, double dt, std::size_t count, const OdeOptions& options = {});
OdeTrajectoryResult integrate_rk4(Rhs rhs, double t0, const std::vector<double>& state0, double dt, std::size_t count, const OdeOptions& options = {});
```

### Iteration

```cpp
double Bisection::solve(Function f, double a, double b, double epsilon = 1e-6, int maxIterations = 1000);
Bisection::Result Bisection::solveDetailed(Function f, double a, double b, double epsilon = 1e-6, int maxIterations = 1000);
double Picard::solve(Function function_varphi, double x_0, double epsilon = 1e-6, int maxIterations = 1000);
Picard::Result Picard::solveDetailed(Function function_varphi, double x_0, double epsilon = 1e-6, int maxIterations = 1000);
Picard::Result Picard::solve_Aitken(Function function_varphi, double x_0, double epsilon = 1e-6, int maxIterations = 1000);
Newton::Result Newton::solve(Function f, Function df_dx, double x_0, double epsilon = 1e-6, int maxIterations = 1000);
Newton::Result Newton::solve_downhill(Function f, Function df_dx, double x_0, double epsilon = 1e-6, int maxIterations = 1000);
Secant::Result Secant::solve_two_point(Function f, double x0, double x1, double epsilon = 1e-6, int maxIterations = 1000);
Secant::Result Secant::solve_one_point(Function f, double x0, double h = 1e-4, double epsilon = 1e-6, int maxIterations = 1000);
NewtonForEquations::Result NewtonForEquations::solve(int n, std::vector<Function>& Funcs, std::vector<double> xs, int maxIterations = 20, double epsilon = 1e-6);
```
