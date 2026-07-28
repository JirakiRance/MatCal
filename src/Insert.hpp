#ifndef INSERT_HPP
#define INSERT_HPP

namespace MatCal {
    namespace Utils {
        //秦九韶结点，用于秦九韶多项式
        class QinJiuShaoNode;
        //秦九韶多项式，提供计算函数值，返回多项式函数(创建通用函数)等服务
        class QinJiuShao;
    }
    namespace Algorithm{
        namespace Insert{
            //Lagrange插值多项式Ln(x),基于秦九韶实现
            class LagrangeInsert;
            //Newton差商插值,基于秦九韶、矩阵实现
            class NewtonInsert_Quotient;
            //Newton差分插值,基于秦九韶、矩阵实现
            class NewtonInsert_Finite;
            //Hermite插值，适用于有导数的数据表
            class Hermite;
        }
    }
}

#include<cstddef>
#include<limits>

#include"QinJiuShao.hpp"
#include"Matrix.hpp"
#include"MatCal/Interpolation/CubicSpline.hpp"
#include"MatCal/Interpolation/LinearInterpolator.hpp"
#include"MatCal/Interpolation/PolynomialInterpolation.hpp"


namespace MatCal::Algorithm::Insert{

namespace detail {

inline MatCal::Utils::QinJiuShao to_legacy_polynomial(const MatCal::Polynomial::Polynomial& polynomial) {
    std::vector<std::pair<int, double>> terms;
    for (const auto& term : polynomial.terms_descending()) {
        if (term.first > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
            throw std::length_error("legacy polynomial degree exceeds int range");
        }
        terms.emplace_back(static_cast<int>(term.first), term.second);
    }
    if (terms.empty()) {
        terms.emplace_back(0, 0.0);
    }
    return MatCal::Utils::QinJiuShao(terms);
}

inline void fill_lower_sheet(MatCal::Utils::LowerTriangularMatrix& lower,
                             const std::vector<std::vector<double>>& table) {
    if (table.empty()) {
        lower.resize(0, 0);
        return;
    }
    if (table.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::length_error("legacy interpolation sheet exceeds int range");
    }
    const int n = static_cast<int>(table.size());
    lower.resize(n, n);
    for (int row = 0; row < n; ++row) {
        for (int col = 0; col <= row; ++col) {
            lower.set(row, col, table[static_cast<std::size_t>(row)][static_cast<std::size_t>(col)]);
        }
    }
}

} // namespace detail

//Lagrange插值多项式Ln(x),基于秦九韶实现
class LagrangeInsert {
    using Poly=MatCal::Utils::QinJiuShao;
//构造和析构
public:
    //默认构造
    LagrangeInsert(int degree_of_poly = 2) {
        if (degree_of_poly < 2) {
            this->degree = 2;
        }
        else {
            this->degree = degree_of_poly;
        }
    }
    //从std::vector<std::pair<double,double>>指定的x,y序列构造
    LagrangeInsert(std::vector<std::pair<double, double>>& data) {
        this->construct_from_vector(data);
    }
    //从std::initializer_list<std::pair<double,double>>指定的x,y序列构造
    LagrangeInsert(std::initializer_list<std::pair<double,double>> pairs) {
        std::vector<std::pair<double, double>> data = pairs;
        this->construct_from_vector(data);
    }
    //析构，不用管
    ~LagrangeInsert() {}
//公共方法
public:
    //外部构造接口，提供给默认构造的对象
    void reconstruct(std::vector<std::pair<double, double>>& data){
        this->construct_from_vector(data);
    }
    //计算插值
    double calculate(double x) {
        return this->poly.calculate(x);
    }
    //getter方法
    int getDegree()const {
        return this->degree;
    }
    const Poly getPoly() const{
        return this->poly;
    }

//内部维护的方法
private:
    void construct_from_vector(std::vector<std::pair<double, double>>& data) {
        if (data.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
            throw std::length_error("legacy Lagrange degree exceeds int range");
        }
        auto polynomial = MatCal::Interpolation::interpolate_lagrange(data);
        this->degree = static_cast<int>(data.size()) - 1;
        this->poly = detail::to_legacy_polynomial(polynomial);
    }

//内部维护的属性
private:
    int degree;   //Ln(x) 次数n
    Poly poly;    //多项式结构，由秦九韶维护
};//class LagrangeInsert


//------------------------------------------------------------
// Newton 差商插值（不等距点）
//------------------------------------------------------------
class NewtonInsert_Quotient {
    //外部using声明
    using Poly  = MatCal::Utils::QinJiuShao;
    using Lower = MatCal::Utils::LowerTriangularMatrix;
    using QinJiuShaoNode=MatCal::Utils::QinJiuShaoNode;
//公共方法
public:
    NewtonInsert_Quotient() = default;
    NewtonInsert_Quotient(std::vector<std::pair<double,double>>& data){
        construct(data);
    }
    NewtonInsert_Quotient(const std::initializer_list<std::pair<double,double>>& data){
        std::vector<std::pair<double,double>> v = data;
        construct(v);
    }
    void reconstruct(std::vector<std::pair<double,double>>& data){
        x.clear();
        lower.resize(0,0);
        poly = Poly();
        construct(data);
    }
    double calculate(double X) const {
        return poly.calculate(X); 
    }
    void insertNewTerm(double X, double Y){
        if (x.empty())
            throw std::out_of_range("insertNewTerm: no existing nodes — call construct() with >=2 points first");
        std::vector<std::pair<double, double>> data;
        data.reserve(x.size() + 1);
        for (std::size_t i = 0; i < x.size(); ++i) {
            data.emplace_back(x[i], lower.get(static_cast<int>(i), 0));
        }
        data.emplace_back(X, Y);
        construct(data);
    }

    //getter方法
    const Poly& getPoly() const { return poly; }
    const Lower& getSheet() const { return lower; }
    const std::vector<double>& getXs() const { return x; }
    int getDegree() const { return degree; }
//内部维护的方法
private:
    void construct(std::vector<std::pair<double,double>>& data){
        if (data.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
            throw std::length_error("legacy Newton divided-difference degree exceeds int range");
        }
        auto result = MatCal::Interpolation::interpolate_newton_divided(data);
        x = result.xs;
        detail::fill_lower_sheet(lower, result.table);
        degree = static_cast<int>(data.size()) - 1;
        poly = detail::to_legacy_polynomial(result.polynomial);
    }
//内部维护的属性
private:
    Poly poly;      //插值多项式
    Lower lower;    //存储表
    std::vector<double> x;  //零点
    int degree = 0;         //次数
};

//------------------------------------------------------------
// Newton 差分插值（等距点，差分 h 固定）
// 版本改动：construct() 只接收 Y 序列，不再需要 x
// x 自动为 x0 + i*h
//------------------------------------------------------------
class NewtonInsert_Finite {
    //外部using声明
    using Poly  = MatCal::Utils::QinJiuShao;
    using Lower = MatCal::Utils::LowerTriangularMatrix;
    using QinJiuShaoNode=MatCal::Utils::QinJiuShaoNode;

//公共方法
public:

    //给定步长 h
    NewtonInsert_Finite(double h):h(h){
        if(h <= QinJiuShaoNode::ZERO_THRESHOLD)
            throw std::invalid_argument("h must be > 0 for finite difference");
    }

    //给定 h 和 y 序列
    NewtonInsert_Finite(double h, double x0, std::vector<double>& yData): h(h), x0(x0){
        if(h <= QinJiuShaoNode::ZERO_THRESHOLD)
            throw std::invalid_argument("h must be > 0 for finite difference");
        construct(yData);
    }

    //求值
    double calculate(double X) const { 
        return poly.calculate(X);
    }

    //插入新 Y（自动使用 x = x0 + degree*h）
    void insertNewTerm(double Y){
        if (degree < 1)
            throw std::out_of_range("insertNewTerm: insufficient base data — call construct(yData) with >=2 values first");
        std::vector<double> yData;
        yData.reserve(static_cast<std::size_t>(degree) + 2);
        for (int i = 0; i <= degree; ++i) {
            yData.push_back(lower.get(i, 0));
        }
        yData.push_back(Y);
        construct(yData);
    }

    void reSet_X0(int x0){
        this->x0=x0;
    }

    const Poly&  getPoly()  const { return poly; }
    const Lower& getSheet() const { return lower; }
    int          getDegree()const { return degree; }
    double       getH()     const { return h; }
    double       getX0()    const { return x0; }

private:

    //------------------------------------------------------------
    // 构造：只传 Y，x = x0 + i*h
    //------------------------------------------------------------
    void construct(std::vector<double>& yData){
        if (yData.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
            throw std::length_error("legacy Newton finite-difference degree exceeds int range");
        }
        auto result = MatCal::Interpolation::interpolate_newton_finite(h, x0, yData);
        detail::fill_lower_sheet(lower, result.table);
        degree = static_cast<int>(yData.size()) - 1;
        poly = detail::to_legacy_polynomial(result.polynomial);
    }

private:
    double h;       //步长 h
    double x0 = 0;  //起始点 x0
    Poly  poly;     //插值多项式
    Lower lower;    //差分表
    int degree = 0; //阶数
};


class Hermite{
public:
    using Poly=MatCal::Utils::QinJiuShao;
    using PolyNode=MatCal::Utils::QinJiuShaoNode;
//构造和析构
public:
    Hermite(int degree_of_poly = 2) {
        if (degree_of_poly < 2) {
            this->degree = 2;
        }
        else {
            this->degree = degree_of_poly;
        }
    }
    //从std::vector<std::pair<double,double>>指定的x,y序列构造
    Hermite(std::vector<double>&xs,std::vector<double>&ys,std::vector<double>&dy_dxs) {
        this->construct(xs,ys,dy_dxs);
    }
    //析构，不用管
    ~Hermite() {}
//公共的方法
public:
    //重新构造（给默认构造用）
    void reconstruct(std::vector<double>& xs,
                     std::vector<double>& ys,
                     std::vector<double>& dy_dxs) {
        construct(xs, ys, dy_dxs);
    }

     //计算插值
    double calculate(double x) {
        return this->poly.calculate(x);
    }
    //getter方法
    int getDegree()const {
        return this->degree;
    }
    const Poly getPoly() const{
        return this->poly;
    }

//内部维护的方法
private:
void construct(std::vector<double>& xs,
                   std::vector<double>& ys,
                   std::vector<double>& dy_dxs) {
        if (xs.size() > (static_cast<std::size_t>(std::numeric_limits<int>::max()) + 1) / 2) {
            throw std::length_error("legacy Hermite degree exceeds int range");
        }
        auto polynomial = MatCal::Interpolation::interpolate_hermite(xs, ys, dy_dxs);
        this->degree = static_cast<int>(2 * xs.size()) - 1;
        this->poly = detail::to_legacy_polynomial(polynomial);
    }

//内部维护的属性
private:
    int degree;   //Ln(x) 次数n
    Poly poly;    //多项式结构，由秦九韶维护
};


//==============================================================
//  Cubic Spline 插值（默认：自然边界样条 Natural Spline）
//==============================================================
class CubicSpline {
public:
    using Matrix = MatCal::Utils::Matrix;
    using Poly   = MatCal::Utils::QinJiuShao;

//--------------------------------------------------------------
// 构造与析构
//--------------------------------------------------------------
public:
    CubicSpline() = default;

    CubicSpline(const std::vector<double>& xs,
                const std::vector<double>& ys)
    {
        construct(xs, ys);
    }

    ~CubicSpline() = default;

//--------------------------------------------------------------
// 公共接口
//--------------------------------------------------------------
public:
    void reconstruct(const std::vector<double>& xs,
                     const std::vector<double>& ys)
    {
        construct(xs, ys);
    }

    // 计算样条
    double calculate(double x) const
    {
        return spline.evaluate(x);
    }

    const std::vector<double>& getXs() const { return X; }
    const std::vector<double>& getYs() const { return Y; }
    const std::vector<double>& getM()  const { return M; }

//--------------------------------------------------------------
// 内部构造
//--------------------------------------------------------------
private:
    void construct(const std::vector<double>& xs,
                   const std::vector<double>& ys)
    {
        spline = MatCal::Interpolation::CubicSpline(xs, ys, MatCal::Interpolation::ExtrapolationPolicy::extrapolate);
        X = spline.xs();
        Y = spline.ys();
        M = spline.second_derivatives();
    }

//--------------------------------------------------------------
// 区间查找（二分）
//--------------------------------------------------------------
    int findInterval(double x) const
    {
        int n = X.size();
        if (x <= X[0]) return 0;
        if (x >= X[n-1]) return n-2;

        int L = 0, R = n-1;
        while (L <= R) {
            int mid = (L + R) / 2;
            if (X[mid] <= x && x <= X[mid+1])
                return mid;

            if (X[mid] < x) L = mid + 1;
            else R = mid - 1;
        }
        return n-2;
    }

//--------------------------------------------------------------
// 内部数据
//--------------------------------------------------------------
private:
    std::vector<double> X;   // 节点 x_i
    std::vector<double> Y;   // 节点 y_i
    std::vector<double> M;   // 二阶导数 M_i
    MatCal::Interpolation::CubicSpline spline;
};


//------------------------------------------------------------
// Linear 插值（分段线性插值）
//------------------------------------------------------------
class LinearInsert {
public:
    // 构造与析构
    LinearInsert() = default;

    LinearInsert(const std::vector<std::pair<double,double>>& data) {
        construct_from_pairs(data);
    }

    LinearInsert(const std::initializer_list<std::pair<double,double>>& data) {
        std::vector<std::pair<double,double>> v = data;
        construct_from_pairs(v);
    }

    LinearInsert(const std::vector<double>& xs,
                 const std::vector<double>& ys)
    {
        construct(xs, ys);
    }

    ~LinearInsert() = default;

public:
    //--------------------------------------------------------
    // 重新构造（与其它插值类保持一致）
    //--------------------------------------------------------
    void reconstruct(const std::vector<double>& xs,
                     const std::vector<double>& ys)
    {
        construct(xs, ys);
    }

    void reconstruct(const std::vector<std::pair<double,double>>& data)
    {
        construct_from_pairs(data);
    }

    //--------------------------------------------------------
    // 插值计算
    //--------------------------------------------------------
    double calculate(double x) const {
        return interpolator.evaluate(x);
    }

    //--------------------------------------------------------
    // Getter
    //--------------------------------------------------------
    const std::vector<double>& getXs() const { return X; }
    const std::vector<double>& getYs() const { return Y; }

private:
    //--------------------------------------------------------
    // 使用 pair 形式构造
    //--------------------------------------------------------
    void construct_from_pairs(const std::vector<std::pair<double,double>>& data)
    {
        interpolator = MatCal::Interpolation::LinearInterpolator(data, MatCal::Interpolation::ExtrapolationPolicy::extrapolate);
        X = interpolator.xs();
        Y = interpolator.ys();
    }

    //--------------------------------------------------------
    // 标准构造
    //--------------------------------------------------------
    void construct(const std::vector<double>& xs,
                   const std::vector<double>& ys)
    {
        interpolator = MatCal::Interpolation::LinearInterpolator(xs, ys, MatCal::Interpolation::ExtrapolationPolicy::extrapolate);
        X = interpolator.xs();
        Y = interpolator.ys();
    }

    //--------------------------------------------------------
    // 校验并保证 X 单调递增（你其它插值函数也做类似检查）
    //--------------------------------------------------------
    void check_and_sort() {
        int n = X.size();
        for (int i = 1; i < n; ++i) {
            if (X[i] <= X[i-1]) {
                throw std::invalid_argument(
                    "LinearInsert: x values must be strictly increasing.");
            }
        }
    }

    //--------------------------------------------------------
    // 二分查找区间 i，使 x ∈ [X[i], X[i+1]]
    //--------------------------------------------------------
    int findInterval(double x) const {
        int L = 0, R = static_cast<int>(X.size()) - 2;
        while (L <= R) {
            int mid = (L + R) / 2;
            if (X[mid] <= x && x <= X[mid+1])
                return mid;
            if (X[mid] < x)
                L = mid + 1;
            else
                R = mid - 1;
        }
        // 理论上不会走到这里
        return std::max(0, L);
    }

    //--------------------------------------------------------
    // 两点间线性插值
    //--------------------------------------------------------
    double interpolate(int i, int j, double x) const {
        double x0 = X[i], x1 = X[j];
        double y0 = Y[i], y1 = Y[j];
        double t = (x - x0) / (x1 - x0);
        return y0 + t * (y1 - y0);
    }

private:
    std::vector<double> X;
    std::vector<double> Y;
    MatCal::Interpolation::LinearInterpolator interpolator;
};



}//namespace MatCal::Algorithm::Insert
#endif//INSERT_HPP
