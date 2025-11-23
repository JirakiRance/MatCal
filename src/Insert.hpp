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

#include"QinJiuShao.hpp"
#include"Matrix.hpp"


namespace MatCal::Algorithm::Insert{

//Lagrange插值多项式Ln(x),基于秦九韶实现
class LagrangeInsert {
    using Poly=MatCal::Utils::QinJiuShao;
    using PolyNode=MatCal::Utils::QinJiuShaoNode;
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
    const int getDegree()const {
        return this->degree;
    }
    const Poly getPoly() const{
        return this->poly;
    }

//内部维护的方法
private:
    void construct_from_vector(std::vector<std::pair<double, double>>& data) {
        if (data.size() < 2) {
            throw std::invalid_argument("poly terms should at least >= 2!");
        }
        this->degree = data.size() - 1;
        //创建n+1个l_i(x)再相加
        std::vector<Poly> li_x(this->degree+1);
        for (int i = 0; i <= this->degree; ++i) {//li
            //创建系数
            double coee = data[i].second;
            //创建主式(数字1)
            Poly li({
                { 0, 1 }
                });
            for (int j = 0; j <= this->degree; ++j) {
                if (i == j)
                    continue;
                //除以系数
                coee /= (data[i].first - data[j].first);
                //乘以上式
                li = li * Poly({
                    { 1 , 1 }, { 0 , - data[j].first }
                    });
            }
            //最后乘以系数
            li = li * coee;
            //加入li_x
            li_x[i] = li;
        }//li

        //最后加起来
        Poly Ln_x = *li_x.begin();
        for (int i = 1; i <= this->degree; ++i) {
            Ln_x = Ln_x + li_x[i];
        }
        this->poly = Ln_x;
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
        if (X == x.back())
            throw std::invalid_argument("insertNewTerm: duplicate x value");
        degree++;
        int newSize = degree + 1;
        x.push_back(X);
        lower.resize(newSize,newSize);

        int r = newSize - 1;
        lower.set(r,0,Y);
        lower.set(r,1,(Y - lower.get(r-1,0)) / (X - x[r-1]));

        for(int j=2;j<newSize;j++){
            lower.set(r,j,(lower.get(r,j-1) - lower.get(r-1,j-1)) / (x[r] - x[r-j]));
        }

        Poly pn(lower.get(r,r), x);
        poly = poly + pn;
    }

    //getter方法
    const Poly& getPoly() const { return poly; }
    const Lower& getSheet() const { return lower; }
    const std::vector<double>& getXs() const { return x; }
    int getDegree() const { return degree; }
//内部维护的方法
private:
    void construct(std::vector<std::pair<double,double>>& data){
        int n=data.size();
        if(n<2) throw std::invalid_argument("Need >=2 points");

        degree = n-1;
        x.resize(n);
        for(int i=0;i<n;i++) x[i]=data[i].first;

        lower.resize(n,n);
        for(int j=0;j<n;j++) lower.set(j,0,data[j].second);
        for(int j=1;j<n;j++) lower.set(j,1,(lower.get(j,0)-lower.get(j-1,0))/(x[j]-x[j-1]));
        for(int j=2;j<n;j++)
            for(int i=j;i<n;i++)
                lower.set(i,j,(lower.get(i,j-1)-lower.get(i-1,j-1))/(x[i]-x[i-j]));

        Poly res = {{0,lower.get(0,0)}};
        for(int i=1;i<n;i++){
            std::vector<double> sub(x.begin(),x.begin()+i);
            Poly pn(lower.get(i,i),sub);
            res = res + pn;
        }
        poly = res;
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
        int newSize = degree + 2;
        lower.resize(newSize,newSize);
        int r = newSize - 1;
        //第一列：新数据点
        lower.set(r,0,Y);
        //差分表更新 Δ, Δ², ...
        for(int j=1;j<newSize;j++)
            lower.set(r,j, lower.get(r,j-1) - lower.get(r-1,j-1));
        //对应系数
        double coe = lower.get(r,r);
        double denom = 1;
        for(int k=1;k<=r;k++) denom *= (h*k);
        coe /= denom;
        if(std::abs(coe)<QinJiuShaoNode::ZERO_THRESHOLD)
            coe = QinJiuShaoNode::ZERO_THRESHOLD*1.00001;
        //构造基函数 (X - x0)(X - (x0+h))...(长度 r)
        Poly pn = makeBasis(r, coe);
        //加入多项式
        poly = poly + pn;
        degree++;
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
        int n = yData.size();
        if(n < 2) throw std::invalid_argument("Need >=2 points");
        degree = n - 1;
        //差分表初始化
        lower.resize(n,n);
        //第一列 = y 值
        for(int i=0;i<n;i++)
            lower.set(i,0,yData[i]);
        //构造前向差分
        for(int j=1;j<n;j++)
            for(int i=j;i<n;i++)
                lower.set(i,j, lower.get(i,j-1) - lower.get(i-1,j-1));
        //--------------------------------------------------------
        //构造多项式
        //--------------------------------------------------------
        Poly res = { {0, lower.get(0,0)} };
        for(int i=1;i<n;i++){
            double coe = lower.get(i,i);
            double denom = 1;
            for(int k=1;k<=i;k++) denom *= (h*k);
            coe /= denom;
            if(std::abs(coe)<QinJiuShaoNode::ZERO_THRESHOLD)
                coe = QinJiuShaoNode::ZERO_THRESHOLD*1.00001;
            Poly pn = makeBasis(i, coe);
            res = res + pn;
        }
        poly = res;
    }

    //------------------------------------------------------------
    // 生成 Newton 等距节点基函数：
    // (X - x0)(X - (x0+h))...(X - (x0+(n-1)h)) * coe
    //------------------------------------------------------------
    Poly makeBasis(int n, double coe) const {
        std::vector<double> xs;
        xs.reserve(n);
        for(int i=0;i<n;i++)
            xs.push_back(x0 + i*h);
        return Poly(coe, xs);
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
    const int getDegree()const {
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
        int n = static_cast<int>(xs.size());
        if (n < 2 || ys.size() != xs.size() || dy_dxs.size() != xs.size()) {
            throw std::invalid_argument(
                "Hermite: vector size should be >=2 "
                "and xs, ys, dy_dxs size must match!");
        }

        // 检查节点是否重复
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                if (std::abs(xs[i] - xs[j]) < PolyNode::ZERO_THRESHOLD) {
                    throw std::invalid_argument(
                        "Hermite: duplicated x nodes detected!");
                }
            }
        }

        // 理论次数：2n - 1
        this->degree = 2 * n - 1;

        // 1. 先构造标准拉格朗日基函数 l_i(x)
        std::vector<Poly> li_x(n);
        std::vector<double> li_prime_at_xi(n);  // l_i'(x_i)

        for (int i = 0; i < n; ++i) {
            // li(x) 从常数 1 开始
            Poly li({{0, 1.0}});
            double denom = 1.0;

            for (int j = 0; j < n; ++j) {
                if (i == j) continue;
                double diff = xs[i] - xs[j];
                if (std::abs(diff) < PolyNode::ZERO_THRESHOLD) {
                    throw std::invalid_argument(
                        "Hermite: xs[i] - xs[j] too small, "
                        "nodes too close or duplicated.");
                }
                denom *= diff;
                // 乘 (x - x_j)
                li = li * Poly({{1, 1.0}, {0, -xs[j]}});
            }

            // li(x) = ∏(x - x_j) / ∏(x_i - x_j)
            li = li * (1.0 / denom);
            li_x[i] = li;

            // 求导并在 x_i 处评价： l_i'(x_i)
            Poly dli = li.derivative();
            li_prime_at_xi[i] = dli.calculate(xs[i]);
        }

        // 2. 用 l_i(x) 和 l_i'(x_i) 构造 Hermite 基函数
        //    H_i(x) = (1 - 2(x - x_i) l_i'(x_i)) [l_i(x)]^2
        //    \hat{H}_i(x) = (x - x_i)[l_i(x)]^2
        Poly result;        // 默认构造应为 0 多项式
        bool firstTerm = true;

        for (int i = 0; i < n; ++i) {
            Poly X_minus_xi({{1, 1.0}, {0, -xs[i]}});
            Poly li_sq = li_x[i] * li_x[i];

            // H_i(x)
            Poly Hi = Poly({{0, 1.0}});                     // 1
            Hi = Hi - (2.0 * li_prime_at_xi[i]) * X_minus_xi;  // 1 - 2 l_i'(x_i)(x - x_i)
            Hi = Hi * li_sq;

            // \hat{H}_i(x)
            Poly Htilde = X_minus_xi * li_sq;

            // y_i * H_i(x) + y'_i * \hat{H}_i(x)
            Poly term = Hi * ys[i] + Htilde * dy_dxs[i];

            if (firstTerm) {
                result = term;
                firstTerm = false;
            } else {
                result = result + term;
            }
        }

        this->poly = result;
    }

//内部维护的属性
private:
    int degree;   //Ln(x) 次数n
    Poly poly;    //多项式结构，由秦九韶维护
};

}//namespace MatCal::Algorithm::Insert
#endif//INSERT_HPP