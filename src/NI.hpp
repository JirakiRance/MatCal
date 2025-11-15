#ifndef NI_HPP
#define NI_HPP

#include "QinJiuShao.hpp"
#include "Matrix.hpp"
#include <vector>
#include <stdexcept>

namespace MatCal {
    namespace Utils {
        //秦九韶结点，用于秦九韶多项式
        class QinJiuShaoNode;
        //秦九韶多项式，提供计算函数值，返回多项式函数(创建通用函数)等服务
        class QinJiuShao;
        //Newton差值以及差分,基于秦九韶、矩阵实现
        class NewtonInsert_Quotient;
        class NewtonInsert_Finite;
    }
}

namespace MatCal::Utils {

//------------------------------------------------------------
// Newton 差商插值（不等距点）
//------------------------------------------------------------
class NewtonInsert_Quotient {
    //外部using声明
    using Poly  = MatCal::Utils::QinJiuShao;
    using Lower = MatCal::Utils::LowerTriangularMatrix;
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
//------------------------------------------------------------
class NewtonInsert_Finite {
    //外部using声明
    using Poly  = MatCal::Utils::QinJiuShao;
    using Lower = MatCal::Utils::LowerTriangularMatrix;

//公共方法
public:
    NewtonInsert_Finite(double h):h(h){
        if(h <= QinJiuShaoNode::ZERO_THRESHOLD)
            throw std::invalid_argument("h must be > 0 for finite difference");
    }
    NewtonInsert_Finite(double h, std::vector<std::pair<double,double>>& data): h(h){
        if(h <= QinJiuShaoNode::ZERO_THRESHOLD)
            throw std::invalid_argument("h must be > 0 for finite difference");
        construct(data);
    }
    double calculate(double X) const { 
        return poly.calculate(X);
    }
    //插入新 Y（自动取 x = x0 + degree*h）
    void insertNewTerm(double Y){
        int newSize = degree + 2;
        lower.resize(newSize,newSize);
        int r = newSize - 1;
        //第一列：新数据点
        lower.set(r,0,Y);
        //差分表更新
        for(int j=1;j<newSize;j++)
            lower.set(r,j, lower.get(r,j-1) - lower.get(r-1,j-1));
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

    const Poly&  getPoly()  const { return poly; }
    const Lower& getSheet() const { return lower; }
    int          getDegree()const { return degree; }
    double       getH()     const { return h; }
    double       getX0()    const { return x0; }

private:

    //构造
    void construct(std::vector<std::pair<double,double>>& data){
        int n = data.size();
        if(n < 2) throw std::invalid_argument("Need >=2 points");

        //起点 x0
        x0 = data[0].first;
        degree = n - 1;

        //差分表初始化
        lower.resize(n,n);

        //第一列
        for(int i=0;i<n;i++)
            lower.set(i,0,data[i].second);

        //差分
        for(int j=1;j<n;j++)
            for(int i=j;i<n;i++)
                lower.set(i,j, lower.get(i,j-1) - lower.get(i-1,j-1));

        //构造多项式
        Poly res = { {0, lower.get(0,0)} };

        for(int i=1;i<n;i++){
            double coe = lower.get(i,i);

            //按你的旧写法：denom = ∏(h*k)
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

} // namespace MatCal::Utils

#endif
