#ifndef LAGRANGEINSERT_HPP
#define LAGRANGEINSERT_HPP


namespace MatCal {
    namespace Utils {
        //秦九韶结点，用于秦九韶多项式
        class QinJiuShaoNode;
        //秦九韶多项式，提供计算函数值，返回多项式函数(创建通用函数)等服务
        class QinJiuShao;
        //Lagrange插值多项式Ln(x),基于秦九韶实现
        class LagrangeInsert;
    }
}
#include"QinJiuShao.hpp"

namespace MatCal::Utils {

//Lagrange插值多项式Ln(x),基于秦九韶实现
    class LagrangeInsert {
        using Poly=MatCal::Utils::QinJiuShao;
        using PolyNode=MatCal::Utils::QinJiuShaoNode;
    //构造和析构
    public:
        //默认构造
        LagrangeInsert(int degree_of_poly = 2) {
            if (degree_of_poly > 2) {
                this->degree = 2;
            }
            else {
                this->degree = 2;
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
        //计算插值
        double calculate(double x) {
            return this->poly.calculate(x);
        }
        //gettet方法
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

}//namespace MatCal::Utils
#endif//LAGRANGEINSERT_HPP