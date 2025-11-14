#ifndef NEWTONINSERT_HPP
#define NEWTONINSERT_HPP

namespace MatCal {
    namespace Utils {
        //秦九韶结点，用于秦九韶多项式
        class QinJiuShaoNode;
        //秦九韶多项式，提供计算函数值，返回多项式函数(创建通用函数)等服务
        class QinJiuShao;
        //Newton差值以及差分,基于秦九韶、矩阵实现
        class NewtonInsert;
    }
}
#include"QinJiuShao.hpp"
#include"Matrix.hpp"


namespace MatCal::Utils{

enum DifferenceType{
    DIFFERENCE_QUOTIENT,    //差商
    FINITE_DIFFERENCE       //差分
};

class NewtonInsert{
    //外部using声明
    using Poly=MatCal::Utils::QinJiuShao;
    using PolyNode=MatCal::Utils::QinJiuShaoNode;
    using Lower=MatCal::Utils::LowerTriangularMatrix;
    using Mat=MatCal::Utils::Matrix;

//构造和析构
public:


//内部维护的属性
private:
    DifferenceType type;    //牛顿插值-差商/差分
    Poly poly;              //插值多项式
    Lower lower;            //差商/分表
};//class NewtonInsert

}//namespace MatCal::Utils



#endif//NEWTONINSERT_HPP