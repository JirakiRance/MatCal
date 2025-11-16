#ifndef BASICS_HPP
#define BASICS_HPP

/*

这个文件用来定义一些基础的算法，如数值求导积分

*/

namespace MatCal{
    namespace Algorithm{
        namespace Basics{
            //数值求导(或者偏导)，y=f(x),y=f(x1,x2...)
            class Derivative;
        }
    }
}

#include<functional>
#include<vector>

namespace MatCal::Algorithm::Basics{

//c++没有静态类，这个类也就写一点静态方法
class Derivative{
public:
    using Func_y = std::function<double(double)>;
    using Func_F = std::function<double(std::vector<double>)>;

    //Func_y求导dy/dx,默认向后取eps
    static double dy_dx(Func_y _func,double x,double eps=1e-6){
        return (
            ( _func(x+eps) - _func(x) ) / eps
        );
    }
    static double dy_dx_center(Func_y _func,double x,double eps=1e-6){
        return (
            ( _func(x+0.5*eps) - _func(x-0.5*eps) ) / eps
        );
    }

    //Func_F求偏导pF/px(指定某一个位置)，本方法不对索引越界负责
    static double pF_px(Func_F _func,std::vector<double>& xs,int i,double eps=1e-6){
        std::vector<double> delta_vector=xs;
        delta_vector[i]+=eps;
        return (
            ( _func(delta_vector) - _func(xs) ) / eps
        );
    }

    //Func_F求全导dF_dx
    static double dF_dx(Func_F _func,std::vector<double>& xs,double eps=1e-6){
        double sum=0;
        for(int i=0;i<xs.size();++i)
            sum+=pF_px(_func,xs,i,eps);
        return sum;
    }


};


}//namespace MatCal::Algorithm::Basics

#endif//BASICS_HPP