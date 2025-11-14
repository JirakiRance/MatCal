#ifndef ITERATION_HPP
#define ITERATION_HPP

/*

本文件实现了一些迭代算法

*/

#include<stdexcept>
#include<cmath>
#include<functional>
#include<vector>


namespace MatCal{
    namespace Algorithm{
        namespace Iteration{
            //二分法，支持传入一个函数，一个区间，一个，求解一个零点，可以按需求传入精度，最大迭代步数，不保证有解，无解会抛出异常
            class Bisection;
            //Picard收敛法（简单迭代法），用于迭代求解形如  x=\varphi(x)的形式，不保证有解，无解会抛出异常
            class Picard;
            //牛顿迭代法
            class Newton;
            //割线法
            class Secant;
        }
    }
}

namespace MatCal::Algorithm::Iteration{

class Bisection{
public:
    //函数指针类型定义
    using Function = std::function<double(double)>;
    //二分求解，基础版,解不收敛则会抛出异常
    static double solve(Function f, double a, double b, double epsilon = 1e-6, int maxIterations = 1000){
        //检查区间有效性
        if(a>=b){
            throw std::invalid_argument("区间 [a, b] 必须满足 a < b");
        }
        //计算两个断电的函数值
        double fa = f(a);
        double fb = f(b);
        //检查根的存在性
        if(fa * fb > 0){
            throw std::invalid_argument("函数在区间端点值同号，无法保证根的存在");
        }
        //如果端点恰好是根
        if(std::abs(fa) < epsilon) return a;
        if(std::abs(fb) < epsilon) return b;
        
        //开始求解
        double x=a;
        for(int i=0;i<maxIterations;++i){
            x=a+(b-a)/2.0;
            double fx = f(x);
            //检查收敛
            if(std::abs(fx)<epsilon||(b-a)/2.0<epsilon)
                return x;
            //更新区间
            if(fa*fx<0){
                b=x;
                fb=fx;
            }else{
                a=x;
                fa=fx;
            }
        }
        throw std::runtime_error("二分法未在最大迭代次数内收敛");
    }

    //二分求解详细版返回结果。{找到的根，迭代次数，最终误差，是否收敛}
    struct Result {
        double root;           // 找到的根
        int iterations;        // 迭代次数
        double error;          // 最终误差 |f(root)|
        bool converged;        // 是否收敛
    };
    //二分求解，详细版。解不收敛会反映在Result里，不会抛出异常
    static Result solveDetailed(Function f, double a, double b, double epsilon = 1e-6, int maxIterations = 1000){
        if(a>=b)
            throw std::invalid_argument("区间 [a, b] 必须满足 a < b");
        double fa = f(a);
        double fb = f(b);
        if(fa*fb>0)
            throw std::invalid_argument("函数在区间端点值同号，无法保证根的存在");
        //检查端点是否为根
        if(std::abs(fa)<epsilon)
            return {a, 0, std::abs(fa), true};
        if(std::abs(fb)<epsilon)
            return {b, 0, std::abs(fb), true};
        //求解
        double x=a;
        int iterations=0;
        for(iterations=0;iterations<maxIterations;++iterations){
            x=a+(b-a)/2.0;
            double fx=f(x);
            // 收敛
            if(std::abs(fx)<epsilon||(b - a)/2.0<epsilon)
                return {x, iterations + 1, std::abs(fx), true};
            // 更新区间
            if(fa*fx<0){
                b=x;
                fb=fx;
            }else{
                a=x;
                fa=fx;
            }
        }
        // 达到最大迭代次数但未收敛
        return {x, maxIterations, std::abs(f(x)), false};
    }
};

class Picard{
public:
    //函数指针类型定义
    using Function = std::function<double(double)>;
    //结果子类
    struct Result {
        double root;                    //找到的根
        int iterations;                 //迭代次数
        double error;                   //最终误差 |f(root)|
        bool converged;                 //是否收敛
        std::vector<double> series;     //求解序列
    };
    //Picard/简单迭代法，基础版，不保证有解，不收敛会抛出异常.需要传入function_varphi和初值x0
    static double solve(Function function_varphi, double x_0, double epsilon = 1e-6, int maxIterations = 1000){
        double x_1=0.0;
        double last_delta=INFINITY;
        while(maxIterations--){
            x_1=function_varphi(x_0);
            double new_delta=std::abs(x_1-x_0);
            if(new_delta<epsilon)
                return x_1;
            last_delta=new_delta;
            x_0=x_1;
        }
        throw std::runtime_error("未在最大迭代次数内收敛");
    }
    //Picard/简单迭代法，详细版，不保证有解，不收敛会抛出异常.需要传入function_varphi和初值x0
    static Result solveDetailed(Function function_varphi, double x_0, double epsilon = 1e-6, int maxIterations = 1000){
        std::vector<double> series{x_0};
        double current=x_0;
        double next=0.0;
        double last_delta=INFINITY;
        double new_delta=0.0;
        int expand_count=0;
        for(int iteration = 1;iteration<=maxIterations;++iteration){
            next=function_varphi(current);
            new_delta=std::abs(next-current);
            if(new_delta<epsilon)
                return {next,iteration,new_delta,true,series};
            if(new_delta>last_delta*1.5){
                ++expand_count;
                if(expand_count>=3){
                    throw std::runtime_error("Piscard发散,请重新选择function_varphi或初值x0");
                }
            }else{
                expand_count=0;
            }
            last_delta=new_delta;
            current=next;
            series.push_back(next);
        }
        //throw std::runtime_error("未在最大迭代次数内收敛");
        return {next,maxIterations,new_delta,false,series};
    }
    //Picard迭代-使用Aitken加速
    static Result solve_Aitken(Function function_varphi, double x_0, double epsilon = 1e-6, int maxIterations = 1000){
        //直接构造新的phi函数，然后复用solveDetailed
        Function aitken_function = [function_varphi](double x){
            double varphi_1 = function_varphi(x);
            double varphi_2 = function_varphi(varphi_1);
            double low = (varphi_2 - 2 * varphi_1 +x);
            if(std::abs(low) < 1e-12)
                return varphi_1;//除零保护，如果除零，就返回varphi_1，即不加速的原形式
            return (
                (x*varphi_2 - varphi_1 * varphi_1) / low
            );
        };
        return solveDetailed(aitken_function,x_0,epsilon,maxIterations);
    }
};

class Newton{
public:
    //函数指针类型定义
    using Function = std::function<double(double)>;
    //结果子类
    struct Result {
        double root;                    //找到的根
        int iterations;                 //迭代次数
        double error;                   //最终误差 |f(root)|
        bool converged;                 //是否收敛
        std::vector<double> series;     //求解序列
    };
    //Newton-Raphson迭代法,需要传入f(x)和f`(x)，本迭代法不负责函数求导
    static Result solve(Function f, Function df_dx, double x_0,double epsilon = 1e-6, int maxIterations = 1000) {
        std::vector<double> series{x_0};
        double current = x_0;
        double next = 0.0;
        double last_delta = INFINITY;
        double new_delta = 0.0;
        int expand_count = 0;

        for(int iteration = 1; iteration <= maxIterations; ++iteration){
            double f_val = f(current);
            double df_val = df_dx(current);
            //防止除以零
            if (std::abs(df_val) < 1e-12)
                throw std::runtime_error("Newton法失败:导数接近零，可能出现奇异点。");

            next = current - f_val / df_val;
            new_delta = std::abs(next - current);

            //收敛判断
            if(new_delta < epsilon || std::abs(f_val) < epsilon){
                return {next, iteration, std::abs(f_val), true, series};
            }
            // 发散检测
            if(new_delta > last_delta){
                ++expand_count;
                if(expand_count>=5)
                    throw std::runtime_error("Newton法发散,请重新选择初值或函数。");
            }else{
                expand_count=0;
            }
            last_delta = new_delta;
            current = next;
            series.push_back(next);
        }

        //throw std::runtime_error("Newton法未在最大迭代次数内收敛。");
        return {next,maxIterations,new_delta,false,series};
    }
    //牛顿下山法
    static Result solve_downhill(Function f, Function df_dx, double x_0,double epsilon = 1e-6, int maxIterations = 1000) {
        std::vector<double> series{x_0};
        double current = x_0;
        double f_current = f(current);
        double f_next = f_current; // 初始化，避免未定义

        for(int iteration = 1; iteration <= maxIterations; ++iteration){
            double df_val = df_dx(current);
            // 防止导数为0
            if (std::abs(df_val) < 1e-12)
                throw std::runtime_error("下山Newton法失败: 导数接近零。");
            double step = -f_current / df_val; // 标准牛顿步长
            double next = current + step;
            f_next = f(next);

            // ===== 下山策略 =====
            // 若 |f(next)| 比 |f(current)| 还大 → 发散 → 缩小步长
            double lambda = 1.0;
            int shrink_count = 0;
            while (std::abs(f_next) > std::abs(f_current) && shrink_count < 20) {
                lambda *= 0.5;  // 缩小步长
                next = current + lambda * step;
                f_next = f(next);
                ++shrink_count;
            }
            //收敛判据
            if(std::abs(next - current) < epsilon || std::abs(f_next) < epsilon){
                series.push_back(next);
                return {next, iteration, std::abs(f_next), true, series};
            }
            //更新状态
            current = next;
            f_current = f_next;
            series.push_back(current);
        }
        //未收敛：安全返回最终状态
        return {current, maxIterations, std::abs(f_next), false, series};
    }
};

class Secant {
public:
    //函数指针类型定义
    using Function = std::function<double(double)>;
    //结果子类
    struct Result {
        double root;
        int iterations;
        double error;
        bool converged;
        std::vector<double> series;
    };
    //双点割线法
    static Result solve_two_point(Function f, double x0, double x1,double epsilon = 1e-6, int maxIterations = 1000) {
        std::vector<double> series{x0, x1};
        double f0 = f(x0);
        double f1 = f(x1);

        if (std::abs(f0) < epsilon)
            return {x0, 0, std::abs(f0), true, series};
        if (std::abs(f1) < epsilon)
            return {x1, 0, std::abs(f1), true, series};

        for (int iteration = 1; iteration <= maxIterations; ++iteration) {
            double denominator = f1 - f0;
            if (std::abs(denominator) < 1e-12)
                throw std::runtime_error("割线法失败：分母接近零。");

            double x2 = x1 - f1 * (x1 - x0) / denominator;
            double f2 = f(x2);

            if (std::abs(x2 - x1) < epsilon || std::abs(f2) < epsilon) {
                series.push_back(x2);
                return {x2, iteration, std::abs(f2), true, series};
            }

            x0 = x1;
            f0 = f1;
            x1 = x2;
            f1 = f2;
            series.push_back(x2);
        }
        return {x1, maxIterations, std::abs(f1), false, series};
    }

    //单点割线法（差分牛顿法）
    static Result solve_one_point(Function f, double x0,double h = 1e-4, double epsilon = 1e-6, int maxIterations = 1000) {
        std::vector<double> series{x0};
        double current = x0;
        double prev = x0 - h;
        double f_current = f(current);
        double f_prev = f(prev);

        for (int iteration = 1; iteration <= maxIterations; ++iteration) {
            double denominator = f_current - f_prev;
            if (std::abs(denominator) < 1e-12)
                throw std::runtime_error("单点割线法失败：差分分母接近零。");

            double next = current - f_current * (current - prev) / denominator;
            double f_next = f(next);

            if (std::abs(next - current) < epsilon || std::abs(f_next) < epsilon) {
                series.push_back(next);
                return {next, iteration, std::abs(f_next), true, series};
            }

            prev = current;
            f_prev = f_current;
            current = next;
            f_current = f_next;
            series.push_back(current);
        }

        return {current, maxIterations, std::abs(f_current), false, series};
    }
};


}//namespace MatCal::Algorithm::Iteration
#endif//ITERATION_HPP