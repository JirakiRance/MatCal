#ifndef ITERATION_HPP
#define ITERATION_HPP

/*

本文件实现了一些迭代算法

*/

#include<stdexcept>
#include<cmath>
#include<functional>
#include<vector>
#include<string>


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
            //Newton求解非线性方程组，需要借助Matrix的方法
            class NewtonForEquations;
        }
    }
}


#include"Matrix.hpp"
#include"Basics.hpp"
#include"MatCal/Roots/Roots.hpp"
#include"MatCal/Nonlinear/Nonlinear.hpp"

namespace MatCal::Algorithm::Iteration{

namespace detail {

inline MatCal::Roots::RootOptions legacy_root_options(double epsilon, int maxIterations) {
    MatCal::Roots::RootOptions options;
    options.absolute_tolerance = epsilon;
    options.relative_tolerance = 0.0;
    options.max_iterations = maxIterations;
    return options;
}

inline void throw_for_legacy_root_failure(const MatCal::Roots::RootResult& result) {
    using MatCal::Roots::RootStatus;
    if (result.success() || result.diagnostic.status == RootStatus::not_converged) {
        return;
    }
    if (result.diagnostic.status == RootStatus::invalid_input ||
        result.diagnostic.status == RootStatus::bracket_error) {
        throw std::invalid_argument(result.diagnostic.message);
    }
    throw std::runtime_error(result.diagnostic.message);
}

} // namespace detail

class Bisection{
public:
    //函数指针类型定义
    using Function = std::function<double(double)>;
    //二分求解，基础版,解不收敛则会抛出异常
    static double solve(Function f, double a, double b, double epsilon = 1e-6, int maxIterations = 1000){
        auto result = solveDetailed(f, a, b, epsilon, maxIterations);
        if(result.converged)
            return result.root;
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
        auto result = MatCal::Roots::solve_bisection(f, a, b, detail::legacy_root_options(epsilon, maxIterations));
        detail::throw_for_legacy_root_failure(result);
        return {result.value,
                result.metrics.iterations,
                result.metrics.residual,
                result.converged};
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
        auto result = solveDetailed(function_varphi, x_0, epsilon, maxIterations);
        if(result.converged)
            return result.root;
        throw std::runtime_error("未在最大迭代次数内收敛");
    }
    //Picard/简单迭代法，详细版，不保证有解，不收敛会抛出异常.需要传入function_varphi和初值x0
    static Result solveDetailed(Function function_varphi, double x_0, double epsilon = 1e-6, int maxIterations = 1000){
        auto result = MatCal::Roots::solve_picard(function_varphi, x_0, detail::legacy_root_options(epsilon, maxIterations));
        detail::throw_for_legacy_root_failure(result);
        return {result.value,
                result.metrics.iterations,
                result.metrics.final_step,
                result.converged,
                result.series};
    }
    //Picard迭代-使用Aitken加速
    static Result solve_Aitken(Function function_varphi, double x_0, double epsilon = 1e-6, int maxIterations = 1000){
        auto result = MatCal::Roots::solve_picard_aitken(function_varphi, x_0, detail::legacy_root_options(epsilon, maxIterations));
        detail::throw_for_legacy_root_failure(result);
        return {result.value,
                result.metrics.iterations,
                result.metrics.final_step,
                result.converged,
                result.series};
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
        auto result = MatCal::Roots::solve_newton(f, df_dx, x_0, detail::legacy_root_options(epsilon, maxIterations));
        detail::throw_for_legacy_root_failure(result);
        return {result.value,
                result.metrics.iterations,
                result.metrics.residual,
                result.converged,
                result.series};
    }
    //牛顿下山法
    static Result solve_downhill(Function f, Function df_dx, double x_0,double epsilon = 1e-6, int maxIterations = 1000) {
        auto result = MatCal::Roots::solve_downhill_newton(f, df_dx, x_0, detail::legacy_root_options(epsilon, maxIterations));
        detail::throw_for_legacy_root_failure(result);
        return {result.value,
                result.metrics.iterations,
                result.metrics.residual,
                result.converged,
                result.series};
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
        auto result = MatCal::Roots::solve_secant_two_point(f, x0, x1, detail::legacy_root_options(epsilon, maxIterations));
        detail::throw_for_legacy_root_failure(result);
        return {result.value,
                result.metrics.iterations,
                result.metrics.residual,
                result.converged,
                result.series};
    }

    //单点割线法（差分牛顿法）
    static Result solve_one_point(Function f, double x0,double h = 1e-4, double epsilon = 1e-6, int maxIterations = 1000) {
        auto result = MatCal::Roots::solve_secant_one_point(f, x0, h, detail::legacy_root_options(epsilon, maxIterations));
        detail::throw_for_legacy_root_failure(result);
        return {result.value,
                result.metrics.iterations,
                result.metrics.residual,
                result.converged,
                result.series};
    }
};

class NewtonForEquations{
public:
    using Function = std::function<double(const std::vector<double>&)>;

    struct Result {
        std::vector<double> root;
        int iterations;
        double error;
        bool converged;
        std::string message;
    };

    static Result solve(int n,std::vector<Function>& Funcs,std::vector<double> xs,int maxIterations = 20,double epsilon = 1e-6){
        if(n<=0)
            throw std::invalid_argument("n must > 0!");
        const auto expected_size = static_cast<std::size_t>(n);
        if(Funcs.size()!=expected_size)
            throw std::invalid_argument("Number of functions must = n");
        if(xs.size()!=expected_size)
            throw std::invalid_argument("Initial xs size must = n");
        if(maxIterations <= 0)
            throw std::invalid_argument("maxIterations must > 0");
        if(!std::isfinite(epsilon) || epsilon <= 0.0)
            throw std::invalid_argument("epsilon must be finite and > 0");

        MatCal::Nonlinear::ResidualFunction residual = [Funcs, expected_size](const std::vector<double>& values) {
            if(values.size() != expected_size)
                throw std::invalid_argument("NewtonForEquations residual input size mismatch");
            std::vector<double> ret(expected_size);
            for(std::size_t i = 0; i < expected_size; ++i)
                ret[i] = Funcs[i](values);
            return ret;
        };

        MatCal::Nonlinear::NonlinearOptions options;
        options.absolute_tolerance = epsilon;
        options.relative_tolerance = 0.0;
        options.finite_difference_step = epsilon;
        options.max_iterations = static_cast<std::size_t>(maxIterations);
        auto result = MatCal::Nonlinear::solve_newton_system_finite_difference(residual, xs, options);

        if(result.success()){
            return {result.solution,
                    static_cast<int>(result.metrics.iterations),
                    result.metrics.residual_norm,
                    true,
                    result.diagnostic.message.empty() ? "Converged successfully" : result.diagnostic.message};
        }
        return {std::vector<double>(),
                static_cast<int>(result.metrics.iterations),
                result.metrics.residual_norm,
                false,
                result.diagnostic.message};
    }
};
}//namespace MatCal::Algorithm::Iteration
#endif//ITERATION_HPP
