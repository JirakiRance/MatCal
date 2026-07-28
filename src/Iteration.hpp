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

/*

1. 原始问题：F(x) = 0
2. 在当前点 xₖ 处线性化：F(x) ≈ F(xₖ) + J(xₖ)·(x - xₖ)
3. 令线性化后的方程 = 0：F(xₖ) + J(xₖ)·(x - xₖ) = 0
4. 整理得到：J(xₖ)·(x - xₖ) = -F(xₖ)
5. 令 Δx = x - xₖ，则：J(xₖ)·Δx = -F(xₖ)
6. 更新：xₖ₊₁ = xₖ + Δx

*/
public:
    //函数声明
    using Function = std::function<double(const std::vector<double>&)>;
public:
    //结果子类
    struct Result {
        std::vector<double> root;       //找到的根
        int iterations;                 //迭代次数
        double error;                   //最终误差 |f(root)|
        bool converged;                 //是否收敛
        std::string message;            //消息
    };

    //本静态方法用于求解非线性方程组的根，需要传入std::vector<Function>,需要用户传入初值（内部设置初值不保证通用）,每一个Function必须是std::vector<double>参数，size必须为n，由于算法复杂度高，默认只迭代20次
    static Result solve(int n,std::vector<Function>& Funcs,std::vector<double> xs,int maxIterations = 20,double epsilon = 1e-6){
        //需要先校验n合法性
        if(n<=0)
            throw std::invalid_argument("n must > 0!");
        if(Funcs.size()!=n)
            throw std::invalid_argument("Number of functions must = n");
        if(xs.size()!=n)
            throw std::invalid_argument("Initial xs size must = n");
        //验证每个函数确实接受 n 个变量
        std::vector<double> test_input(n, 1.0);
        for (int i = 0; i < n; ++i) {
            try {
                Funcs[i](test_input);
            }catch(std::exception & e){
                std::cout<<e.what();
                throw std::invalid_argument("Function " + std::to_string(i) + " does not properly handle n variables");
            }
        }

        //待返回参数
        bool converged =false;
        double error=INFINITY;
        std::string message = "Maximum iterations reached";

        for(int iteration = 1;iteration<=maxIterations;++iteration){
            //计算函数值 F(x)
            std::vector<double> F_val(n);
            for(int i=0;i<n;++i)
                F_val[i] = Funcs[i](xs);
            //计算Jacobian矩阵
            MatCal::Utils::Matrix J(n,n);
            for(int i=0;i<n;++i)
                for (int j = 0; j < n; ++j)
                    J.set(i, j, MatCal::Algorithm::Basics::Derivative::pF_px(Funcs[i], xs, j, epsilon));
            //创建右边向量 -F(x)
            MatCal::Utils::Matrix F_mat(n, 1);
            for(int i = 0; i < n; ++i)
                F_mat.set(i, 0, -F_val[i]);

            //使用列主元消去法求解 J · Δx = -F(x)，如果失败则使用Gauss-Seidel迭代法作为备选方案
            std::unique_ptr<MatCal::Utils::AbstractMatrix> delta_x_ptr;
            bool use_backup_solver = false;
            try{
                delta_x_ptr = MatCal::Algorithm::Matrix::solve_columnElimination(J,F_mat);
            }catch(const std::exception& e){
                //列主元消去法失败，尝试使用Gauss-Seidel迭代法作为备选方案
                try{
                    auto gauss_result = MatCal::Algorithm::Matrix::Gauss_Seidel(J, F_mat, epsilon*0.1, 50);
                    if(gauss_result.converged){
                        //创建临时矩阵来存储Gauss-Seidel结果
                        auto temp_matrix = std::make_unique<MatCal::Utils::Matrix>(n, 1);
                        for(int i=0;i<n;++i)
                            temp_matrix->set(i, 0, gauss_result.root.get(i, 0));
                        delta_x_ptr = std::move(temp_matrix);
                        use_backup_solver = true;
                        message = "Used Gauss-Seidel backup solver at iteration " + std::to_string(iteration);
                    }else{
                        message = "Both direct and iterative linear solvers failed: " + std::string(e.what());
                        return {xs,iteration,error,false,message};
                    }
                }catch(const std::exception& e2){
                    message = "Linear system solve failed: " + std::string(e.what()) + ", backup also failed: " + std::string(e2.what());
                    return {xs,iteration,error,false,message};
                }
            }
            //提取 Δx 向量
            std::vector<double> delta_x(n);
            auto* delta_x_matrix = dynamic_cast<MatCal::Utils::Matrix*>(delta_x_ptr.get());
            if(!delta_x_matrix || delta_x_matrix->getCols() != 1){
                message = "Unexpected solution format from linear solver";
                return {xs,iteration,error,false,message};
            }
            for(int i = 0; i < n; ++i)
                delta_x[i] = delta_x_matrix->get(i, 0);
            //更新变量: xₖ₊₁ = xₖ + Δx
            std::vector<double> xs_new = xs;
            double max_delta = 0.0;
            for (int i = 0; i < n; ++i) {
                xs_new[i] += delta_x[i];
                max_delta = std::max(max_delta, std::abs(delta_x[i]));
            }
            //计算新误差 ||F(x_new)||∞
            error = 0.0;
            for(int i=0; i<n;++i)
                error = std::max(error, std::abs(Funcs[i](xs_new)));

            //收敛判断
            if(error < epsilon){
                converged = true;
                //如果使用了备选求解器，在成功消息中注明
                if(use_backup_solver){
                    message = "Converged successfully (with Gauss-Seidel backup)";
                }else{
                    message = "Converged successfully";
                }
                return {xs_new,iteration,error,converged,message};
            }
            //防止数值发散
            if(error > 1e10 || !std::isfinite(error)){
                message = "Numerical divergence detected";
                return {xs,iteration,error,false,message};
            }
            //更新解
            xs = xs_new;
        }
        return {xs,maxIterations,error,converged,message};
    }
};

}//namespace MatCal::Algorithm::Iteration
#endif//ITERATION_HPP
