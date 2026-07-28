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
#include<cmath>
#include<algorithm>
#include<climits>
#include"Matrix.hpp"
#include"QinJiuShao.hpp"
#include"Insert.hpp"
#include"MatCal/Calculus/Calculus.hpp"
#include"MatCal/ODE/ODE.hpp"
#include"MatCal/LeastSquares/LeastSquares.hpp"

namespace MatCal::Algorithm::Basics{
    using QinJiuShao = MatCal::Utils::QinJiuShao;

//圆周率pi
const static double PI =  3.14159265358979323846;
//自然底数e
const static double E = 2.71828182845904523536;

//求解线性系统Ax=b。使用方法：Gauss-Sedel,ColmunElimination,LU
static std::pair<MatCal::Utils::Matrix,std::string> solve_Linear_System(
    MatCal::Utils::Matrix& A,
    MatCal::Utils::Matrix& b){
    std::string msg="";

    try{    //Gauss_Seidel
        msg+="method:Gauss_Seidel:";
        auto ret1 = MatCal::Algorithm::Matrix::Gauss_Seidel(A,b);
        if(ret1.converged){
            msg+="success!";
            return std::make_pair(ret1.root,msg);
        }else{
            msg+="not converged!\n";
        }
    }catch(std::exception& e){
        msg+="failed---->";
        msg+=e.what();
        msg+="\n";
    }
    
    try{    //ColmunElimination
        msg+="method:ColmunElimination:";
        auto ret2 = MatCal::Algorithm::Matrix::solve_columnElimination(A,b);
        MatCal::Utils::Matrix* copy_child_2 = dynamic_cast<MatCal::Utils::Matrix*>(ret2.get());
        msg+="success!";
        return std::make_pair(*copy_child_2,msg);
    }catch(std::exception& e){
        msg+="failed!---->";
        msg+=e.what();
        msg+="\n";
    }
    
    try{    //LU Decompose
        msg+="method:LU Decompose:";
        auto ret3 = MatCal::Algorithm::Matrix::LU_Decompose(A).solve(b);
        MatCal::Utils::Matrix* copy_child_3 = dynamic_cast<MatCal::Utils::Matrix*>(ret3.get());
        msg+="success!";
        return std::make_pair(*copy_child_3,msg);
    }catch(std::exception& e){
        msg+="failed---->";
        msg+=e.what();
        msg+="\n";
    }
    msg+="all methods failed!";
    return std::make_pair(MatCal::Utils::Matrix(),msg);
}


//指定函数求导
class Derivative{
public:
    using Func_y = std::function<double(double)>;
    using Func_F = std::function<double(std::vector<double>)>;

    //Func_y求导dy/dx,默认向后取eps
    static double dy_dx(Func_y _func,double x,double eps=1e-6){
        auto result = MatCal::Calculus::forward_difference(_func, x, eps);
        if(!result.success())
            throw std::invalid_argument(result.diagnostic.message);
        return result.value;
    }
    static double dy_dx_center(Func_y _func,double x,double eps=1e-6){
        auto result = MatCal::Calculus::central_difference(_func, x, eps);
        if(!result.success())
            throw std::invalid_argument(result.diagnostic.message);
        return result.value;
    }

    //Func_F求偏导pF/px(指定某一个位置)，本方法不对索引越界负责
    static double pF_px(Func_F _func,std::vector<double>& xs,int i,double eps=1e-6){
        if(i < 0)
            throw std::out_of_range("partial derivative coordinate is out of range");
        MatCal::Calculus::MultivariateFunction adapter = [_func](const std::vector<double>& values) {
            return _func(values);
        };
        auto result = MatCal::Calculus::partial_difference(adapter, xs, static_cast<std::size_t>(i), eps);
        if(!result.success())
            throw std::invalid_argument(result.diagnostic.message);
        return result.value;
    }

    //Func_F求全导dF_dx
    static double dF_dx(Func_F _func,std::vector<double>& xs,double eps=1e-6){
        MatCal::Calculus::MultivariateFunction adapter = [_func](const std::vector<double>& values) {
            return _func(values);
        };
        auto result = MatCal::Calculus::gradient(adapter, xs, eps);
        if(!result.success())
            throw std::invalid_argument(result.diagnostic.message);
        double sum = 0.0;
        for(double value : result.values)
            sum += value;
        return sum;
    }
};

//最小二乘法
class Least_Square{
public:
    class Result_least_square{
        public:
        std::vector<double> coee;
        MatCal::Utils::Matrix A;
        MatCal::Utils::Matrix b;
        MatCal::Utils::QinJiuShao poly;
        std::string msg;
    };

private:
    static MatCal::Utils::Matrix to_legacy_matrix(const MatCal::Linalg::DenseMatrix& matrix){
        MatCal::Utils::Matrix result(static_cast<int>(matrix.rows()), static_cast<int>(matrix.cols()));
        for(std::size_t r = 0; r < matrix.rows(); ++r)
            for(std::size_t c = 0; c < matrix.cols(); ++c)
                result.set(static_cast<int>(r), static_cast<int>(c), matrix(r, c));
        return result;
    }

    static MatCal::Utils::Matrix to_legacy_matrix(const MatCal::Linalg::Vector& vector){
        MatCal::Utils::Matrix result(static_cast<int>(vector.size()), 1);
        for(std::size_t i = 0; i < vector.size(); ++i)
            result.set(static_cast<int>(i), 0, vector[i]);
        return result;
    }

    static MatCal::Utils::QinJiuShao to_legacy_polynomial(const MatCal::Polynomial::Polynomial& polynomial){
        std::vector<std::pair<int,double>> terms;
        for(const auto& term : polynomial.terms_descending()){
            if(term.first > static_cast<std::size_t>(std::numeric_limits<int>::max()))
                throw std::length_error("legacy polynomial degree exceeds int range");
            terms.emplace_back(static_cast<int>(term.first), term.second);
        }
        if(terms.empty())
            terms.emplace_back(0, 0.0);
        return MatCal::Utils::QinJiuShao(terms);
    }

    static Result_least_square map_result(const MatCal::LeastSquares::LeastSquaresResult& result){
        MatCal::Utils::Matrix A = to_legacy_matrix(result.normal_matrix);
        MatCal::Utils::Matrix b = to_legacy_matrix(result.rhs);
        MatCal::Utils::QinJiuShao poly = to_legacy_polynomial(result.polynomial);
        if(result.success())
            return {result.coefficients, A, b, poly, "success!"};
        return {std::vector<double>(), A, b, poly, result.diagnostic.message};
    }

public:
    static Result_least_square solve(int degree,std::vector<double>& x,
    std::vector<double>& y,
    std::vector<double>& weights ){
        if(degree <= 0 )
            throw std::invalid_argument("degree should > 0!");
        if(y.size()!=x.size()||weights.size()!=x.size())
            throw std::runtime_error("vector_size not match!");
        auto result = MatCal::LeastSquares::fit_polynomial_degree(degree, x, y, weights);
        return map_result(result);
    }

    static Result_least_square solve(int degree,std::vector<double>& x,
    std::vector<double>& y){
        std::vector<double> weights(x.size(),1.0);
        return solve(degree,x,y,weights);
    }

    static Result_least_square solve(int degree, std::vector<double>& x,
        std::vector<double>& y, std::vector<double>& weights,
        std::vector<bool>& selects) {
        if(degree <= 0)
            throw std::invalid_argument("degree should > 0!");
        if(selects.size() != static_cast<std::size_t>(degree) + 1)
            throw std::runtime_error("selects size should be degree+1!");
        if(y.size() != x.size() || weights.size() != x.size())
            throw std::runtime_error("vector_size not match!");
        bool any_selected = false;
        for(bool selected : selects)
            any_selected = any_selected || selected;
        if(!any_selected)
            throw std::runtime_error("at least one term should be selected!");
        auto result = MatCal::LeastSquares::fit_polynomial_selected(degree, x, y, weights, selects);
        return map_result(result);
    }

    static Result_least_square solve(int degree, std::vector<double>& x,
        std::vector<double>& y, std::vector<bool>& selects) {
        std::vector<double> weights(x.size(), 1.0);
        return solve(degree, x, y, weights, selects);
    }

};

//正交多项式
class OrthogonalPolynomials{
public:
    //Chebyshev
    static QinJiuShao Chebyshev(int n, bool second = false){
        if(n<0)
            throw std::invalid_argument("Orthogonal poly n should at least >= 0");
        int u1 = 1;
        if(second) u1=2;
        //dp数组
        QinJiuShao dp[2] = {
            QinJiuShao({{0,1}}),
            QinJiuShao({{1,u1}}),
        };
        if(n<2) return dp[n];
        QinJiuShao _2x({{1,2}});
        QinJiuShao ret = _2x;
        int target = 2;
        while(target<=n){
            ret = _2x * dp[1]- dp[0];
            dp[0] = dp[1];
            dp[1] = ret;
            ++target;
        }
        return ret;
    }

    //第一类Chebyshev多项式的零点  [-1,1]  之间
    static std::vector<double> ChebyshevZeros(int n){
        if (n <= 0)
            throw std::invalid_argument("to calc ChebyshevZeros, n should > 0 (now is " + std::to_string(n) + ")");
        std::vector<double> ret(n,0);
        for(int i=1;i<=n;++i){
            ret[i-1]=std::cos(
                PI * (2 * i -1) / (2 * n)
            );
        }
        return ret;
    }

    //Legendre
    static QinJiuShao Legendre(int n){
        if(n<0)
            throw std::invalid_argument("Orthogonal poly n should at least >= 0");
        //dp数组
        QinJiuShao dp[2] = {
            QinJiuShao({{0,1}}),
            QinJiuShao({{1,1}}),
        };
        if(n<2) return dp[n];
        QinJiuShao _x({{1,1}});
        QinJiuShao ret = _x;
        int target = 2;
        while(target<=n){
            ret = _x * (static_cast<double>(2*target-1)/(target)) *dp[1]- (static_cast<double>(target-1)/(target)) * dp[0];
            dp[0] = dp[1];
            dp[1] = ret;
            ++target;
        }
        return ret;
    }
};//class OrthogonalPolynomials


//数值积分
class NumericalIntegration{
public:
    //指定函数，按定义积分
    static double Instant(std::function<double(double)> _func,double a,double b,double eps=1e-6){
        auto result = MatCal::Calculus::integrate_instant(_func, a, b, eps);
        if(!result.success())
            throw std::invalid_argument(result.diagnostic.message);
        return result.value;
    }

    //NewtonCotes
    static double NewtonCotes(std::function<double(double)> func, double a, double b, int n = 4) {
        auto result = MatCal::Calculus::integrate_newton_cotes(func, a, b, n);
        if(!result.success())
            throw std::invalid_argument(result.diagnostic.message);
        return result.value;
    }

    //CompositeNewtonCotes
    static double CompositeNewtonCotes(std::function<double(double)> func, double a, double b, int segments = 10, int n = 4) {
        auto result = MatCal::Calculus::integrate_composite_newton_cotes(func, a, b, segments, n);
        if(!result.success())
            throw std::invalid_argument(result.diagnostic.message);
        return result.value;
    }

    //Romberg
    static std::pair<double,MatCal::Utils::Matrix> Romberg(std::function<double(double)> func,double a,double b,double eps=1e-6,int maxIterations =20){
        MatCal::Calculus::IntegrationOptions options;
        options.tolerance = eps;
        options.max_iterations = maxIterations;
        auto result = MatCal::Calculus::integrate_romberg(func, a, b, options);
        int rows = static_cast<int>(std::max<std::size_t>(result.table.size(), 1));
        MatCal::Utils::Matrix sheet(rows, 4);
        for(int r = 0; r < rows; ++r){
            for(int c = 0; c < 4 && r < static_cast<int>(result.table.size()) && c < static_cast<int>(result.table[r].size()); ++c)
                sheet.set(r, c, result.table[r][c]);
        }
        if(!result.success())
            throw std::runtime_error(result.diagnostic.message);
        return std::make_pair(result.value, sheet);
    }//Romberg

    //离散数据表
    static double NewtonCotes(std::vector<std::pair<double, double>>& data, double a, double b, int n = 4) {
        MatCal::Algorithm::Insert::LagrangeInsert lagrange(data);
        auto interp_func = lagrange.getPoly().toFunction();
        return NewtonCotes(interp_func, a, b, n);
    }
    static double CompositeNewtonCotes(std::vector<std::pair<double, double>>& data, double a, double b, int segments = 10, int n = 4) {
        MatCal::Algorithm::Insert::LagrangeInsert lagrange(data);
        auto interp_func = lagrange.getPoly().toFunction();
        return CompositeNewtonCotes(interp_func, a, b, segments, n);
    }
    static std::pair<double,MatCal::Utils::Matrix> Romberg(std::vector<std::pair<double, double>>& data,double a,double b,double eps=1e-6,int maxIterations =20){
        MatCal::Algorithm::Insert::LagrangeInsert lagrange(data);
        auto interp_func = lagrange.getPoly().toFunction();
        return Romberg(interp_func, a, b, eps , maxIterations);
    }

public:
    const static std::vector<std::vector<double>> CotesSheet;
};
inline const std::vector<std::vector<double>> NumericalIntegration::CotesSheet = {
    {1.0/2,     1.0/2},
    {1.0/6,     4.0/6,      1.0/6},
    {1.0/8,     3.0/8,      3.0/8,      1.0/8},
    {7.0/90,    16.0/45,    2.0/15,     16.0/45,      7.0/90},
    {19.0/288,  25.0/96,    25.0/144,   25.0/144,     25.0/96,    19.0/288},
    {41.0/840,  9.0/35,     9.0/280,    34.0/105,     9.0/280,    9.0/35,     41.0/840},
    {751.0/17280,3577.0/17280,1323.0/17280,2989.0/17280,2989.0/17280,1323.0/17280,3577.0/17280,751.0/17280}
};

//常微分方程求解(ODE)
class ODE{
private:
    static MatCal::Utils::Matrix matrix_from_trajectory(const std::vector<std::vector<double>>& trajectory){
        if(trajectory.empty())
            return MatCal::Utils::Matrix();
        MatCal::Utils::Matrix result(static_cast<int>(trajectory.size()), static_cast<int>(trajectory[0].size()));
        for(int r = 0; r < static_cast<int>(trajectory.size()); ++r)
            for(int c = 0; c < static_cast<int>(trajectory[r].size()); ++c)
                result[r][c] = trajectory[r][c];
        return result;
    }

    static MatCal::ODE::Rhs adapt_legacy_rhs(int n, std::vector<std::function<double(std::vector<double>&)>>& funcs){
        return [n, funcs](double t, const std::vector<double>& state) {
            std::vector<double> row(static_cast<std::size_t>(n) + 1);
            row[0] = t;
            for(int i = 0; i < n; ++i)
                row[static_cast<std::size_t>(i) + 1] = state[static_cast<std::size_t>(i)];
            std::vector<double> values(static_cast<std::size_t>(n));
            for(int i = 0; i < n; ++i){
                auto mutable_row = row;
                values[static_cast<std::size_t>(i)] = funcs[static_cast<std::size_t>(i)](mutable_row);
            }
            return values;
        };
    }

public:
    //SimpleEuler
    static MatCal::Utils::Matrix SimpleEuler(int n,std::vector<std::function<double(std::vector<double>&)>>& funcs,std::vector<double>& inits,double h=1e-2,int count=100){
        if(n < 1)
            throw std::invalid_argument("n should be >= 1 !");
        const auto expected_state_size = static_cast<std::size_t>(n);
        const auto expected_row_size = expected_state_size + 1;
        if(funcs.size()!=expected_state_size)
            throw std::invalid_argument("func.size not match n!   func.size:" + std::to_string(funcs.size())+ "  n:" + std::to_string(n) );
        if(inits.size()!=expected_row_size)
            throw std::invalid_argument("inits.size not match n+1!   inits.size:" + std::to_string(inits.size())+ "  n+1:" + std::to_string(n+1) );
        if(h <= 0)
            throw std::invalid_argument("h should be > 0 !");
        if(count < 1)
            throw std::invalid_argument("count should be >= 1 !");
        std::vector<double> state(inits.begin() + 1, inits.end());
        auto result = MatCal::ODE::integrate_euler(adapt_legacy_rhs(n, funcs), inits[0], state, h, count);
        if(!result.success())
            throw std::runtime_error(result.diagnostic.message);
        return matrix_from_trajectory(result.trajectory);
    }//SimpleEuler

    //imporovedEuler,前一个为算出来的值，后面一个是中间的临时值，供手写用(没什么卵用,仅供手写做题用)
    static std::pair<MatCal::Utils::Matrix,MatCal::Utils::Matrix> Euler(int n,std::vector<std::function<double(std::vector<double>&)>>& funcs,std::vector<double>& inits,double h=1e-2,int count=100){
        if(n < 1)
            throw std::invalid_argument("n should be >= 1 !");
        const auto expected_state_size = static_cast<std::size_t>(n);
        const auto expected_row_size = expected_state_size + 1;
        if(funcs.size()!=expected_state_size)
            throw std::invalid_argument("func.size not match n!   func.size:" + std::to_string(funcs.size())+ "  n:" + std::to_string(n) );
        if(inits.size()!=expected_row_size)
            throw std::invalid_argument("inits.size not match n+1!   inits.size:" + std::to_string(inits.size())+ "  n+1:" + std::to_string(n+1) );
        if(h <= 0)
            throw std::invalid_argument("h should be > 0 !");
        if(count < 1)
            throw std::invalid_argument("count should be >= 1 !");
        MatCal::Utils::Matrix ret_true(count+1,n+1);
        MatCal::Utils::Matrix ret_temp(count+1,n+1);
        for(int i=0;i<=n;++i){
            ret_temp[0][i] = inits[i];
            ret_true[0][i] = inits[i];
        }
        auto rhs = adapt_legacy_rhs(n, funcs);
        double t = inits[0];
        std::vector<double> state(inits.begin() + 1, inits.end());
        for(int cnt = 1; cnt <= count; ++cnt){
            auto step = MatCal::ODE::improved_euler_step(rhs, t, state, h);
            if(!step.first.success())
                throw std::runtime_error(step.first.diagnostic.message);
            t = step.first.t;
            state = step.first.state;
            ret_true[cnt][0] = t;
            ret_temp[cnt][0] = t;
            for(int k = 1; k <= n; ++k){
                ret_true[cnt][k] = state[static_cast<std::size_t>(k - 1)];
                ret_temp[cnt][k] = step.second[static_cast<std::size_t>(k - 1)];
            }
        }
        return std::make_pair(ret_true,ret_temp);
    }//Euler

    //RungeKutta_44,返回结果的每一代的矩阵表示，一行为一代
    static MatCal::Utils::Matrix RungeKutta_44(int n,std::vector<std::function<double(std::vector<double>&)>>& funcs,std::vector<double>& inits,double h=1e-2,int count=100){
        if(n < 1)
            throw std::invalid_argument("n should be >= 1 !");
        const auto expected_state_size = static_cast<std::size_t>(n);
        const auto expected_row_size = expected_state_size + 1;
        if(funcs.size()!=expected_state_size)
            throw std::invalid_argument("func.size not match n!   func.size:" + std::to_string(funcs.size())+ "  n:" + std::to_string(n) );
        if(inits.size()!=expected_row_size)
            throw std::invalid_argument("inits.size not match n+1!   inits.size:" + std::to_string(inits.size())+ "  n+1:" + std::to_string(n+1) );
        if(h <= 0)
            throw std::invalid_argument("h should be > 0 !");
        if(count < 1)
            throw std::invalid_argument("count should be >= 1 !");

        std::vector<double> state(inits.begin() + 1, inits.end());
        auto result = MatCal::ODE::integrate_rk4(adapt_legacy_rhs(n, funcs), inits[0], state, h, count);
        if(!result.success())
            throw std::runtime_error(result.diagnostic.message);
        return matrix_from_trajectory(result.trajectory);
    }//RungeKutta_44

};//class ODE

// ======================================================
// 数值积分（RK4）//PT项目专用方法以及命名空间，勿动!
// ======================================================
namespace Integrate {

class RK4 {
public:
    // 通用向量版本：y(n) -> dy/dt(n)
    using RHS = std::function<void(const std::vector<double>& y,
                                   std::vector<double>& dydt)>;

    // 单步 RK4：给定 y, dt 和 RHS，计算 y_next
    static void step(const RHS& f,
                     const std::vector<double>& y,
                     double dt,
                     std::vector<double>& y_out)
    {
        if(!f)
            throw std::invalid_argument("RK4::step: RHS cannot be empty");
        MatCal::ODE::Rhs rhs = [&f](double, const std::vector<double>& state) {
            std::vector<double> dydt(state.size());
            f(state, dydt);
            return dydt;
        };
        auto result = MatCal::ODE::rk4_step(rhs, 0.0, y, dt);
        if(!result.success())
            throw std::runtime_error(result.diagnostic.message);
        y_out = result.state;
    }

    // 2 维系统专用版本： (theta, omega)
    // f(theta, omega, dtheta, domega) 填写导数
    using RHS2 = std::function<void(double theta, double omega,
                                    double& dtheta, double& domega)>;

    static void step2(const RHS2& f,
                      double theta,
                      double omega,
                      double dt,
                      double& theta_out,
                      double& omega_out)
    {
        if(!f)
            throw std::invalid_argument("RK4::step2: RHS cannot be empty");
        MatCal::ODE::Rhs rhs = [&f](double, const std::vector<double>& state) {
            double dtheta = 0.0;
            double domega = 0.0;
            f(state[0], state[1], dtheta, domega);
            return std::vector<double>{dtheta, domega};
        };
        auto result = MatCal::ODE::rk4_step(rhs, 0.0, std::vector<double>{theta, omega}, dt);
        if(!result.success())
            throw std::runtime_error(result.diagnostic.message);
        theta_out = result.state[0];
        omega_out = result.state[1];
    }
};

}//namespace MatCal::Algorithm::Basics::Integrate

}//namespace MatCal::Algorithm::Basics

#endif//BASICS_HPP
