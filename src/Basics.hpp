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

    //计算线性的基于多项式的最小二乘法
    //版本1：带权重的完整多项式
    static Result_least_square solve(int degree,std::vector<double>& x,
    std::vector<double>& y,
    std::vector<double>& weights ){
        if(degree <= 0 )
            throw std::invalid_argument("degree should > 0!");
        //先进行inner_cal
        MatCal::Utils::Matrix A(degree + 1,degree +1);
        MatCal::Utils::Matrix b(degree + 1,1);
        int i=0;
        int vector_size=x.size();
        if(y.size()!=vector_size||weights.size()!=vector_size)
            throw std::runtime_error("vector_size not match!");
        //创建A  b
        for(;i<=degree;++i){
            for(int j=0;j<=i;++j){
                double sum_a=0;
                double sum_b=0;
                double tmp=0;
                for(int k=0;k<vector_size;++k){
                    tmp=weights[k] * std::pow(x[k],i+j);
                    sum_a+=tmp;
                    if(j==0) 
                        sum_b+=(tmp * y[k]);
                }
                A.set(i,j,sum_a);
                A.set(j,i,sum_a);
                if(j==0)
                    b.set(i,0,sum_b);
            }
        }
        //求解Ax=b
        try{
            try{
                auto ret = MatCal::Algorithm::Matrix::Gauss_Seidel(A,b);
                std::vector<double> coee(degree+1);
                std::vector<std::pair<int,double>> for_qinjiushao(degree+1);
                for(int i=0;i<=degree;++i){
                    coee[i] = ret.root.get(i,0);
                    for_qinjiushao[i] = std::make_pair(i,coee[i]);
                }
                QinJiuShao poly(for_qinjiushao);
                return {coee,A,b,poly,"success!"};
            }catch(std::exception&e_gauss){
                try{
                    auto ret = MatCal::Algorithm::Matrix::solve_columnElimination(A,b);
                    std::vector<double> coee(degree+1);
                    std::vector<std::pair<int,double>> for_qinjiushao(degree+1);
                    for(int i=0;i<=degree;++i){
                        coee[i] = ret->get(i,0);
                        for_qinjiushao[i] = std::make_pair(i,coee[i]);
                    }
                    QinJiuShao poly(for_qinjiushao);
                    std::string msg = "tried guass iteration to solve,but failed:\n";
                    msg+=e_gauss.what();
                    msg+="\ntry columnElimination:success!";
                    return {coee,A,b,poly,msg};
                }catch(std::exception & e_elimination){
                    std::string msg="solve linear failed!\n";
                    msg+=e_gauss.what();
                    msg+="\n";
                    msg+=e_elimination.what();
                    return{std::vector<double>(),A,b,QinJiuShao({0,1}),msg};
                }
            }
        }catch(std::exception& e){
            return {std::vector<double>(),A,b,QinJiuShao({0,1}),e.what()};
        }
    }

    //版本2：无权重的完整多项式
    static Result_least_square solve(int degree,std::vector<double>& x,
    std::vector<double>& y){
        std::vector<double> weights(x.size(),1.0);
        return solve(degree,x,y,weights);
    }

    // 版本3：带权重和选择项的多项式拟合
    static Result_least_square solve(int degree, std::vector<double>& x,
        std::vector<double>& y, std::vector<double>& weights, 
        std::vector<bool>& selects) {
        
        if(degree <= 0)
            throw std::invalid_argument("degree should > 0!");
        
        if(selects.size() != degree + 1)
            throw std::runtime_error("selects size should be degree+1!");
        
        // 计算实际使用的项数
        int actual_degree = 0;
        for(bool selected : selects) {
            if(selected) actual_degree++;
        }
        if(actual_degree == 0)
            throw std::runtime_error("at least one term should be selected!");
        
        int vector_size = x.size();
        if(y.size() != vector_size || weights.size() != vector_size)
            throw std::runtime_error("vector_size not match!");
        
        // 创建索引映射：原次数 -> 新索引
        std::vector<int> index_map(degree + 1, -1);
        int current_index = 0;
        for(int i = 0; i <= degree; ++i) {
            if(selects[i]) {
                index_map[i] = current_index++;
            }
        }
        
        MatCal::Utils::Matrix A(actual_degree, actual_degree);
        MatCal::Utils::Matrix b(actual_degree, 1);
        
        // 创建A和b矩阵
        for(int i = 0; i <= degree; ++i) {
            if(!selects[i]) continue;
            int new_i = index_map[i];
            
            for(int j = 0; j <= degree; ++j) {
                if(!selects[j]) continue;
                int new_j = index_map[j];
                
                // 只计算下三角或对角线
                if(new_j > new_i) continue;
                
                double sum_a = 0;
                for(int k = 0; k < vector_size; ++k) {
                    sum_a += weights[k] * std::pow(x[k], i + j);
                }
                A.set(new_i, new_j, sum_a);
                A.set(new_j, new_i, sum_a);  // 对称矩阵
            }
            
            // 计算b向量
            double sum_b = 0;
            for(int k = 0; k < vector_size; ++k) {
                sum_b += weights[k] * y[k] * std::pow(x[k], i);
            }
            b.set(new_i, 0, sum_b);
        }
        
        // 求解Ax=b
        try{
            try{
                auto ret = MatCal::Algorithm::Matrix::Gauss_Seidel(A,b);
                std::vector<double> coee(actual_degree);
                std::vector<std::pair<int,double>> for_qinjiushao;
                
                for(int i = 0; i <= degree; ++i) {
                    if(selects[i]) {
                        int new_index = index_map[i];
                        coee[new_index] = ret.root.get(new_index, 0);
                        for_qinjiushao.push_back(std::make_pair(i, coee[new_index]));
                    }
                }
                
                QinJiuShao poly(for_qinjiushao);
                return {coee, A, b, poly, "success!"};
            }catch(std::exception& e_gauss){
                try{
                    auto ret = MatCal::Algorithm::Matrix::solve_columnElimination(A,b);
                    std::vector<double> coee(actual_degree);
                    std::vector<std::pair<int,double>> for_qinjiushao;
                    
                    for(int i = 0; i <= degree; ++i) {
                        if(selects[i]) {
                            int new_index = index_map[i];
                            coee[new_index] = ret->get(new_index, 0);
                            for_qinjiushao.push_back(std::make_pair(i, coee[new_index]));
                        }
                    }
                    
                    QinJiuShao poly(for_qinjiushao);
                    std::string msg = "tried guass iteration to solve,but failed:\n";
                    msg += e_gauss.what();
                    msg += "\ntry columnElimination:success!";
                    return {coee, A, b, poly, msg};
                }catch(std::exception & e_elimination){
                    std::string msg = "solve linear failed!\n";
                    msg += e_gauss.what();
                    msg += "\n";
                    msg += e_elimination.what();
                    return {std::vector<double>(), A, b, QinJiuShao({0,1}), msg};
                }
            }
        }catch(std::exception& e){
            return {std::vector<double>(), A, b, QinJiuShao({0,1}), e.what()};
        }
    }

    // 版本4：无权重的选择项多项式拟合
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
