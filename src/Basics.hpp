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
#include"Matrix.hpp"
#include"QinJiuShao.hpp"

namespace MatCal::Algorithm::Basics{
    using QinJiuShao = MatCal::Utils::QinJiuShao;

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
    static Result_least_square least_square(int degree,std::vector<double>& x,
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
    static Result_least_square least_square(int degree,std::vector<double>& x,
    std::vector<double>& y){
        std::vector<double> weights(x.size(),1.0);
        return least_square(degree,x,y,weights);
    }

    // 版本3：带权重和选择项的多项式拟合
    static Result_least_square least_square(int degree, std::vector<double>& x,
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
    static Result_least_square least_square(int degree, std::vector<double>& x,
        std::vector<double>& y, std::vector<bool>& selects) {
        
        std::vector<double> weights(x.size(), 1.0);
        return least_square(degree, x, y, weights, selects);
    }

};


}//namespace MatCal::Algorithm::Basics

#endif//BASICS_HPP