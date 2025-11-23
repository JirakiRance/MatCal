#include <iostream>
#include <limits>
#include <vector>
#include <functional>
#include<windows.h>

#include "Basics.hpp"
#include "QinJiuShao.hpp"
#include "Insert.hpp"
#include "Iteration.hpp"
#include "Matrix.hpp"

using MatCal::Utils::QinJiuShao;
using MatCal::Utils::Matrix;
using MatCal::Utils::AbstractMatrix;

namespace AlgoBasics  = MatCal::Algorithm::Basics;
namespace AlgoInsert  = MatCal::Algorithm::Insert;
namespace AlgoIter    = MatCal::Algorithm::Iteration;
namespace AlgoMatrix  = MatCal::Algorithm::Matrix;

// ========= 通用输入 =========

void clearInput(){
    std::cin.clear();
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}

double askDouble(const std::string &msg){
    double x;
    while(true){
        std::cout << msg;
        if(std::cin >> x) return x;
        clearInput();
        std::cout << "输入错误，请输入数字。\n";
    }
}

int askInt(const std::string &msg){
    int x;
    while(true){
        std::cout << msg;
        if(std::cin >> x) return x;
        clearInput();
        std::cout << "输入错误，请输入整数。\n";
    }
}

Matrix inputMatrix(const std::string& name){
    int r = askInt(name + " 行数 = ");
    int c = askInt(name + " 列数 = ");
    Matrix m(r,c);
    for(int i=0;i<r;i++){
        for(int j=0;j<c;j++){
            double v=askDouble(name+"("+std::to_string(i)+","+std::to_string(j)+")=");
            m.set(i,j,v);
        }
    }
    return m;
}

// ========= 多项式 =========

void polynomialMenu(QinJiuShao &poly, bool &hasPoly){
    while(true){
        std::cout << "\n--- 多项式菜单 ---\n";
        std::cout << "1. 构造多项式\n";
        std::cout << "2. 显示多项式\n";
        std::cout << "3. 计算 f(x)\n";
        std::cout << "4. 求导\n";
        std::cout << "5. 不定积分\n";
        std::cout << "6. 定积分\n";
        std::cout << "0. 返回\n";
        int op=askInt("选择：");

        if(op==0) return;

        if(op==1){
            int n=askInt("输入项数：");
            std::vector<MatCal::Utils::QinJiuShaoNode> nodes;
            for(int i=0;i<n;i++){
                int deg=askInt("  次数：");
                double a=askDouble("  系数：");
                nodes.emplace_back(deg,a);
            }
            poly=QinJiuShao(nodes);
            hasPoly=true;
        }
        else if(!hasPoly){
            std::cout<<"请先构造多项式。\n";
        }
        else if(op==2){
            poly.show();
        }
        else if(op==3){
            double x=askDouble("x=");
            std::cout<<"f("<<x<<") = "<<poly.calculate(x)<<"\n";
        }
        else if(op==4){
            auto d=poly.derivative();
            std::cout<<"f'(x)="; d.show();
            char c; std::cout<<"是否替换当前多项式？(y/n) ";
            std::cin>>c;
            if(c=='y'||c=='Y') poly=d;
        }
        else if(op==5){
            double c=askDouble("积分常数 C=");
            auto I=poly.integral(c);
            std::cout<<"F(x)="; I.show();
        }
        else if(op==6){
            double a=askDouble("a=");
            double b=askDouble("b=");
            std::cout<<"定积分="<<poly.definiteIntegral(a,b)<<"\n";
        }
    }
}

// ========= 插值 =========

void interpolationMenu(){
    while(true){
        std::cout<<"\n--- 插值菜单 ---\n";
        std::cout<<"1. Lagrange\n";
        std::cout<<"2. Newton 差商(不等距)\n";
        std::cout<<"3. Newton 差分(等距)\n";
        std::cout<<"4. Hermite\n";
        std::cout<<"0. 返回\n";
        int op=askInt("选择：");
        if(op==0) return;

        try{
            if(op==1){
                int n=askInt("点数=");
                std::vector<std::pair<double,double>> data(n);
                for(int i=0;i<n;i++){
                    data[i].first = askDouble("x=");
                    data[i].second= askDouble("y=");
                }
                AlgoInsert::LagrangeInsert L(data);
                double X=askDouble("X=");
                std::cout<<"f(X)="<<L.calculate(X)<<"\n";
            }
            else if(op==2){
                int n=askInt("点数=");
                std::vector<std::pair<double,double>> data(n);
                for(int i=0;i<n;i++){
                    data[i].first=askDouble("x=");
                    data[i].second=askDouble("y=");
                }
                AlgoInsert::NewtonInsert_Quotient N(data);
                double X=askDouble("X=");
                std::cout<<"f(X)="<<N.calculate(X)<<"\n";
            }
            else if(op==3){
                double h=askDouble("步长 h=");
                double x0=askDouble("起点 x0=");
                int n=askInt("点数=");
                std::vector<double> ys(n);
                for(int i=0;i<n;i++) ys[i]=askDouble("y="+std::to_string(i)+"=");
                AlgoInsert::NewtonInsert_Finite F(h,x0,ys);
                double X=askDouble("X=");
                std::cout<<"f(X)="<<F.calculate(X)<<"\n";
            }
            else if(op==4){
                int n=askInt("点数=");
                std::vector<double> xs(n), ys(n), dys(n);
                for(int i=0;i<n;i++){
                    xs[i]=askDouble("x=");
                    ys[i]=askDouble("y=");
                    dys[i]=askDouble("y'=");
                }
                AlgoInsert::Hermite H(xs,ys,dys);
                double X=askDouble("X=");
                std::cout<<"f(X)="<<H.calculate(X)<<"\n";
            }
        }
        catch(const std::exception& e){
            std::cout<<"异常："<<e.what()<<"\n";
        }
    }
}

// ========= 迭代求根 =========

void iterationMenu(QinJiuShao &poly, bool hasPoly){
    while(true){
        std::cout<<"\n--- 迭代菜单 ---\n";
        std::cout<<"1. 二分法\n";
        std::cout<<"2. Picard\n";
        std::cout<<"3. Newton\n";
        std::cout<<"4. Newton(下山)\n";
        std::cout<<"5. 割线法(双点)\n";
        std::cout<<"6. 割线法(单点)\n";
        std::cout<<"7. Newton 非线性方程组\n";
        std::cout<<"0. 返回\n";
        int op=askInt("选择：");
        if(op==0) return;

        AlgoIter::Bisection::Function f;
        if(hasPoly) f = poly.toFunction();
        else f = [](double x){ return x*x*x-2*x-5; };

        try{
            if(op==1){
                double a=askDouble("a=");
                double b=askDouble("b=");
                auto r=AlgoIter::Bisection::solveDetailed(f,a,b);
                std::cout<<"root="<<r.root<<"  iter="<<r.iterations<<"  err="<<r.error<<"\n";
            }
            else if(op==2){
                auto phi=[](double x){ return std::cos(x); };
                double x0=askDouble("初值 x0=");
                auto r=AlgoIter::Picard::solveDetailed(phi,x0);
                std::cout<<"root="<<r.root<<" iter="<<r.iterations<<" err="<<r.error<<"\n";
            }
            else if(op==3){
                auto df=AlgoBasics::Derivative::dy_dx;
                double x0=askDouble("初值=");
                auto ff=f;
                auto d =[&](double x){return AlgoBasics::Derivative::dy_dx(ff,x,1e-6);};
                auto r=AlgoIter::Newton::solve(ff,d,x0);
                std::cout<<"root="<<r.root<<" iter="<<r.iterations<<" err="<<r.error<<"\n";
            }
            else if(op==4){
                auto df=AlgoBasics::Derivative::dy_dx;
                double x0=askDouble("初值=");
                auto ff=f;
                auto d =[&](double x){return AlgoBasics::Derivative::dy_dx(ff,x,1e-6);};
                auto r=AlgoIter::Newton::solve_downhill(ff,d,x0);
                std::cout<<"root="<<r.root<<" iter="<<r.iterations<<" err="<<r.error<<"\n";
            }
            else if(op==5){
                double x0=askDouble("x0=");
                double x1=askDouble("x1=");
                auto r=AlgoIter::Secant::solve_two_point(f,x0,x1);
                std::cout<<"root="<<r.root<<" iter="<<r.iterations<<" err="<<r.error<<"\n";
            }
            else if(op==6){
                double x0=askDouble("初值 x0=");
                auto r=AlgoIter::Secant::solve_one_point(f,x0);
                std::cout<<"root="<<r.root<<" iter="<<r.iterations<<" err="<<r.error<<"\n";
            }
            else if(op==7){
                int n=askInt("变量个数 n=");
                std::vector<AlgoIter::NewtonForEquations::Function> funcs(n);
                for(int i=0;i<n;i++){
                    std::cout<<"输入 F"<<i<<"(x1..xn) 的表达方式：此处只能使用示例函数。\n";
                    funcs[i]=[i](const std::vector<double>& v){
                        if(i==0) return v[0]*v[0]+v[1]*v[1]-1.0; // 示例
                        return v[0]-v[1];                     // 示例
                    };
                }
                std::vector<double> x0(n);
                for(int i=0;i<n;i++) x0[i]=askDouble("初值 x["+std::to_string(i)+"]=");

                auto r=AlgoIter::NewtonForEquations::solve(n,funcs,x0);
                std::cout<<"iter="<<r.iterations<<" err="<<r.error<<" conv="<<r.converged<<"\n";
                std::cout<<"root=[";
                for(double v:r.root) std::cout<<v<<" ";
                std::cout<<"]\n";
            }
        }
        catch(const std::exception& e){
            std::cout<<"异常："<<e.what()<<"\n";
        }
    }
}

// ========= 矩阵运算 =========

void matrixMenu(){
    Matrix A;
    bool has=false;

    while(true){
        std::cout<<"\n--- 矩阵菜单 ---\n";
        std::cout<<"1. 输入矩阵 A\n";
        std::cout<<"2. 打印 A\n";
        std::cout<<"3. A + B\n";
        std::cout<<"4. A * B\n";
        std::cout<<"5. 解方程组 Ax=b（列主元）\n";
        std::cout<<"0. 返回\n";

        int op=askInt("选择：");
        if(op==0) return;

        if(op==1){
            A=inputMatrix("A");
            has=true;
        }
        else if(!has){
            std::cout<<"请先输入 A\n";
        }
        else if(op==2){
            A.show();
        }
        else if(op==3){
            Matrix B=inputMatrix("B");
            auto C=A.add(B);
            C->show();
        }
        else if(op==4){
            Matrix B=inputMatrix("B");
            auto C=A.multiply(B);
            C->show();
        }
        else if(op==5){
            Matrix b=inputMatrix("b");
            auto x=AlgoMatrix::solve_columnElimination(A,b);
            x->show();
        }
    }
}

// ========= 数值求导 =========

void derivativeMenu(QinJiuShao &poly,bool hasPoly){
    auto f = hasPoly ? poly.toFunction()
                     : [](double x){return x*x*x-2*x-5;};

    while(true){
        std::cout<<"\n--- 数值求导 ---\n";
        std::cout<<"1. 前向差分\n";
        std::cout<<"2. 中心差分\n";
        std::cout<<"0. 返回\n";

        int op=askInt("选择：");
        if(op==0) return;

        double x=askDouble("x=");
        if(op==1){
            std::cout<<"dy/dx="<<AlgoBasics::Derivative::dy_dx(f,x)<<"\n";
        }
        else if(op==2){
            std::cout<<"dy/dx_center="<<AlgoBasics::Derivative::dy_dx_center(f,x)<<"\n";
        }
    }
}

// ========= 主程序 =========

int main(){
    SetConsoleOutputCP(CP_UTF8);
    std::cout<<"===== MatCal 命令行程序 =====\n";

    QinJiuShao poly;
    bool hasPoly=false;

    while(true){
        std::cout<<"\n=== 主菜单 ===\n"
                 <<"1. 多项式\n"
                 <<"2. 插值\n"
                 <<"3. 迭代求根\n"
                 <<"4. 矩阵运算\n"
                 <<"5. 数值求导\n"
                 <<"0. 退出\n";

        int op=askInt("选择：");
        if(op==0) break;

        if(op==1) polynomialMenu(poly,hasPoly);
        else if(op==2) interpolationMenu();
        else if(op==3) iterationMenu(poly,hasPoly);
        else if(op==4) matrixMenu();
        else if(op==5) derivativeMenu(poly,hasPoly);
    }

    return 0;
}
