#ifndef QINJIUSHAO_HPP
#define QINJIUSHAO_HPP
/*
秦九韶算法
latex or md:
P_n(x)=a_nx^n+a_{n-1}x^{n-1}+...+a_1x+a_0;  \\
P_n(x)=((a_nx^{n-2}+a_{n-1}x^{n-3}+a_2)x+a_1)x+a_0  \\
... \\
P_n(x)=(((a_{n}x+a_{n-1})x+a_{n-2})x+...+a_1)x+a_0

使用std::vector<QinJiuShaoNode>构建链表，实现多项式求值
Node结点采用{次数，系数}构建
需要传入一个非多元函数构建秦九韶计算式，然后可以进行非多元函数快速计算----O(n)复杂度的乘法
秦九韶计算式由(n，a)定义
本秦九韶类计算快速，修改不便

仅支持n为int类型，且不为负

*/
#include<vector>
#include<algorithm>
#include<iostream>
#include<stdexcept>
#include<string>
#include<functional>
#include<map>
#include<limits>

#include "MatCal/Polynomial/Polynomial.hpp"

namespace MatCal{
    namespace Utils{
        //秦九韶结点，用于秦九韶多项式
        class QinJiuShaoNode;
        //秦九韶多项式，提供计算函数值，返回多项式函数(创建通用函数)等服务
        class QinJiuShao;
    }
}

namespace MatCal::Utils{

class QinJiuShaoNode{
public:
    //构造函数，aa为系数，nn为次数。次数仅支持正整数或0
    QinJiuShaoNode(int nn=0,double aa=0){
        if(nn<0){
            std::string errInfo="次数需要大于等于0的整数!";
            this->a=0;
            this->n=0;
            throw std::invalid_argument(errInfo);
        }
        this->a=aa;
        this->n=nn;
    }
    //析构函数，空实现
    ~QinJiuShaoNode(){}
    //无效结点,当系数小于QinJiuShaoNode::ZERO_THRESHOLD认为无效(1e-12)
    bool isFake()const{
        return std::abs(this->a)<QinJiuShaoNode::ZERO_THRESHOLD;
    }
    //toString方法
    std::string toString()const{
        if(this->isFake()){
            return "0";
        }else if(this->n==0){
            return std::to_string(this->a);
        }else if(this->n==1){
            return std::to_string(this->a)+"X";
        }else{
            return std::to_string(this->a) + "X^" + std::to_string(this->n);
        }
    }
    //比较函数
    bool operator>(const QinJiuShaoNode& node) const {
        return (this->n > node.n) || (this->n == node.n && this->a > node.a);
    }
    bool operator<(const QinJiuShaoNode& node) const {
        return (this->n < node.n) || (this->n == node.n && this->a < node.a);
    }
    bool operator==(const QinJiuShaoNode& node) const {
        return this->n == node.n && this->a == node.a;
    }
    bool operator>=(const QinJiuShaoNode& node) const {
        return !(*this < node);
    }
    bool operator<=(const QinJiuShaoNode& node) const {
        return !(*this > node);
    }
    bool operator!=(const QinJiuShaoNode& node) const {
        return !(*this == node);
    }
public:
    //系数
    double a;
    //次数
    int n;
    //静态常量，定义接近0的阈值
    static constexpr double ZERO_THRESHOLD = 1e-12;
};//class QinJiuShaoNode


class QinJiuShao{
//外界可访问方法
public:
//******************构造和析构*****************
    //构造函数，可传入std::vector<QinJiuShaoNode> parameters,可以为空
    QinJiuShao(std::vector<QinJiuShaoNode> par={}){
        this->parameters = par;
        this->cleanup();//自动快速清理数据，维护多项式结构
    }
    //构造函数，可传入std::vector<std::pair<int,double>> parameters，不可为空
    QinJiuShao(std::vector<std::pair<int,double>>& vec){
        //先转换为QinJiuShaoNode
        std::vector<QinJiuShaoNode> temp;
        for(const auto &val:vec)
            temp.emplace_back(val.first, val.second);
        this->parameters = temp;
        this->cleanup();
    }
    //构造函数,支持初始化列表，如：{{2, 3.0}, {1, 1.0}, {0, 5.0}}
    QinJiuShao(std::initializer_list<std::pair<int, double>> pairs){
        //先转换为QinJiuShaoNode
        std::vector<QinJiuShaoNode> temp;
        for(const auto& pair : pairs)
            temp.emplace_back(pair.first, pair.second);
        this->parameters = temp;
        this->cleanup();
    }
    //构造函数,从系数A和零点序列std::vector<double> zeros构造
    QinJiuShao(int A,std::vector<double>& zeros){
        if(std::abs(A)>QinJiuShaoNode::ZERO_THRESHOLD){
            QinJiuShao res={{0,1}};
            for(std::size_t i = 0; i < zeros.size(); ++i){
                QinJiuShao zero={{1,1},{0,-zeros[i]}};
                res=res * zero;
            }
            res = res * A;
            this->parameters = res.getParameters();
        }else{
            QinJiuShao q;
            this->parameters = q.getParameters();
        }
    }
    //拷贝构造函数
    QinJiuShao(const QinJiuShao& qin){
        this->parameters = qin.getParameters();
        this->cleanup();
    }
    QinJiuShao(QinJiuShao&& qin) noexcept = default;
    QinJiuShao& operator=(const QinJiuShao& qin) = default;
    QinJiuShao& operator=(QinJiuShao&& qin) noexcept = default;
    //析构函数
    ~QinJiuShao(){
        this->clear();
    }
//******************增删改查*************
    //插入新项,n-次数，a-系数,insert会自动维护多项式结构
    void insert(int n,double a){
        //旧逻辑，实现快，很稳定，构造时可以不保证有序。就是慢了点---O(n)
        // if(std::abs(a)<QinJiuShaoNode::ZERO_THRESHOLD) return;//如果系数接近0，直接返回
        // auto it = std::find_if(this->parameters.begin(),this->parameters.end(),
        // [n](const QinJiuShaoNode& node){
        //     return node.n == n;
        // });
        // if(it==this->parameters.end()){
        //     this->parameters.push_back(QinJiuShaoNode(n,a));
        //     this->sortByDegree();
        // }else{
        //     it->a+=a;
        //     if(it->isFake())
        //         this->parameters.erase(it);
        // }

        //新逻辑，现在弃用remove_fake使用cleanup,新的构造函数已经保证有序了，那insert可以用二分。---O(logn)
        if(std::abs(a)<QinJiuShaoNode::ZERO_THRESHOLD) return;
        auto it = std::lower_bound(parameters.begin(), parameters.end(), n,
            [](const QinJiuShaoNode& node, int degree){
                return node.n > degree;
            });
        if(it == parameters.end()||it->n!=n){
            parameters.insert(it, QinJiuShaoNode(n, a));
        }else{
            it->a += a;
            if(it->isFake()){
                parameters.erase(it);
            }
        }
    }
    //插入新项,n-次数，a-系数,insert会自动维护多项式结构
    void insert(QinJiuShaoNode node){
        insert(node.n,node.a);
    }
    //移除项，按照次数移除
    bool remove(int n) {
        auto it = std::find_if(this->parameters.begin(),this->parameters.end(),
        [n](const QinJiuShaoNode& node) {
            return node.n == n;
        });
        if(it != this->parameters.end()){
            this->parameters.erase(it);
            return true;
        }
        return false;
    }
    //移除项，按照传入函数
    bool remove_if(std::function<bool(const QinJiuShaoNode&)> _pFunc){
        auto it = std::find_if(this->parameters.begin(),this->parameters.end(),_pFunc);
        if(it != this->parameters.end()){
            this->parameters.erase(it);
            return true;
        }
        return false;
    }
    //秦九韶计算，可传入参数X，计算秦九韶多项式
    double calculate(double x)const{
        return this->toPolynomial().evaluate(x);
    }
    //显示当前秦九韶多项式
    void show()const{
        if(this->parameters.empty()){
            std::cout<<"0"<< std::endl;
            return;
        }
        for(auto it=this->parameters.begin();it!=this->parameters.end();++it){
            if(it->isFake()) continue;
            if(it!=this->parameters.begin()&&it->a>0){
                std::cout<< " + ";
            }
            std::cout<<it->toString();
        }
        std::cout<<std::endl;
    }
    // toString方法 - 基于show方法的逻辑
    std::string toString() const {
        if(this->parameters.empty()){
            return "0";
        }
        std::string result;
        for(auto it=this->parameters.begin();it!=this->parameters.end();++it){
            if(it->isFake()) continue;
            if(it!=this->parameters.begin()&&it->a>0){
                result += " + ";
            }
            result += it->toString();
        }
        return result;
    }
    //获取所有参数,不可修改
    const std::vector<QinJiuShaoNode>& getParameters() const {
        return this->parameters;
    }
    //将多项式转换为拥有系数副本的函数对象，便于数值计算使用
    std::function<double(double)> toFunction() const {
        return this->toPolynomial().to_function();
    }
//************************实用方法*******************
    //获取最高次数
    int getHighestDegree()const{
        if(this->parameters.empty()){
            throw std::runtime_error("Empty polynomial");
        }
        return this->parameters.front().n;
    }
    //获取指定次数的系数
    double getParameter(int n) const {
        auto it = std::find_if(this->parameters.begin(), this->parameters.end(),
            [n](const QinJiuShaoNode& node){
                return node.n == n;
            });
        return (it != this->parameters.end()) ? it->a : 0.0;
    }
    //获取项数
    size_t size()const{
        return this->parameters.size();
    }
    //判断是否为空
    bool empty()const{
        return this->parameters.empty();
    }
    //清空多项式
    void clear(){
        this->parameters.clear();
    }
//************************式级运算********************
    //多项式加法,返回新的秦九韶多项式对象，不改变原来的秦九韶对象
    QinJiuShao operator+(const QinJiuShao& other)const{
        return QinJiuShao::fromPolynomial(this->toPolynomial() + other.toPolynomial());
    }
    //多项式减法,返回新的秦九韶多项式对象，不改变原来的秦九韶对象
    QinJiuShao operator-(const QinJiuShao& other)const{
        return QinJiuShao::fromPolynomial(this->toPolynomial() - other.toPolynomial());
    }
    //多项式乘法，复杂度较高，请谨慎使用---O(n^2)
    QinJiuShao operator*(const QinJiuShao& other) const {
        return QinJiuShao::fromPolynomial(this->toPolynomial() * other.toPolynomial());
    }
    //标量乘法
    QinJiuShao operator*(double scalar)const{
        return QinJiuShao::fromPolynomial(this->toPolynomial() * scalar);
    }
    //标量乘法的友元函数，支持 scalar * polynomial
    friend QinJiuShao operator*(double scalar, const QinJiuShao& poly){
        return poly*scalar;
    }
    //标量除法
    QinJiuShao operator/(double scalar)const{
        return QinJiuShao::fromPolynomial(this->toPolynomial() / scalar);
    }
    //求导运算
    QinJiuShao derivative()const{
        return QinJiuShao::fromPolynomial(this->toPolynomial().derivative());
    }
    //积分运算，请指定积分常数
    QinJiuShao integral(double constant=0)const{
        return QinJiuShao::fromPolynomial(this->toPolynomial().integral(constant));
    }
    //定积分
    double definiteIntegral(double a,double b)const{
        return this->toPolynomial().definite_integral(a, b);
    }

    //求内积,默认权系数为1
    double product(QinJiuShao&other, double a, double b){
        QinJiuShao tmp = *this * other;
        return tmp.definiteIntegral(a,b);
    }


//内部维护方法
private:
    MatCal::Polynomial::Polynomial toPolynomial() const {
        std::vector<MatCal::Polynomial::Polynomial::term> terms;
        terms.reserve(this->parameters.size());
        for(const auto& node : this->parameters){
            if(node.n < 0){
                throw std::invalid_argument("negative polynomial degree");
            }
            terms.emplace_back(static_cast<std::size_t>(node.n), node.a);
        }
        return MatCal::Polynomial::Polynomial::from_terms(terms);
    }

    static QinJiuShao fromPolynomial(const MatCal::Polynomial::Polynomial& polynomial) {
        std::vector<QinJiuShaoNode> nodes;
        for(const auto& [degree, coefficient] : polynomial.terms_descending()){
            if(degree > static_cast<std::size_t>(std::numeric_limits<int>::max())){
                throw std::length_error("legacy polynomial degree exceeds int range");
            }
            nodes.emplace_back(static_cast<int>(degree), coefficient);
        }
        return QinJiuShao(nodes);
    }

    //内部维护，在构造时使用。强制清理所有无效结点并强制排序
    void cleanup(){
        if (this->parameters.empty()) return;
        //使用map自动合并同类项并按键排序
        std::map<int, double, std::greater<int>> coeff_map;//降序排列
        //合并同类项
        for(const auto& node : this->parameters)
                coeff_map[node.n] += node.a;
        //重建vector，过滤合并后接近0的项
        this->parameters.clear();
        for(const auto& [degree, coeff] : coeff_map){
            if (std::abs(coeff) >= QinJiuShaoNode::ZERO_THRESHOLD) {
                this->parameters.emplace_back(degree, coeff);
            }
        }
    }
    //内部维护的秦九韶次数降序排列函数，用户不关心,现在由于已经改了构造结构，这个不需要了
    [[deprecated("这个没啥用了")]]
    void sortByDegree(){
        std::sort(this->parameters.begin(),this->parameters.end(), 
                 [](const QinJiuShaoNode& a, const QinJiuShaoNode& b) {
                     return a.n > b.n; //按次数降序排列
                 });
    }
    //内部维护的无效结点剔除函数，用户不关心。现在没用了
    [[deprecated("用cleanup()方法更好，这个方法是老早之前为了构造方便用的")]]
    void remove_fake(){
        this->parameters.erase(
            std::remove_if(this->parameters.begin(), this->parameters.end(),
                [](const QinJiuShaoNode& node){
                    return node.isFake();
                }),
            this->parameters.end()
        );
        std::sort(this->parameters.begin(), this->parameters.end(),
                 [](const QinJiuShaoNode& a, const QinJiuShaoNode& b) {
                     return a.n > b.n;
                 });//重新排序
    }
//内部维护属性
private:
    //多项式参数
    std::vector<QinJiuShaoNode> parameters;
};//class QinJiuShao
}//namespace MatCal::Utils
#endif//本文件的头文件定义
