#ifndef MATRIX_HPP
#define MATRIX_HPP

/*

矩阵工具类，实现了常用矩阵


*/

#include<iostream>
#include<memory>
#include<string>
#include<vector>
#include<stdexcept>
#include<iomanip>
#include<algorithm>
#include<random>
#include<ctime>
#include<sstream>
#include<cstring>
#include<tuple>

namespace MatCal{
    namespace Utils{
        //抽象矩阵
        class AbstractMatrix;
        //普通/稠密矩阵
        class Matrix;
        //稀疏矩阵
        class SparseMatrix;
        //抽象三角矩阵
        class AbstractTriangularMatrix;
        //特殊矩阵:上三角矩阵
        class UpperTriangularMatrix;
        //特殊矩阵:下三角矩阵
        class LowerTriangularMatrix;
        //特殊矩阵:三对角矩阵
        class TridiagonalMatrix;
    }
    namespace Algorithm{
        //初等行变换
        std::unique_ptr<Utils::AbstractMatrix> swapRows(const Utils::AbstractMatrix& A, int r1, int r2);
        std::unique_ptr<Utils::AbstractMatrix> scaleRow(const Utils::AbstractMatrix& A, int r, double scalar);
        std::unique_ptr<Utils::AbstractMatrix> addScaledRow(const Utils::AbstractMatrix& A, int r_target, int r_source, double scalar);

        //列初等变换
        std::unique_ptr<Utils::AbstractMatrix> swapCols(const Utils::AbstractMatrix& A, int c1, int c2);
        std::unique_ptr<Utils::AbstractMatrix> scaleCol(const Utils::AbstractMatrix& A, int c, double scalar);
        std::unique_ptr<Utils::AbstractMatrix> addScaledCol(const Utils::AbstractMatrix& A, int c_target, int c_source, double scalar);
    
    }
}

namespace MatCal::Utils{

//抽象矩阵，这是一个抽象类
class AbstractMatrix{
protected:
    int rows;   //行数
    int cols;   //列数
public:
    //抽象的父类构造，矩阵构造需要行数和列数
    AbstractMatrix(int row=0,int col=0):rows(row),cols(col){}
    //虚析构
    virtual ~AbstractMatrix()=0;
public:
//抽象矩阵公有方法
    //获取行数
    virtual int getRows()const{return rows;}
    //获取列数
    virtual int getCols()const{return cols;}
    //设置矩阵大小（纯虚函数）
    virtual void resize(int newRows,int newCols)=0;
    //获取元素（纯虚函数）
    virtual double get(int row,int col)const=0;
    //设置元素（纯虚函数）
    virtual void set(int row, int col,double value)=0;
    //展示矩阵
    virtual void show()const{
        std::cout<<"Matrix:"<<rows<<"x"<<cols<<"\n";
    }
    //toString方法
    virtual std::string toString()const{
        return "AbstractMatrix "+std::to_string(rows)+"x"+std::to_string(cols);
    }
    //转换普通矩阵方法,返回智能指针（纯虚函数）
    virtual std::unique_ptr<AbstractMatrix> toNormalMatrix()const=0;
    //矩阵加法（纯虚函数）
    virtual std::unique_ptr<AbstractMatrix> add(const AbstractMatrix& other)const=0;
    //矩阵乘法（纯虚函数）
    virtual std::unique_ptr<AbstractMatrix> multiply(const AbstractMatrix& other)const=0;
    virtual std::unique_ptr<AbstractMatrix> operator*(const AbstractMatrix& other)const=0;
    //标量乘法（纯虚函数）
    virtual std::unique_ptr<AbstractMatrix> scalarMultiply(double scalar)const=0;
    virtual std::unique_ptr<AbstractMatrix> operator*(double scalar)const=0;
    //转置矩阵（纯虚函数）
    //virtual std::unique_ptr<AbstractMatrix> transpose()const=0;
    //判断是否是方阵
    virtual bool isSquare()const{return rows==cols;}
};//class AbstractMatrix
AbstractMatrix::~AbstractMatrix() = default;//纯许析构实现

//普通矩阵
class Matrix:public AbstractMatrix{
private:    //私有属性
    std::vector<std::vector<double>> data;

public:     //构造函数和析构
    //普通构造
    Matrix(int row=0,int col=0):AbstractMatrix(row,col){
        if(rows<0||cols<0){
            throw std::invalid_argument("Matrix dimensions cannot be negative");
        }
        resize(rows, cols);
    }

    //从二维向量构造
    Matrix(const std::vector<std::vector<double>>& input):AbstractMatrix(input.size(),input.empty()?0:input[0].size()){
        //验证所有行长度一致
        for(size_t i = 1;i<input.size();++i){
            if(input[i].size()!=input[0].size()){
                throw std::invalid_argument("All rows should have the same length");
            }
        }
        data=input;
    }

    //从二维数组构造
    Matrix(const double** input,int row,int col):AbstractMatrix(row, col){
        if(input == nullptr)
            throw std::invalid_argument("Input array cannot be null");
        // 验证所有行指针非空
        for(int i=0;i<row; ++i){
            if(input[i]==nullptr){
                throw std::invalid_argument("Row pointer cannot be null");
            }
        }
        // 初始化
        data.resize(row);
        for(int i=0;i<row;++i){
            data[i].resize(col);
            for(int j=0;j<col;++j){
                data[i][j] = input[i][j];
            }
        }
    }

    //拷贝构造函数
    Matrix(const Matrix& other):AbstractMatrix(other.getRows(),other.getCols()){
        data=other.getData();
    }

    Matrix(std::initializer_list<std::initializer_list<double>> init)
        : AbstractMatrix(init.size(), init.size() > 0 ? init.begin()->size() : 0) {
        
        data.resize(rows);
        int i = 0;
        for (const auto& row : init) {
            data[i].resize(cols);
            int j = 0;
            for (const auto& elem : row) {
                data[i][j] = elem;
                ++j;
            }
            // 填充剩余列为零
            for (; j < cols; ++j) {
                data[i][j] = 0.0;
            }
            ++i;
        }
    }

    //从父级指针构造
    Matrix(AbstractMatrix& A){
        std::unique_ptr<AbstractMatrix> temp = A.toNormalMatrix();
        Matrix* copy_child = dynamic_cast<Matrix*>(temp.get());
        if(!copy_child)
            throw std::bad_cast();
        this->data = copy_child->data;
    }
    Matrix& operator=(AbstractMatrix& A){
        if(this == &A) return *this;
        std::unique_ptr<AbstractMatrix> temp = A.toNormalMatrix();
        Matrix* copy_child = dynamic_cast<Matrix*>(temp.get());
        if(!copy_child)
            throw std::bad_cast();
        this->data.clear();
        this->data = copy_child->data;
        return *this;
    }

    //赋值运算符
    // Matrix& operator=(const std::vector<std::vector<double>>& input){
    //     //验证所有行长度一致
    //     for(size_t i = 1;i<input.size();++i){
    //         if(input[i].size()!=input[0].size()){
    //             throw std::invalid_argument("All rows should have the same length");
    //         }
    //     }
    //     data=input;
    //     return *this;
    // }
    Matrix& operator=(const Matrix& other) {
        if(this != &other){
            rows = other.getRows();
            cols = other.getCols();
            data = other.getData();
        }
        return *this;
    }

    //析构
    ~Matrix(){}
public:     //抽象方法实现
    //重新设置大小，请注意，如果设置比之前小了会丢失数据，增大按照右下角扩大
    void resize(int newRows,int newCols)override{
        if(newRows < 0||newCols < 0)
            throw std::invalid_argument("Matrix dimensions cannot be negative");
        this->rows = newRows;
        this->cols = newCols;
        data.resize(this->rows);
        for(int i = 0; i < this->rows; ++i){
            data[i].resize(this->cols, 0.0);
        }
    }

    //获取元素
    double get(int row,int col)const override{
        if(row < 0 || row >= rows || col < 0 || col >= cols)
            throw std::out_of_range("Matrix index out of range");
        return data[row][col];
    }

    //设置元素
    void set(int row, int col, double value)override{
        if(row < 0 || row >= rows || col < 0 || col >= cols)
            throw std::out_of_range("Matrix index out of range");
        data[row][col] = value;
    }

    //展示普通矩阵
    void show() const override{
        std::cout << "\nMatrix:"<<rows <<"x"<<cols<<":\n";
        for(int i = 0; i < rows; ++i){
            for(int j = 0; j < cols; ++j)
                std::cout<<std::setw(10)<<std::setprecision(6)<<data[i][j] << " ";
            std::cout<<"\n";
        }
    }

    //toString方法，转化成Json格式
    std::string toString()const override{
        std::ostringstream json;
        json << "{\n";
        json << "  \"type\": \"Matrix\",\n";
        json << "  \"rows\": " << rows << ",\n";
        json << "  \"cols\": " << cols << ",\n";
        json << "  \"isSquare\": " << (isSquare() ? "true" : "false") << ",\n";
        json << "  \"data\": [\n";
        
        for (int i = 0; i < rows; ++i) {
            json << "    [";
            for (int j = 0; j < cols; ++j) {
                json << std::fixed << std::setprecision(8) << data[i][j];
                if (j < cols - 1) {
                    json << ", ";
                }
            }
            json << "]";
            if (i < rows - 1) {
                json << ",";
            }
            json << "\n";
        }
        
        json << "  ]\n";
        json << "}";
        return json.str();
    }

    //转换成普通矩阵,普通矩阵直接返回自己,你不用调用了
    std::unique_ptr<AbstractMatrix> toNormalMatrix() const override{
        //修改，需要拷贝一份，不能直接返回*this
        Matrix copy=*this;
        return std::make_unique<Matrix>(copy);
    }

    //内置的矩阵加法，时间复杂度O(n^2)*(+)
    std::unique_ptr<AbstractMatrix> add(const AbstractMatrix& other)const override{
        if(rows != other.getRows() || cols != other.getCols())
            throw std::invalid_argument("Matrix dimensions not match!");
        auto result = std::make_unique<Matrix>(rows, cols);
        for(int i = 0; i < rows; ++i){
            for(int j = 0; j < cols; ++j){
                (*result)[i][j] = data[i][j] + other.get(i, j);
            }
        }
        return result;
    }

    //内置的矩阵乘法（右乘）,时间复杂度O(n^3)*(*)
    std::unique_ptr<AbstractMatrix> multiply(const AbstractMatrix& other) const override {
        if (cols != other.getRows())
            throw std::invalid_argument("Matrix dimensions not match! try (m x n),(n x k)");

        auto result = std::make_unique<Matrix>(rows, other.getCols());
        for(int i = 0; i < rows; ++i){
            for(int j = 0; j < other.getCols(); ++j){
                double sum = 0.0;
                for(int k = 0; k < cols; ++k){
                    sum += data[i][k] * other.get(k, j);
                }
                (*result)[i][j] = sum;
            }
        }
        return result;
    }
    std::unique_ptr<AbstractMatrix> operator*(const AbstractMatrix& other) const override{
        return this->multiply(other);
    }

    //内置的标量乘法,时间复杂度O(n^2)*(*)
    std::unique_ptr<AbstractMatrix> scalarMultiply(double scalar) const override{
        auto result = std::make_unique<Matrix>(rows, cols);
        for(int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                (*result)[i][j] = data[i][j] * scalar;
            }
        }
        return result;
    }
    std::unique_ptr<AbstractMatrix> operator*(double scalar) const override{
        return this->scalarMultiply(scalar);
    }
    friend std::unique_ptr<AbstractMatrix> operator*(double scalar,const AbstractMatrix&matrix){
        return matrix.scalarMultiply(scalar);
    }

    //内置的转置运算
    std::unique_ptr<AbstractMatrix> transpose()const {
        auto result = std::make_unique<Matrix>(cols, rows);
        for(int i = 0; i < rows; ++i){
            for(int j = 0; j < cols; ++j){
                (*result)[j][i]=data[i][j];
            }
        }
        return result;
    }

public:     //Matrix实例对象独有的的实用方法
    //重载运算符[],可以快速访问元素
    std::vector<double>& operator[](int row){
        if(row < 0 || row >= rows) {
            throw std::out_of_range("Matrix row index out of range");
        }
        return data[row];
    }

    //重载运算符[],可以快速访问元素,const版本
    const std::vector<double>& operator[](int row) const {
        if(row < 0 || row >= rows) {
            throw std::out_of_range("Matrix row index out of range");
        }
        return data[row];
    }

    //获得内部数据引用
    const std::vector<std::vector<double>>& getData()const{
        return this->data;
    }

public:     //Matrix类级别静态方法
    //获取单位矩阵
    static Matrix identity(int size){
        Matrix result(size, size);
        for (int i = 0; i < size; ++i) {
            result[i][i] = 1.0;
        }
        return result;
    }

    //获得零矩阵
    static Matrix zeros(int rows, int cols){
        return Matrix(rows, cols);
    }

    //获得一矩阵
    static Matrix ones(int rows, int cols){
        Matrix result(rows, cols);
        for (int i = 0;i < rows; ++i)
            for (int j = 0;j < cols;++j)
                result[i][j]=1.0;
        return result;
    }

};//class Matrix

//普通稀疏矩阵
class SparseMatrix:public AbstractMatrix{
public:     //内部类
    struct Element{
        int row;
        int col;
        double value;
        Element(int r=0,int c=0,double v=0.0):row(r),col(c),value(v){}
    };
private:    //私有属性
    std::vector<Element> elements;  //非零元素
    int nonZeroCount;               //非零元素个数
public:     //构造和析构
    //普通构造
    SparseMatrix(int rows=0, int cols=0):AbstractMatrix(rows, cols),nonZeroCount(0){}
    
    //从普通矩阵构造
    SparseMatrix(const Matrix& matrix):AbstractMatrix(matrix.getRows(),matrix.getCols()),nonZeroCount(0){
        auto data = matrix.getData();
        for(int i = 0; i < rows; ++i){
            for (int j = 0; j < cols; ++j) {
                if(data[i][j]!=0.0){
                    elements.emplace_back(i, j, data[i][j]);
                    nonZeroCount++;
                }
            }
        }
    }

    //从三元组构造 (行, 列, 值)
    SparseMatrix(int rows,int cols,const std::vector<std::tuple<int,int,double>>& triplets):AbstractMatrix(rows,cols),nonZeroCount(0){
        for(const auto& triplet:triplets){
            int r = std::get<0>(triplet);
            int c = std::get<1>(triplet);
            double v = std::get<2>(triplet);
            if (r >= 0 && r < rows && c >= 0 && c < cols && v != 0.0) {
                elements.emplace_back(r, c, v);
                nonZeroCount++;
            }
        }
        //按行列排序以便高效访问
        std::sort(elements.begin(), elements.end(), 
            [](const Element& a, const Element& b){
                return (a.row == b.row) ? (a.col < b.col) : (a.row < b.row);
            });
    }

    //析构
    ~SparseMatrix(){}

public:     //抽象方法实现
    //重设大小
    void resize(int newRows, int newCols) override{
        //移除超出新范围的非零元素
        elements.erase(
            std::remove_if(elements.begin(), elements.end(),
                [newRows, newCols](const Element& e){
                    return e.row >= newRows || e.col >= newCols;
                }),
            elements.end()
        );
        nonZeroCount = elements.size();
        rows = newRows;
        cols = newCols;
    }

    //get方法
    double get(int row, int col) const override{
        if(row < 0 || row >= rows || col < 0 || col >= cols)
            throw std::out_of_range("SparseMatrix index out of range");
        auto it = std::find_if(elements.begin(), elements.end(),
            [row, col](const Element& e){
                return e.row == row && e.col == col;
            });
        return (it != elements.end()) ? it->value : 0.0;
    }

    //set方法
    void set(int row, int col, double value) override{
        if (row < 0 || row >= rows || col < 0 || col >= cols)
            throw std::out_of_range("SparseMatrix index out of range");
        
        //查找是否已存在
        auto it = find(row,col);
        if(it != elements.end()){   //已经存在
            if(value == 0.0){       //零元素直接删掉
                elements.erase(it);
                nonZeroCount--;
            }else{
                it->value = value;
            }
        }else if(value != 0.0){     //不存在，插入非零数值
            elements.emplace_back(row, col, value);
            nonZeroCount++;
            inner_sort();
        }
    }

    //展示矩阵
    void show() const override{
        std::cout << "SparseMatrix: "<<rows<<"x"<< cols
                  <<" ; (NonZero Elements: " << nonZeroCount << "):\n";
        
        auto dense =this->toNormalMatrix();//转换为稠密矩阵显示
        dense->show();
    }

    //toString方法，转换成Json
    std::string toString() const override {
        std::ostringstream json;
        json << "{\n";
        json << "  \"type\": \"SparseMatrix\",\n";
        json << "  \"rows\": " << rows << ",\n";
        json << "  \"cols\": " << cols << ",\n";
        json << "  \"nonZeroCount\": " << nonZeroCount << ",\n";
        json << "  \"sparsity\": " << std::fixed << std::setprecision(8) 
             << (1.0 - static_cast<double>(nonZeroCount) / (rows * cols)) << ",\n";
        json << "  \"elements\": [\n";
        
        for (size_t i = 0; i < elements.size(); ++i) {
            const auto& elem = elements[i];
            json << "    {\"row\": " << elem.row 
                 << ", \"col\": " << elem.col 
                 << ", \"value\": " << std::fixed << std::setprecision(8) << elem.value << "}";
            if (i < elements.size() - 1) {
                json << ",";
            }
            json << "\n";
        }
        
        json << "  ]\n";
        json << "}";
        return json.str();
    }

    //转换成普通矩阵(稠密矩阵)
    std::unique_ptr<AbstractMatrix> toNormalMatrix() const override {
        auto dense = std::make_unique<Matrix>(rows, cols);
        for (const auto& elem : elements)
            dense->set(elem.row, elem.col, elem.value);
        return dense;
    }

    //内置的稀疏矩阵加法
    std::unique_ptr<AbstractMatrix> add(const AbstractMatrix& other) const override {
        if (rows != other.getRows() || cols != other.getCols())
            throw std::invalid_argument("Matrix dimensions must match!");
        auto result = std::make_unique<SparseMatrix>(rows, cols);
        //拷贝当前矩阵的非零元素
        for(const auto& elem : elements)
            result->set(elem.row, elem.col, elem.value);
        //开加!
        for(int i = 0; i < rows; ++i){
            for(int j = 0; j < cols; ++j){
                double otherVal = other.get(i, j);
                if(otherVal != 0.0)
                    result->set(i, j, result->get(i, j) + otherVal);//交给set方法来维护结构，这里只管放进去
            }
        }
        return result;
    }

    //内置的稀疏矩阵乘法
    std::unique_ptr<AbstractMatrix> multiply(const AbstractMatrix& other) const override{
        if (cols != other.getRows())
            throw std::invalid_argument("Matrix dimensions not match! try (m x n),(n x k)");
        auto result = std::make_unique<SparseMatrix>(rows, other.getCols());
        for(const auto& elem : elements){
            int i = elem.row;
            int k = elem.col;
            double val = elem.value;
            
            for(int j = 0; j < other.getCols(); ++j){
                double otherVal = other.get(k, j);
                if(otherVal != 0.0){
                    double current = result->get(i, j);
                    result->set(i, j, current + val * otherVal);//要逐步累加
                }
            }
        }
        return result;
    }
    std::unique_ptr<AbstractMatrix> operator*(const AbstractMatrix& other) const override{
        return this->multiply(other);
    }

    //内置的标量乘法
    std::unique_ptr<AbstractMatrix> scalarMultiply(double scalar) const override{
        auto result = std::make_unique<SparseMatrix>(rows, cols);
        if (scalar == 0.0)//0直接返回
            return result;
        for (const auto& elem : elements)
            result->set(elem.row, elem.col, elem.value*scalar);
        return result;
    }
    std::unique_ptr<AbstractMatrix> operator*(double scalar) const override{
        return this->scalarMultiply(scalar);
    }

    //内置的矩阵转置
    std::unique_ptr<AbstractMatrix> transpose() const {
        auto result = std::make_unique<SparseMatrix>(cols, rows);
        for (const auto& elem : elements)
            result->set(elem.col, elem.row, elem.value);
        return result;
    }

public:     //SparseMatrix实例对象独有方法
    //内部排序
    void inner_sort(){
        std::sort(elements.begin(), elements.end(),
                [](const Element& a, const Element& b) {
                    return (a.row == b.row) ? (a.col < b.col) : (a.row < b.row);
                });
    }

    //查找元素,返回迭代器
    std::vector<Element>::iterator find(int row,int col){
        return std::find_if(elements.begin(), elements.end(),
            [row, col](const Element& e){
                return e.row == row && e.col == col;
            });
    }

    //获取非零元素个数
    int getNonZeroCount() const{
        return nonZeroCount;
    }

    //获取稀疏度，1-非0个数除以总个数(0-1之间，1表示完全稀疏)
    double getSparsity() const {
        if(rows==0||cols==0)
            throw std::invalid_argument("no elements!");
        return 1.0 - static_cast<double>(nonZeroCount) / (rows * cols);
    }

    //获取所有非零元素的三元组
    std::vector<std::tuple<int, int, double>> getTriplets()const{
        std::vector<std::tuple<int, int, double>> triplets;
        triplets.reserve(elements.size());
        for (const auto& elem : elements)
            triplets.emplace_back(elem.row, elem.col, elem.value);
        return triplets;
    }

    //获取内部elements的引用
    const std::vector<Element>& getElements() const{
        return this->elements;
    }

public:     //SparseMatrix类级静态方法
    //创建单位稀疏矩阵
    static SparseMatrix identity(int size){
        return SparseMatrix(Matrix::identity(size));
    }

    //创建对角稀疏矩阵
    static SparseMatrix diagonal(const std::vector<double>& diag){
        int size = diag.size();
        SparseMatrix result(size, size);
        for(int i = 0; i < size; ++i)
            if (diag[i] != 0.0)
                result.set(i, i, diag[i]);
        return result;
    }

    //创建随机稀疏矩阵
    static SparseMatrix random(int rows, int cols, double density){
        SparseMatrix result(rows, cols);
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> valueDist(0.0, 1.0);
        std::uniform_int_distribution<int> rowDist(0, rows-1);
        std::uniform_int_distribution<int> colDist(0, cols-1);
        
        int targetNonZero = static_cast<int>(rows * cols * density);
        for (int i = 0; i < targetNonZero; ++i)
            result.set(rowDist(gen), colDist(gen), valueDist(gen));
        return result;
    }

};//class SparseMatrix


//三角矩阵抽象基类
class AbstractTriangularMatrix:public AbstractMatrix{
protected:  //私有属性
    std::vector<double> data;
protected:  //纯虚函数 - 子类必须实现
    //索引
    virtual int calculateIndex(int row, int col) const = 0;
    //判断区域
    virtual bool isInStoredRegion(int row, int col) const = 0;
    //创建空的同类型三角矩阵    
    virtual std::unique_ptr<AbstractTriangularMatrix> createEmptyClone() const = 0;
public:     //构造析构
    //普通构造
    AbstractTriangularMatrix(int size=0) : AbstractMatrix(size, size){
        int elementCount = size * (size + 1) / 2;
        data.resize(elementCount, 0.0);
    }

    //析构
    virtual ~AbstractTriangularMatrix()=0;

public:     //实现父类AbstractMatrix的方法
    //重设大小
    void resize(int newRows, int newCols) override{
        if (newRows != newCols)
            throw std::invalid_argument("Triangular matrix must be square");
        rows = newRows;
        cols = newCols;
        int elementCount = newRows * (newRows + 1) / 2;
        data.resize(elementCount, 0.0);
    }

    //get方法
    double get(int row, int col) const override{
        if(row < 0 || row >= rows || col < 0 || col >= cols)
            throw std::out_of_range("Matrix index out of range");
        if(!isInStoredRegion(row, col))
            return 0.0;//非存储区域为零
        int index = calculateIndex(row, col);
        return data[index];
    }

    //set方法
    void set(int row, int col, double value) override {
        if(row < 0 || row >= rows || col < 0 || col >= cols)
            throw std::out_of_range("Matrix index out of range");
        if(!isInStoredRegion(row, col) && value != 0.0)
            throw std::invalid_argument("Cannot set non-zero value outside triangular region");
        if(isInStoredRegion(row, col)){
            int index = calculateIndex(row, col);
            data[index] = value;
        }
    }

    //转换成普通矩阵
    std::unique_ptr<AbstractMatrix> toNormalMatrix() const override {
        auto dense = std::make_unique<Matrix>(rows, cols);
        for(int i = 0; i < rows; ++i){
            for (int j = 0; j < cols; ++j)
                dense->set(i, j, get(i, j));
        }
        return dense;
    }

    //展示矩阵
    void show()const override{
        //转换成普通矩阵再show
        auto matrix=this->toNormalMatrix();
        matrix->show();
    }

    //内置加法
    std::unique_ptr<AbstractMatrix> add(const AbstractMatrix& other) const override {
        if (rows != other.getRows() || cols != other.getCols())
            throw std::invalid_argument("Matrix dimensions must match!");
        // 返回普通矩阵，因为相加结果可能不再是三角矩阵
        auto result = std::make_unique<Matrix>(rows, cols);
        for(int i = 0; i < rows; ++i)
            for (int j = 0; j < cols; ++j)
                result->set(i, j, get(i, j) + other.get(i, j));
        return result;
    }

    //内置乘法
    std::unique_ptr<AbstractMatrix> multiply(const AbstractMatrix& other) const override {
        if (cols != other.getRows())
            throw std::invalid_argument("Matrix dimensions not match! try (m x n),(n x k)");
        //三角矩阵相乘结果可能不再是三角矩阵
        auto result = std::make_unique<Matrix>(rows, other.getCols());
        for(int i = 0; i < rows; ++i){
            for(int j = 0; j < other.getCols(); ++j){
                double sum = 0.0;
                for (int k = 0; k < cols; ++k)
                    sum += get(i, k) * other.get(k, j);
                result->set(i, j, sum);
            }
        }
        return result;
    }
    std::unique_ptr<AbstractMatrix> operator*(const AbstractMatrix& other) const override{
        return this->multiply(other);
    }

    //内置标量乘法
    std::unique_ptr<AbstractMatrix> scalarMultiply(double scalar) const override{
        //标量乘法保持三角性质
        auto result = createEmptyClone();
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                if (isInStoredRegion(i, j)) {
                    result->set(i, j, get(i, j) * scalar);
                }
            }
        }
        return result;
    }
    std::unique_ptr<AbstractMatrix> operator*(double scalar) const override{
        return this->scalarMultiply(scalar);
    }


    //转置，继续往下抛
    //std::unique_ptr<AbstractMatrix> transpose() const override=0;

public:     //三角矩阵实例对象特有的方法
    const std::vector<double> getData()const{
        return this->data;
    }

    //计算行列式（对角线相乘）
    double determinant() const{
        double det = 1.0;
        for (int i = 0; i < rows; ++i)
            det *= get(i, i);
        return det;
    }

    //是否奇异
    bool isSingular() const {
        for(int i = 0; i < rows; ++i)
            if (get(i, i) == 0.0)
                return true;
        return false;
    }



};//class AbstractTriangularMatrix
AbstractTriangularMatrix::~AbstractTriangularMatrix()=default;//纯许析构实现

//特殊矩阵:上三角矩阵
class UpperTriangularMatrix : public AbstractTriangularMatrix {
protected:      //父类虚函数实现
    //计算顺序vector里的索引，行优先存储上三角元素:
    int calculateIndex(int row, int col) const override {
        return col*(col+1)/2+row;
    }
    //是否在上三角区域
    bool isInStoredRegion(int row, int col) const override {
        return row <= col;
    }
    //空的结构拷贝
    std::unique_ptr<AbstractTriangularMatrix> createEmptyClone() const override {
        return std::make_unique<UpperTriangularMatrix>(rows);
    }

public:     //构造析构
    UpperTriangularMatrix(int size=0) : AbstractTriangularMatrix(size) {}
    UpperTriangularMatrix(UpperTriangularMatrix&u){
        this->data=u.getData();
    }
    UpperTriangularMatrix(const Matrix& dense) : AbstractTriangularMatrix(dense.getRows()){
        if(!dense.isSquare())
            throw std::invalid_argument("Triangular matrix must be square");
        for(int i = 0; i < rows; ++i){
            for(int j = 0; j < cols; ++j){
                double v = dense.get(i, j);
                // if(!isInStoredRegion(i, j) && v != 0.0)
                //     throw std::invalid_argument("Matrix does not match triangular type for UpperTriangularMatrix");
                if(isInStoredRegion(i, j))
                    set(i, j, v);
            }
        }
    }


public:     //实现父类AbstractTriangularMatrix未实现的爷爷类的方法
    //toString方法,Json
    std::string toString() const override {
        std::ostringstream json;
        json << "{\n";
        json << "  \"type\": \"UpperTriangularMatrix\",\n";
        json << "  \"size\": " << rows << ",\n";
        json << "  \"storageSize\": " << data.size() << ",\n";
        json << "  \"determinant\": " << std::fixed << std::setprecision(8) << determinant() << ",\n";
        json << "  \"data\": [\n";
        
        for (int i = 0; i < rows; ++i) {
            json << "    [";
            for (int j = 0; j < cols; ++j) {
                json << std::fixed << std::setprecision(8) << get(i, j);
                if (j < cols - 1) json << ", ";
            }
            json << "]";
            if (i < rows - 1) json << ",";
            json << "\n";
        }
        
        json << "  ]\n";
        json << "}";
        return json.str();
    }

    //转置，上三角矩阵转置后变成下三角矩阵
    std::unique_ptr<LowerTriangularMatrix> transpose()const ;

public:    // 上三角矩阵特有的方法
    //求解Ax=b(b支持多组解)
    std::unique_ptr<AbstractMatrix> solve(const AbstractMatrix& b) const {
        // 回代法求解上三角系统 Ux = b （Ux=y）
        if (b.getRows() != rows) {
            throw std::invalid_argument("Dimension mismatch");
        }
        
        auto x = std::make_unique<Matrix>(rows, b.getCols());
        
        for (int k = 0; k < b.getCols(); ++k) {     //第一层循环，遍历不同的b_1,b_2 ......
            for (int i = rows - 1; i >= 0; --i) {   //第二层循环,从最后一行开始往上遍历
                double sum = 0.0;
                for (int j = i + 1; j < cols; ++j) {
                    sum += get(i, j) * x->get(j, k);
                }
                double diag = get(i, i);
                if (diag == 0.0) {
                    throw std::runtime_error("Matrix is singular");
                }
                x->set(i, k, (b.get(i, k) - sum) / diag);
            }
        }
        
        return x;
    }

    //上三角矩阵与上三角矩阵相乘仍然是上三角矩阵
    std::unique_ptr<AbstractMatrix> multiplyUpper(const UpperTriangularMatrix& other) const {
        if (cols != other.getRows())
            throw std::invalid_argument("Matrix dimensions incompatible for multiplication");
        auto result = std::make_unique<UpperTriangularMatrix>(rows);
        for (int i = 0; i < rows; ++i) {
            for (int j = i; j < other.getCols(); ++j){  //只计算上三角部分
                double sum = 0.0;
                for (int k = i; k <= j; ++k)            //优化计算范围
                    if (k < cols && k < other.getRows())
                        sum += get(i, k) * other.get(k, j);
                result->set(i, j, sum);
            }
        }
        return result;
    }
};//class UpperTriangularMatrix

//特殊矩阵:下三角矩阵
class LowerTriangularMatrix : public AbstractTriangularMatrix {
protected:      //父类虚函数实现
    //计算顺序vector里的索引，行优先存储上三角元素
    int calculateIndex(int row, int col) const override {
         return row * (row + 1) / 2 + col;
    }
    //是否在上三角区域
    bool isInStoredRegion(int row, int col) const override {
        return row >= col;
    }
    //空的结构拷贝
    std::unique_ptr<AbstractTriangularMatrix> createEmptyClone() const override {
        return std::make_unique<LowerTriangularMatrix>(rows);
    }

public:     //构造析构
    LowerTriangularMatrix(int size=0) : AbstractTriangularMatrix(size) {}
    LowerTriangularMatrix(LowerTriangularMatrix&l){
        this->data=l.getData();
    }
    LowerTriangularMatrix(const LowerTriangularMatrix& l){
        this->data=l.getData();
    }
    LowerTriangularMatrix(const Matrix& dense) : AbstractTriangularMatrix(dense.getRows()) {
        if(!dense.isSquare())
            throw std::invalid_argument("Triangular matrix must be square");
        for(int i = 0; i < rows; ++i){
            for(int j = 0; j < cols; ++j){
                double v = dense.get(i, j);
                // if(!isInStoredRegion(i, j) && v != 0.0)
                //     throw std::invalid_argument("Matrix does not match triangular type for LowerTriangularMatrix");
                if(isInStoredRegion(i, j))
                    set(i, j, v);
            }
        }
    }


public:     //实现父类AbstractTriangularMatrix未实现的爷爷类的方法
    //toString方法,Json
    std::string toString() const override {
        std::ostringstream json;
        json << "{\n";
        json << "  \"type\": \"LowerTriangularMatrix\",\n";
        json << "  \"size\": " << rows << ",\n";
        json << "  \"storageSize\": " << data.size() << ",\n";
        json << "  \"determinant\": " << std::fixed << std::setprecision(8) << determinant() << ",\n";
        json << "  \"data\": [\n";
        
        for (int i = 0; i < rows; ++i) {
            json << "    [";
            for (int j = 0; j < cols; ++j) {
                json << std::fixed << std::setprecision(8) << get(i, j);
                if (j < cols - 1) json << ", ";
            }
            json << "]";
            if (i < rows - 1) json << ",";
            json << "\n";
        }
        
        json << "  ]\n";
        json << "}";
        return json.str();
    }

    //转置，下三角矩阵转置后变成上三角矩阵
    std::unique_ptr<UpperTriangularMatrix> transpose() const ;

public:    //下三角矩阵特有的方法
    //下三角的求解，常用于LU-->Ly=b,同样支持多组解
    std::unique_ptr<AbstractMatrix> solve(const AbstractMatrix& b) const{
        // 回代法求解上三角系统 Ly=b
        if (b.getRows() != rows) {
            throw std::invalid_argument("Dimension mismatch");
        }
        
        auto y = std::make_unique<Matrix>(rows, b.getCols());
        
        for (int k = 0; k < b.getCols(); ++k) {     //第一层循环，遍历不同的b_1,b_2 ......
            for (int i = 0; i < rows; ++i) {   //第二层循环,从第一行开始往下遍历
                double sum = 0.0;
                for (int j = 0; j < i; ++j) {
                    sum += get(i, j) * y->get(j, k);
                }
                double diag = get(i, i);
                if (diag == 0.0) {
                    throw std::runtime_error("Matrix is singular");
                }
                y->set(i, k, (b.get(i, k) - sum) / diag);
            }
        }
        
        return y;
    }

    //下三角矩阵与下三角矩阵相乘仍然是下三角矩阵
    std::unique_ptr<LowerTriangularMatrix> multiplyLower(const LowerTriangularMatrix& other) const {
        if (cols != other.getRows())
            throw std::invalid_argument("Matrix dimensions incompatible for multiplication");
        auto result = std::make_unique<LowerTriangularMatrix>(rows);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j <= i; ++j) {              //只计算下三角部分
                double sum = 0.0;
                for (int k = j; k <= i; ++k)            //优化计算范围
                    if (k < cols && k < other.getRows())
                        sum += get(i, k) * other.get(k, j);
                result->set(i, j, sum);
            }
        }
        return result;
    }
};//class LowerTriangularMatrix


//特殊矩阵:三对角矩阵
class TridiagonalMatrix : public AbstractMatrix {
protected:
    std::vector<double> l;   // 下对角线，长度 n-1
    std::vector<double> d;   // 主对角线，长度 n
    std::vector<double> u;   // 上对角线，长度 n-1
public:
    /*-------- 构造 / 析构 --------*/
    TridiagonalMatrix(int n = 0) : AbstractMatrix(n, n) {
        if (n < 0) throw std::invalid_argument("n must be non-negative");
        resize(n, n);
    }

    // 从三条对角线构造
    TridiagonalMatrix(const std::vector<double>& lower,
                      const std::vector<double>& diag,
                      const std::vector<double>& upper)
        : AbstractMatrix(diag.size(), diag.size()),
          l(lower), d(diag), u(upper) {
        int n = static_cast<int>(d.size());
        if (static_cast<int>(l.size()) != n - 1 ||
            static_cast<int>(u.size()) != n - 1)
            throw std::invalid_argument("Diagonal size mismatch");
        rows = cols = n;
    }

    // 从普通 Matrix 提取三对角部分
    explicit TridiagonalMatrix(const Matrix& mat)
        : AbstractMatrix(mat.getRows(), mat.getCols()) {
        if (!mat.isSquare())
            throw std::invalid_argument("Tridiagonal matrix must be square");
        int n = mat.getRows();
        l.resize(n - 1);
        d.resize(n);
        u.resize(n - 1);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                double v = mat.get(i, j);
                if (i == j) d[i] = v;
                else if (i == j + 1) l[i - 1] = v;
                else if (i + 1 == j) u[j - 1] = v;
                else if (std::abs(v) > 1e-14)
                    throw std::invalid_argument("Matrix is not tridiagonal");
            }
        }
    }

    virtual ~TridiagonalMatrix() = default;

    /*-------- AbstractMatrix 接口 --------*/
    void resize(int newRows, int newCols) override {
        if (newRows != newCols)
            throw std::invalid_argument("Tridiagonal matrix must be square");
        if (newRows < 0)
            throw std::invalid_argument("Size must be non-negative");
        rows = cols = newRows;
        int n = newRows;
        l.resize(n - 1);
        d.resize(n);
        u.resize(n - 1);
        std::fill(l.begin(), l.end(), 0.0);
        std::fill(d.begin(), d.end(), 0.0);
        std::fill(u.begin(), u.end(), 0.0);
    }

    double get(int row, int col) const override {
        checkIndex(row, col);
        if (row == col) return d[row];
        if (row == col + 1) return l[row - 1];
        if (row + 1 == col) return u[col - 1];
        return 0.0;
    }

    void set(int row, int col, double value) override {
        checkIndex(row, col);
        if (row == col) d[row] = value;
        else if (row == col + 1) l[row - 1] = value;
        else if (row + 1 == col) u[col - 1] = value;
        else if (std::abs(value) > 1e-14)
            throw std::invalid_argument("Cannot set non-zero outside tridiagonal");
    }

    std::unique_ptr<AbstractMatrix> toNormalMatrix() const override {
        auto dense = std::make_unique<Matrix>(rows, cols);
        for (int i = 0; i < rows; ++i) {
            dense->set(i, i, d[i]);
            if (i > 0) dense->set(i, i - 1, l[i - 1]);
            if (i + 1 < cols) dense->set(i, i + 1, u[i]);
        }
        return dense;
    }

    void show() const override {
        auto tmp = toNormalMatrix();
        tmp->show();
    }

    std::string toString() const override {
        std::ostringstream os;
        os << "{\n"
           << "  \"type\": \"TridiagonalMatrix\",\n"
           << "  \"size\": " << rows << ",\n"
           << "  \"lower\": " << vec2str(l) << ",\n"
           << "  \"diag\": " << vec2str(d) << ",\n"
           << "  \"upper\": " << vec2str(u) << "\n"
           << "}";
        return os.str();
    }

    /*-------- 运算接口（返回普通 Matrix） --------*/
    std::unique_ptr<AbstractMatrix> add(const AbstractMatrix& other) const override {
        if (rows != other.getRows() || cols != other.getCols())
            throw std::invalid_argument("Dimension mismatch");
        auto res = std::make_unique<Matrix>(rows, cols);
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < cols; ++j)
                res->set(i, j, get(i, j) + other.get(i, j));
        return res;
    }

    std::unique_ptr<AbstractMatrix> multiply(const AbstractMatrix& other) const override {
        if (cols != other.getRows())
            throw std::invalid_argument("Dimension mismatch");
        auto res = std::make_unique<Matrix>(rows, other.getCols());
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < other.getCols(); ++j) {
                double s = 0.0;
                // 三对角每行最多 3 个非零
                for (int k = std::max(0, i - 1); k <= std::min(cols - 1, i + 1); ++k)
                    s += get(i, k) * other.get(k, j);
                res->set(i, j, s);
            }
        }
        return res;
    }

    std::unique_ptr<AbstractMatrix> operator*(const AbstractMatrix& other) const override {
        return multiply(other);
    }

    std::unique_ptr<AbstractMatrix> scalarMultiply(double scalar) const override {
        auto T = std::make_unique<TridiagonalMatrix>(rows);
        for (size_t i = 0; i < l.size(); ++i) T->l[i] = l[i] * scalar;
        for (size_t i = 0; i < d.size(); ++i) T->d[i] = d[i] * scalar;
        for (size_t i = 0; i < u.size(); ++i) T->u[i] = u[i] * scalar;
        return T;
    }

    std::unique_ptr<AbstractMatrix> operator*(double scalar) const override {
        return scalarMultiply(scalar);
    }

    /*-------- 三对角专用求解 --------*/
    // Thomas(chase) 算法解  T * X = B，X 与 B 同尺寸
    std::unique_ptr<AbstractMatrix> solve(const AbstractMatrix& B) const {
        if (B.getRows() != rows)
            throw std::invalid_argument("Dimension not match");
        int n = rows;
        int nrhs = B.getCols();

        //拷贝一份
        std::vector<double> dl = l;
        std::vector<double> dd = d;
        std::vector<double> du = u;
        auto X=B.toNormalMatrix();

        //前消,构建LU
        for (int i = 1; i < n; ++i) {
            if (std::abs(dd[i - 1]) < 1e-14)
                throw std::runtime_error("Zero pivot in Thomas(chase)");
            dl[i - 1] = dl[i - 1] / dd[i - 1];
            dd[i] -= dl[i - 1] * du[i - 1];
            for (int k = 0; k < nrhs; ++k)
                X->set(i, k, X->get(i, k) - dl[i - 1] * X->get(i - 1, k));
        }
        if (std::abs(dd.back()) < 1e-14)
            throw std::runtime_error("Singular matrix in Thomas(chase)");

        //回代,求解LU(不直接使用LU算法，那样就体现不出追赶法降低复杂度的特点了)
        for (int k = 0; k < nrhs; ++k)
            X->set(n - 1, k, X->get(n - 1, k) / dd.back());
        for (int i = n - 2; i >= 0; --i)
            for (int k = 0; k < nrhs; ++k) {
                double tmp = (X->get(i, k) - du[i] * X->get(i + 1, k)) / dd[i];
                X->set(i, k, tmp);
            }
        return X;
    }

    /*-------- 工具 --------*/
    const std::vector<double>& lower() const { return l; }
    const std::vector<double>& diag()  const { return d; }
    const std::vector<double>& upper() const { return u; }

private:
    void checkIndex(int r, int c) const {
        if (r < 0 || r >= rows || c < 0 || c >= cols)
            throw std::out_of_range("Index out of range");
    }
    static std::string vec2str(const std::vector<double>& v) {
        std::ostringstream o;
        o << "[";
        for (size_t i = 0; i < v.size(); ++i) {
            o << std::fixed << std::setprecision(8) << v[i];
            if (i + 1 < v.size()) o << ", ";
        }
        o << "]";
        return o.str();
    }
};

}//namespace MatCal::Utils

namespace MatCal::Utils {
// UpperTriangularMatrix::transpose() 的实现
std::unique_ptr<LowerTriangularMatrix> UpperTriangularMatrix::transpose() const {
    auto result = std::make_unique<LowerTriangularMatrix>(rows);
    for (int i = 0; i < rows; ++i)
        for (int j = i; j < cols; ++j)
            result->set(j, i, get(i, j)); 
    return result; // 移出后不需要 std::move
}

// LowerTriangularMatrix::transpose() 的实现
std::unique_ptr<UpperTriangularMatrix> LowerTriangularMatrix::transpose() const {
    auto result = std::make_unique<UpperTriangularMatrix>(rows);
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j <= i; ++j)
            result->set(j, i, get(i, j)); 
    return result; // 移出后不需要 std::move
}
}
namespace MatCal::Algorithm::Matrix{
    using Utils::AbstractMatrix;
    using Utils::Matrix;
    using Utils::SparseMatrix;

//**************初等变换*********************
    //检查边界的辅助函数
    inline void checkRowBounds(int rows, int r) {
        if (r < 0 || r >= rows) {
            throw std::out_of_range("Row index out of bounds during elementary transformation.");
        }
    }
    inline void checkColBounds(int cols, int c) {
        if (c < 0 || c >= cols) {
            throw std::out_of_range("Column index out of bounds during elementary transformation.");
        }
    }
    
    //行初等变换1: 行交换 R_r1 <-> R_r2
    void swapRows(AbstractMatrix& A, int r1, int r2) {
        checkRowBounds(A.getRows(), r1);
        checkRowBounds(A.getRows(), r2);
        if (r1 == r2) return;

        if (auto dense = dynamic_cast<Matrix*>(&A)) {
            std::swap((*dense)[r1], (*dense)[r2]);
            return;
        }

        auto result = A.toNormalMatrix();
        auto dense = dynamic_cast<Matrix*>(result.get());
        if (!dense)
            throw std::runtime_error("toNormalMatrix() did not return a Matrix!");
        std::swap((*dense)[r1], (*dense)[r2]);
        A = *dense;
    }


    //行初等变换2: 行乘法 R_r <- scalar * R_r
    void scaleRow(AbstractMatrix& A, int r, double scalar) {
        if (scalar == 1.0) return;
        checkRowBounds(A.getRows(), r);

        if (auto dense = dynamic_cast<Matrix*>(&A)) {
            for (int j = 0; j < dense->getCols(); ++j)
                (*dense)[r][j] *= scalar;
            return;
        }

        auto result = A.toNormalMatrix();
        auto dense = dynamic_cast<Matrix*>(result.get());
        if (!dense)
            throw std::runtime_error("toNormalMatrix() did not return a Matrix!");

        for (int j = 0; j < dense->getCols(); ++j)
            (*dense)[r][j] *= scalar;

        A = *dense;
    }


    //行初等变换3: 行加法 R_r_target <- R_r_target + scalar * R_r_source
    void addScaledRow(AbstractMatrix& A, int r_target, int r_source, double scalar) {
        if (scalar == 0.0 || r_target == r_source) return;
        checkRowBounds(A.getRows(), r_target);
        checkRowBounds(A.getRows(), r_source);

        if (auto dense = dynamic_cast<Matrix*>(&A)) {
            for (int j = 0; j < dense->getCols(); ++j)
                (*dense)[r_target][j] += scalar * (*dense)[r_source][j];
            return;
        }

        auto result = A.toNormalMatrix();
        auto dense = dynamic_cast<Matrix*>(result.get());
        if (!dense)
            throw std::runtime_error("toNormalMatrix() did not return a Matrix!");

        for (int j = 0; j < dense->getCols(); ++j)
            (*dense)[r_target][j] += scalar * (*dense)[r_source][j];

        A = *dense;
    }


    //列初等变换1: 列交换 C_c1 <-> C_c2
    void swapCols(AbstractMatrix& A, int c1, int c2) {
        if (c1 == c2) return;
        checkColBounds(A.getCols(), c1);
        checkColBounds(A.getCols(), c2);

        if (auto dense = dynamic_cast<Matrix*>(&A)) {
            for (int i = 0; i < dense->getRows(); ++i)
                std::swap((*dense)[i][c1], (*dense)[i][c2]);
            return;
        }

        auto result = A.toNormalMatrix();
        auto dense = dynamic_cast<Matrix*>(result.get());
        if (!dense)
            throw std::runtime_error("toNormalMatrix() did not return a Matrix!");

        for (int i = 0; i < dense->getRows(); ++i)
            std::swap((*dense)[i][c1], (*dense)[i][c2]);

        A = *dense;
    }


    //列初等变换2: 列乘法 C_c <- scalar * C_c
    void scaleCol(AbstractMatrix& A, int c, double scalar) {
        if (scalar == 1.0) return;
        checkColBounds(A.getCols(), c);

        if (auto dense = dynamic_cast<Matrix*>(&A)) {
            for (int i = 0; i < dense->getRows(); ++i)
                (*dense)[i][c] *= scalar;
            return;
        }

        auto result = A.toNormalMatrix();
        auto dense = dynamic_cast<Matrix*>(result.get());
        if (!dense)
            throw std::runtime_error("toNormalMatrix() did not return a Matrix!");

        for (int i = 0; i < dense->getRows(); ++i)
            (*dense)[i][c] *= scalar;

        A = *dense;
    }


    //列初等变换3: 列加法 C_c_target <- C_c_target + scalar * C_c_source
    void addScaledCol(AbstractMatrix& A, int c_target, int c_source, double scalar) {
        if (scalar == 0.0 || c_target == c_source) return;
        checkColBounds(A.getCols(), c_target);
        checkColBounds(A.getCols(), c_source);

        if (auto dense = dynamic_cast<Matrix*>(&A)) {
            for (int i = 0; i < dense->getRows(); ++i)
                (*dense)[i][c_target] += scalar * (*dense)[i][c_source];
            return;
        }

        auto result = A.toNormalMatrix();
        auto dense = dynamic_cast<Matrix*>(result.get());
        if (!dense)
            throw std::runtime_error("toNormalMatrix() did not return a Matrix!");

        for (int i = 0; i < dense->getRows(); ++i)
            (*dense)[i][c_target] += scalar * (*dense)[i][c_source];

        A = *dense;
    }

    //按照记录的行交换操作对矩阵 M 进行多次行变换
    // 参数 reverse=true 表示按相反顺序执行，用于“还原”矩阵
    void applyRowSwaps(AbstractMatrix& M, const std::vector<std::pair<int,int>>& swaps, bool reverse = false) {
        int rows = M.getRows();

        if (swaps.empty()) return;

        // 决定迭代方向
        if (!reverse) {
            // 正向重放：重现高斯消元过程中的行交换
            for (const auto& [r1, r2] : swaps) {
                checkRowBounds(rows, r1);
                checkRowBounds(rows, r2);
                if (r1 != r2) swapRows(M, r1, r2);
            }
        } else {
            // 逆向重放：按相反顺序还原为原始顺序
            for (auto it = swaps.rbegin(); it != swaps.rend(); ++it) {
                int r1 = it->first;
                int r2 = it->second;
                checkRowBounds(rows, r1);
                checkRowBounds(rows, r2);
                if (r1 != r2) swapRows(M, r1, r2);
            }
        }
    }


//求解矩阵范数
    //一范数（列范数）
    double norm_one(AbstractMatrix&matrix){
        //需要先转化成普通矩阵
        if(typeid(matrix)!=typeid(Matrix)){
            auto ret=matrix.toNormalMatrix();
            auto dense=dynamic_cast<Matrix*>(ret.get());

            double max=0;
            for(int j=0;j<dense->getCols();++j){
                double sum=0;
                for(int i=0;i<dense->getRows();++i)
                    sum+=std::abs(dense->get(i,j));
                if(sum>max)
                    max=sum;
            }
            return max;
        }else{
            double max=0;
            for(int j=0;j<matrix.getCols();++j){
                double sum=0;
                for(int i=0;i<matrix.getRows();++i)
                    sum+=std::abs(matrix.get(i,j));
                if(sum>max)
                    max=sum;
            }
            return max;
        }
    }

    //无穷范数（行范数）
    double norm_infinite(AbstractMatrix&matrix){
        //需要先转化成普通矩阵
        if(typeid(matrix)!=typeid(Matrix)){
            auto ret=matrix.toNormalMatrix();
            auto dense=dynamic_cast<Matrix*>(ret.get());

            double max=0;
            for(int i=0;i<dense->getRows();++i){
                double sum=0;
                for(int j=0;j<dense->getCols();++j)
                    sum+=std::abs(dense->get(i,j));
                if(sum>max)
                    max=sum;
            }
            return max;
        }else{
            double max=0;
            for(int i=0;i<matrix.getRows();++i){
                double sum=0;
                for(int j=0;j<matrix.getCols();++j)
                    sum+=std::abs(matrix.get(i,j));
                if(sum>max)
                    max=sum;
            }
            return max;
        }
    }

    // F范数
    double norm_Frobenius(AbstractMatrix& matrix){
        double sum = 0.0;
        if(auto dense = dynamic_cast<Matrix*>(&matrix)){
            for(int i=0;i<dense->getRows();++i){
                for (int j = 0; j < dense->getCols(); ++j){
                    double val = dense->get(i, j);
                    sum+=val*val;
                }
            }
        }else{
            auto converted = matrix.toNormalMatrix();
            auto dense_ptr = dynamic_cast<Matrix*>(converted.get());
            if (!dense_ptr)
                throw std::runtime_error("Failed to convert to normal matrix");
            
            for(int i=0;i<dense_ptr->getRows();++i){
                for(int j=0;j<dense_ptr->getCols();++j){
                    double val = dense_ptr->get(i,j);
                    sum+=val*val;
                }
            }
        }
        return std::sqrt(sum);
    }




//**********************线性方程组求解****************

    //列主元消去法求解Ax=b
    //列主元解一次方程(支持一次多组求解),求解交给上三角矩阵来做，所以不会除以对角元素
    std::unique_ptr<AbstractMatrix> solve_columnElimination(AbstractMatrix& A,AbstractMatrix& b){
        //一共rows轮次
        int rows=A.getRows();
        int col_a=A.getCols();
        int col_b=b.getCols();
        auto A_copy=A.toNormalMatrix();
        auto b_copy=b.toNormalMatrix();
        for(int i=0;i<rows;++i){
            //每一轮选取col==row的列，在row>=i区域进行行初等变换
            int index=i;
            for(int j=i;j<rows;++j)
                if(std::abs(A_copy->get(j, i))>std::abs(A_copy->get(index, i)))
                    index=j;
            swapRows(*A_copy,index,i);
            swapRows(*b_copy,index,i);
            //高斯消元
            for(int j = i + 1; j < rows; ++j){  //从当前行之后的所有行
                double factor = A_copy->get(j, i) / A_copy->get(i, i);
                for (int k = i; k < col_a; ++k)
                    A_copy->set(j,k,A_copy->get(j,k)-factor*A_copy->get(i, k));   //更新矩阵 A
                for(int cb=0;cb<col_b;++cb)
                    b_copy->set(j,cb,b_copy->get(j,cb)-factor*b_copy->get(i,cb)); //更新向量 b
            }
        }
        auto dense = dynamic_cast<Matrix*>(A_copy.get());
        MatCal::Utils::UpperTriangularMatrix U(*dense);
        return U.solve(*b_copy);
    }

    //列主元解多次方程，返回的是变换后的上三角矩阵和初等变换的记录,不会对A变换（如果要变换，请applyRowSwaps）
    std::pair<std::unique_ptr<AbstractMatrix>,std::vector<std::pair<int,int>>> columnElimination_Transformation(AbstractMatrix& A) {
        int rows = A.getRows();
        int colsA = A.getCols();
        std::vector<std::pair<int,int>> swaps; //记录所有行交换

        auto A_copy=A.toNormalMatrix();

        for (int i = 0; i < rows; ++i) {
            // 1. 选择当前列的主元
            int pivot = i;
            for (int j = i + 1; j < rows; ++j)
                if (std::abs(A_copy->get(j, i)) > std::abs(A_copy->get(pivot, i)))
                    pivot = j;

            // 2. 行交换（如果必要）
            if (pivot != i) {
                swapRows(*A_copy, pivot, i);
                swaps.emplace_back(i, pivot);
            }

            double pivotVal = A_copy->get(i, i);
            if (std::abs(pivotVal) < 1e-12)
                throw std::runtime_error("Matrix is singular or nearly singular");

            // 3. 高斯消元
            for (int j = i + 1; j < rows; ++j) {
                double factor = A_copy->get(j, i) / pivotVal;
                for (int k = i; k < colsA; ++k)
                    A_copy->set(j, k, A_copy->get(j, k) - factor * A_copy->get(i, k));
            }
        }

        // 4. 提取上三角矩阵
        auto dense = dynamic_cast<Matrix*>(&(*A_copy));
        auto U = std::make_unique<MatCal::Utils::UpperTriangularMatrix>(*dense);

        return {std::move(U), swaps};
    }

    //借助列主元的快速行列式求解
    double determinant(AbstractMatrix& A){
        if(!A.isSquare())
            throw std::invalid_argument("Matrix should be squre!");
        auto [upper,log]=columnElimination_Transformation(A);
        int n=A.getCols();
        double ret=1.0;
        for(int i=0;i<n;++i)
            ret*=upper->get(i,i);
        return ret;
    }

    //LU分解结果类
    struct LUresult{
        MatCal::Utils::LowerTriangularMatrix L;
        MatCal::Utils::UpperTriangularMatrix U;
        LUresult(MatCal::Utils::LowerTriangularMatrix&l,MatCal::Utils::UpperTriangularMatrix&u){
            L=l;
            U=u;
        }
        std::unique_ptr<AbstractMatrix> solve(AbstractMatrix&x){
            return U.solve(*L.solve(x));
        }
    };
    //LU分解,返回LU分解结果类
    LUresult LU_Decompose(AbstractMatrix& A, double eps = 1e-12) {
        int rows = A.getRows();
        int cols = A.getCols();
        
        if (rows != cols) {
            throw std::invalid_argument("LU Decomposition requires a square matrix.");
        }

        // 初始化 L 和 U
        // L (LowerTriangularMatrix) 默认下三角存储，非存储区域为 0
        // U (UpperTriangularMatrix) 默认上三角存储，非存储区域为 0
        auto L = std::make_unique<MatCal::Utils::LowerTriangularMatrix>(rows);
        auto U = std::make_unique<MatCal::Utils::UpperTriangularMatrix>(rows);
        
        // 对 A 进行 LU 分解
        for (int i = 0; i < rows; ++i) {
            // 1. 设置 L 矩阵的对角线 L_ii = 1.0 (Doolittle)
            L->set(i, i, 1.0); 

            // 2. 检查主元，避免除零
            double pivot = A.get(i, i);
            if (std::abs(pivot) < eps) {
                throw std::runtime_error("Matrix is singular or near singular, cannot perform A=LU decomposition without pivoting.");
            }

            // 3. 计算 U 的元素 (U_ij = A_ij - sum(L_ik * U_kj))
            // 在此处，我们将 A 的当前行作为 U 的行，但需要减去 L U 的点积贡献
            for (int j = i; j < cols; ++j) { // 从对角线开始，计算 U 的行
                double sum = 0.0;
                for (int k = 0; k < i; ++k) {
                    sum += L->get(i, k) * U->get(k, j);
                }
                double u_ij = A.get(i, j) - sum;
                U->set(i, j, u_ij);
            }

            // 4. 计算 L 的元素 (L_ji = (A_ji - sum(L_jk * U_ki)) / U_ii)
            // 从下一行开始，计算 L 的列
            for (int j = i + 1; j < rows; ++j) {
                double sum = 0.0;
                for (int k = 0; k < i; ++k) {
                    sum += L->get(j, k) * U->get(k, i);
                }
                double l_ji = (A.get(j, i) - sum) / U->get(i, i);
                L->set(j, i, l_ji);
            }
        }

        // 返回 LU 结果
        return LUresult(*L, *U);
    }

    //迭代法求解线性方程组
    //迭代结果类,包含  {找到的解，迭代次数，最终误差，是否收敛}，考虑到解是向量，存储使用Matrix，不返回迭代序列（空间占用太大）
    struct Iteration_Result{
        Matrix root;        //找到的解
        int iterations;     //迭代次数
        double error;       //最终误差,默认采用一范数（每一列都要成立）
        bool converged;     //是否收敛
    };
    //Jacobi法
    Iteration_Result Jacobi(AbstractMatrix&A,AbstractMatrix&b,double epsilon=1e-6,int max_iterations=100){
        //转换为普通矩阵
        auto converted_A = A.toNormalMatrix();
        auto converted_b = b.toNormalMatrix();
        auto dense_A = dynamic_cast<Matrix*>(converted_A.get());
        auto dense_b = dynamic_cast<Matrix*>(converted_b.get());
        if (!dense_A||!dense_b)
            throw std::runtime_error("Failed to convert matrices to dense format");
        
        int n = dense_A->getRows();
        int sets=dense_b->getCols();
        //提取对角矩阵D
        auto D = std::make_unique<Matrix>(n, n);
        for(int i = 0; i < n; ++i)
            D->set(i,i,dense_A->get(i,i));
        
        //构造L+U (A = D + L + U, 所以 L+U = A - D)
        auto LU = std::make_unique<Matrix>(n, n);
        for (int i = 0; i < n; ++i){
            for(int j = 0; j < n; ++j){
                if(i == j){
                    LU->set(i,j,0.0);
                }else{
                    LU->set(i,j, dense_A->get(i, j));
                }
            }
        }
        
        //计算迭代矩阵 J = -D^(-1)(L+U)  以及g=D^(-1)b
        auto J = std::make_unique<Matrix>(n,n);
        auto g = std::make_unique<Matrix>(n,sets);
        for(int i = 0; i < n; ++i){
            double diag_val = D->get(i, i);
            if(std::abs(diag_val) < 1e-12)
                throw std::runtime_error("Zero diagonal element found, Jacobi method fails");
            for(int j = 0; j < n; ++j){
                J->set(i, j, -LU->get(i, j) / diag_val);
            }
            for(int k=0;k<sets;++k)
                g->set(i,k, dense_b->get(i,k) / diag_val);
        }

        //初始化迭代
        auto x = std::make_unique<Matrix>(n, sets);  // 初始解向量，全零
        Iteration_Result result;

        //Jacobi迭代: x^(k+1) = J * x^(k) + g
        for(int iteration=1;iteration<=max_iterations;++iteration){
            //将AbstractMatrix转换为Matrix
            auto x_new_abstract = J->multiply(*x)->add(*g);
            auto x_new_matrix = dynamic_cast<Matrix*>(x_new_abstract.get());
            if (!x_new_matrix)
                throw std::runtime_error("Matrix operation returned unexpected type");
            auto x_new = std::make_unique<Matrix>(*x_new_matrix);
            //计算误差（一范数：每一列的最大绝对误差）
            auto delta_x=x_new->add(*(x->scalarMultiply(-1)));
            double eps=norm_one(*delta_x);
            //更新解
            *x=*x_new;
            //检查收敛
            if(eps < epsilon){
                result.converged=true;
                result.iterations=iteration;
                result.error=eps;
                result.root=*x;
                return result;
            }
        }
        //未收敛
        result.converged=false;
        result.iterations=max_iterations;
        result.error=-1;//表示未收敛
        result.root=*x;
        return result;
    }

    //Gauss-Seidel(不使用矩阵方法了，因为求矩阵的逆不好求且不稳定，与其求逆不如直接对A求逆)
    Iteration_Result Gauss_Seidel(AbstractMatrix&A,AbstractMatrix&b,double epsilon=1e-6,int max_iterations=100){
        //转换为普通矩阵
        auto converted_A = A.toNormalMatrix();
        auto converted_b = b.toNormalMatrix();
        auto dense_A = dynamic_cast<Matrix*>(converted_A.get());
        auto dense_b = dynamic_cast<Matrix*>(converted_b.get());
        if (!dense_A||!dense_b)
            throw std::runtime_error("Failed to convert matrices to dense format");
        
        int n = dense_A->getRows();
        int sets=dense_b->getCols();
       
        //初始化迭代
        auto x = std::make_unique<Matrix>(n, sets);
        Iteration_Result result;

        //开始迭代
        for(int iteration=1;iteration<=max_iterations;++iteration){
            auto x_new=std::make_unique<Matrix>(*x);
            //每一组解
            for(int k=0;k<sets;++k){
                //每个方程
               for(int i=0;i<n;++i){
                    double sum = 0.0;
                    for(int j=0;j<i;++j)
                        sum += dense_A->get(i,j)*x_new->get(j,k);
                    for(int j=i+1;j<n;++j)
                        sum += dense_A->get(i,j)*x->get(j,k);
                    double diag_val = dense_A->get(i,i);
                    if(std::abs(diag_val)<1e-12)
                        throw std::runtime_error("Zero diagonal element found, Gauss-Seidel method fails");
                    double new_val=(dense_b->get(i,k)-sum)/diag_val;
                    x_new->set(i,k, new_val);
                }//for i
            }//for k
            //计算误差
            auto delta_x_abstract = x->add(*(x_new->scalarMultiply(-1)));
            auto delta_x_ptr = dynamic_cast<Matrix*>(delta_x_abstract.get());
            if (!delta_x_ptr)
                throw std::runtime_error("Matrix operation returned unexpected type");
            double eps = norm_one(*delta_x_ptr);
            //更新解
            *x=*x_new;
            //检查收敛
            if(eps < epsilon) {
                result.converged = true;
                result.iterations = iteration;
                result.error = eps;
                result.root = *x;
                return result;
            }

        }//for iteration
        // 未收敛
        result.converged = false;
        result.iterations = max_iterations;
        result.error = -1;
        result.root = *x;
        return result;
    }

    //SOR方法，需要传入omega参数，0<omega<2
    Iteration_Result SOR(AbstractMatrix& A, AbstractMatrix& b, int omega, double epsilon = 1e-6, int max_iterations = 100) {
        if (omega >= 2 || omega <= 0) {
            throw std::invalid_argument("omega should be in (0,2) !");
        }
        
        //转换为普通矩阵
        auto converted_A = A.toNormalMatrix();
        auto converted_b = b.toNormalMatrix();
        auto dense_A = dynamic_cast<Matrix*>(converted_A.get());
        auto dense_b = dynamic_cast<Matrix*>(converted_b.get());
        if (!dense_A || !dense_b)
            throw std::runtime_error("Failed to convert matrices to dense format");

        int n = dense_A->getRows();
        int sets = dense_b->getCols();

        //初始化迭代
        auto x = std::make_unique<Matrix>(n, sets);
        Iteration_Result result;

        //开始迭代
        for (int iteration = 1; iteration <= max_iterations; ++iteration) {
            auto x_new = std::make_unique<Matrix>(*x);
            //每一组解
            for (int k = 0; k < sets; ++k) {
                //每个方程
                for (int i = 0; i < n; ++i) {
                    double sum = 0.0;
                    for (int j = 0; j < i; ++j)
                        sum += dense_A->get(i, j) * x_new->get(j, k);
                    for (int j = i + 1; j < n; ++j)
                        sum += dense_A->get(i, j) * x->get(j, k);
                    double diag_val = dense_A->get(i, i);
                    if (std::abs(diag_val) < 1e-12)
                        throw std::runtime_error("Zero diagonal element found, SOR method fails");
                    double gauss = (dense_b->get(i, k) - sum) / diag_val;//相对于Gauss,只改了这两行
                    double new_val = gauss * omega + (1 - omega) * dense_A->get(i, k);
                    x_new->set(i, k, new_val);
                }//for i
            }//for k
            //计算误差
            auto delta_x_abstract = x->add(*(x_new->scalarMultiply(-1)));
            auto delta_x_ptr = dynamic_cast<Matrix*>(delta_x_abstract.get());
            if (!delta_x_ptr)
                throw std::runtime_error("Matrix operation returned unexpected type");
            double eps = norm_one(*delta_x_ptr);
            //更新解
            *x = *x_new;
            //检查收敛
            if (eps < epsilon) {
                result.converged = true;
                result.iterations = iteration;
                result.error = eps;
                result.root = *x;
                return result;
            }

        }//for iteration
        // 未收敛
        result.converged = false;
        result.iterations = max_iterations;
        result.error = -1;
        result.root = *x;
        return result;
    
    }


//***************求解特征值***************** */

    //幂法结果类
    struct PowerMethod_Result{
        double eigenvalue;
        Matrix eigenvector;
        int iterations;
        PowerMethod_Result(double val,Matrix& vec,int iter = 0):eigenvalue(val),eigenvector(vec),iterations(iter){};
    };
    //规范化幂法
    PowerMethod_Result PowerMethod(AbstractMatrix&A,double eps = 1e-6,int max_iter = 1000){
        if(!A.isSquare())
            throw std::invalid_argument("using power method,matrix should be square!");
        if(A.getCols()==0)
            throw std::invalid_argument("matrix empty!");
        int n = A.getCols();
        //std::function<double(AbstractMatrix& vec)>
        auto _powerMethod_find_mk = [n](AbstractMatrix& vec){
            double mk = vec.get(0,0);
            for(int i=0;i<n;++i){
                if(std::abs(vec.get(i,0)) > std::abs(mk)){
                    mk = vec.get(i,0);
                }
            }
            return mk;
        };

        Matrix raw_vec(n,1);
        for(int i=0;i<n;++i)
            raw_vec.set(i,0,1);
        
        Matrix u = raw_vec;
        Matrix v = raw_vec;
        double m = 0;
        double last_m = 1;
        int iter = 0;
        //主循环
        while(std::abs(m - last_m)>eps && iter < max_iter){
            last_m = m;
            v = *A.multiply(u);
            m = _powerMethod_find_mk(v);
            u = *v.scalarMultiply(1/m);
            ++iter;
        }
        return PowerMethod_Result(m,u,iter);
    }
    //规范化反幂法
    PowerMethod_Result PowerMethod_reverse(AbstractMatrix&A,double near_num = 0,double eps = 1e-6,int max_iter = 1000){
        if(!A.isSquare())
            throw std::invalid_argument("using power method,matrix should be square!");
        if(A.getCols()==0)
            throw std::invalid_argument("matrix empty!");
        int n = A.getCols();
        auto _powerMethod_find_mk = [n](AbstractMatrix& vec){
            double mk = vec.get(0,0);
            for(int i=0;i<n;++i){
                if(std::abs(vec.get(i,0)) > std::abs(mk)){
                    mk = vec.get(i,0);
                }
            }
            return mk;
        };
        
        Matrix raw_vec(n,1);
        for(int i=0;i<n;++i)
            raw_vec.set(i,0,1);
        
        Matrix I = Matrix::identity(A.getCols());
        Matrix B(n,n);
        B = *(A.add( *I.scalarMultiply(-near_num)));
        auto lu = LU_Decompose(B);

        Matrix u = raw_vec;
        Matrix v = raw_vec;
        double m = 0;
        double last_m = 1;
        int iter = 0;
        //主循环
        while(std::abs(m - last_m)>eps && iter < max_iter){
            last_m = m;
            v = *lu.solve(u);
            m = _powerMethod_find_mk(v);
            u = *v.scalarMultiply(1/m);
            ++iter;
        }
        return PowerMethod_Result(near_num+1/m,u,iter);
    }

}//namespace MatCal::Algorithm::Matrix


#endif//MATRIX_HPP