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
        //class TridiagonalMatrix;
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
    AbstractMatrix(int row,int col):rows(row),cols(col){}
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
    //标量乘法（纯虚函数）
    virtual std::unique_ptr<AbstractMatrix> scalarMultiply(double scalar)const=0;
    //转置矩阵（纯虚函数）
    virtual std::unique_ptr<AbstractMatrix> transpose()const=0;
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

    //赋值运算符
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
                std::cout<<std::setw(10)<<std::setprecision(4)<<data[i][j] << " ";
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
                json << std::fixed << std::setprecision(6) << data[i][j];
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
        return std::make_unique<Matrix>(*this);
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

    //内置的矩阵乘法,时间复杂度O(n^3)*(*)
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

    //内置的转置运算
    std::unique_ptr<AbstractMatrix> transpose()const override{
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
        json << "  \"sparsity\": " << std::fixed << std::setprecision(4) 
             << (1.0 - static_cast<double>(nonZeroCount) / (rows * cols)) << ",\n";
        json << "  \"elements\": [\n";
        
        for (size_t i = 0; i < elements.size(); ++i) {
            const auto& elem = elements[i];
            json << "    {\"row\": " << elem.row 
                 << ", \"col\": " << elem.col 
                 << ", \"value\": " << std::fixed << std::setprecision(6) << elem.value << "}";
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

    //内置的标量乘法
    std::unique_ptr<AbstractMatrix> scalarMultiply(double scalar) const override{
        auto result = std::make_unique<SparseMatrix>(rows, cols);
        if (scalar == 0.0)//0直接返回
            return result;
        for (const auto& elem : elements)
            result->set(elem.row, elem.col, elem.value*scalar);
        return result;
    }

    //内置的矩阵转置
    std::unique_ptr<AbstractMatrix> transpose() const override{
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
    AbstractTriangularMatrix(int size) : AbstractMatrix(size, size){
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

    //转置，继续往下抛
    std::unique_ptr<AbstractMatrix> transpose() const override=0;

public:     //三角矩阵实例对象特有的方法
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
    UpperTriangularMatrix(int size) : AbstractTriangularMatrix(size) {}
    UpperTriangularMatrix(const Matrix& dense) : AbstractTriangularMatrix(dense.getRows()){
        if(!dense.isSquare())
            throw std::invalid_argument("Triangular matrix must be square");
        for(int i = 0; i < rows; ++i){
            for(int j = 0; j < cols; ++j){
                double v = dense.get(i, j);
                if(!isInStoredRegion(i, j) && v != 0.0)
                    throw std::invalid_argument("Matrix does not match triangular type for UpperTriangularMatrix");
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
        json << "  \"determinant\": " << std::fixed << std::setprecision(6) << determinant() << ",\n";
        json << "  \"data\": [\n";
        
        for (int i = 0; i < rows; ++i) {
            json << "    [";
            for (int j = 0; j < cols; ++j) {
                json << std::fixed << std::setprecision(6) << get(i, j);
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
    std::unique_ptr<AbstractMatrix> transpose() const override {
        auto result = std::make_unique<LowerTriangularMatrix>(rows);
        for (int i = 0; i < rows; ++i)
            for (int j = i; j < cols; ++j)  //只遍历上三角
                result->set(j, i, get(i, j)); //行列互换
        return result;
    }

public:    // 上三角矩阵特有的方法
    //求解Ax=b
    std::unique_ptr<AbstractMatrix> solve(const AbstractMatrix& b) const {
        // 回代法求解上三角系统 Ux = b
        if (b.getRows() != rows) {
            throw std::invalid_argument("Dimension mismatch");
        }
        
        auto x = std::make_unique<Matrix>(rows, b.getCols());
        
        for (int k = 0; k < b.getCols(); ++k) {
            for (int i = rows - 1; i >= 0; --i) {
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
    LowerTriangularMatrix(int size) : AbstractTriangularMatrix(size) {}
    LowerTriangularMatrix(const Matrix& dense) : AbstractTriangularMatrix(dense.getRows()) {
        if(!dense.isSquare())
            throw std::invalid_argument("Triangular matrix must be square");
        for(int i = 0; i < rows; ++i){
            for(int j = 0; j < cols; ++j){
                double v = dense.get(i, j);
                if(!isInStoredRegion(i, j) && v != 0.0)
                    throw std::invalid_argument("Matrix does not match triangular type for LowerTriangularMatrix");
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
        json << "  \"determinant\": " << std::fixed << std::setprecision(6) << determinant() << ",\n";
        json << "  \"data\": [\n";
        
        for (int i = 0; i < rows; ++i) {
            json << "    [";
            for (int j = 0; j < cols; ++j) {
                json << std::fixed << std::setprecision(6) << get(i, j);
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
    std::unique_ptr<AbstractMatrix> transpose() const override {
        auto result = std::make_unique<UpperTriangularMatrix>(rows);
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j <= i; ++j)        //只遍历下三角
                result->set(j, i, get(i, j));   //行列互换
        return result;
    }

public:    //下三角矩阵特有的方法
    //下三角矩阵与下三角矩阵相乘仍然是下三角矩阵
    std::unique_ptr<AbstractMatrix> multiplyLower(const LowerTriangularMatrix& other) const {
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



}//namespace MatCal::Utils

namespace MatCal::Algorithm{
    using Utils::AbstractMatrix;
    using Utils::Matrix;
    using Utils::SparseMatrix;

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
    std::unique_ptr<AbstractMatrix> swapRows(const AbstractMatrix& A, int r1, int r2){
        if (r1 == r2)   //相同行，返回副本
            return A.toNormalMatrix();
        checkRowBounds(A.getRows(), r1);
        checkRowBounds(A.getRows(), r2);
        //转换为稠密矩阵副本进行操作，确保操作的通用性
        auto result = dynamic_cast<Matrix*>(A.toNormalMatrix().release());
        std::unique_ptr<AbstractMatrix> safe_result(result);
        //执行行交换
        std::swap((*result)[r1], (*result)[r2]);
        return safe_result;
    }

    //行初等变换2:行乘法 R_r <- scalar * R_r
    std::unique_ptr<AbstractMatrix> scaleRow(const AbstractMatrix& A, int r, double scalar) {
        if (scalar == 1.0)      //不变，返回副本
            return A.toNormalMatrix();
        checkRowBounds(A.getRows(), r);
        //转换为稠密矩阵副本
        auto result = dynamic_cast<Matrix*>(A.toNormalMatrix().release()); 
        std::unique_ptr<AbstractMatrix> safe_result(result);
        //执行行乘法
        for(int j = 0; j < A.getCols(); ++j)
            (*result)[r][j] *= scalar;
        return safe_result;
    }

    //行初等变换2:行加法 R_r_target <- R_r_target + scalar * R_r_source
    std::unique_ptr<AbstractMatrix> addScaledRow(const AbstractMatrix& A, int r_target, int r_source, double scalar) {
        if (scalar == 0.0 || r_target == r_source)
            return A.toNormalMatrix();  //不变，返回副本
        checkRowBounds(A.getRows(), r_target);
        checkRowBounds(A.getRows(), r_source);
        //转换为稠密矩阵副本
        auto result = dynamic_cast<Matrix*>(A.toNormalMatrix().release()); 
        std::unique_ptr<AbstractMatrix> safe_result(result);
        //执行加法
        for (int j = 0; j < A.getCols(); ++j)
            (*result)[r_target][j] += scalar * (*result)[r_source][j];
        return safe_result;
    }


    //初等列变换1: 列交换 C_c1 <-> C_c2
    std::unique_ptr<AbstractMatrix> swapCols(const AbstractMatrix& A, int c1, int c2) {
        if(c1 == c2)
            return A.toNormalMatrix();
        checkColBounds(A.getCols(), c1);
        checkColBounds(A.getCols(), c2);

        auto result = dynamic_cast<Matrix*>(A.toNormalMatrix().release()); 
        std::unique_ptr<AbstractMatrix> safe_result(result);
        
        int rows = A.getRows();
        //遍历每一行，交换目标列的元素
        for (int i = 0; i < rows; ++i)
            std::swap((*result)[i][c1], (*result)[i][c2]);
        return safe_result;
    }

    //初等列变换 2: 列乘法 C_c <- scalar * C_c
    std::unique_ptr<AbstractMatrix> scaleCol(const AbstractMatrix& A, int c, double scalar) {
        if (scalar == 1.0) {
            return A.toNormalMatrix();
        }
        checkColBounds(A.getCols(), c);

        auto result = dynamic_cast<Matrix*>(A.toNormalMatrix().release()); 
        std::unique_ptr<AbstractMatrix> safe_result(result);

        int rows = A.getRows();
        //遍历每一行，缩放目标列的元素
        for (int i = 0; i < rows; ++i)
            (*result)[i][c] *= scalar;
        return safe_result;
    }

    //初等列变换 3: 列加法 C_c_target <- C_c_target + scalar * C_c_source
    std::unique_ptr<AbstractMatrix> addScaledCol(const AbstractMatrix& A, int c_target, int c_source, double scalar) {
        if (scalar == 0.0 || c_target == c_source)
            return A.toNormalMatrix();
        checkColBounds(A.getCols(), c_target);
        checkColBounds(A.getCols(), c_source);

        auto result = dynamic_cast<Matrix*>(A.toNormalMatrix().release()); 
        std::unique_ptr<AbstractMatrix> safe_result(result);

        int rows = A.getRows();
        //遍历每一行，执行列加法
        for (int i = 0; i < rows; ++i)
            (*result)[i][c_target] += scalar * (*result)[i][c_source];
        return safe_result;
    }


}//namespace MatCal::Algorithm


#endif//MATRIX_HPP