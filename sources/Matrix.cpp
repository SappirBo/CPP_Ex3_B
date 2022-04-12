#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include "Matrix.hpp"

using namespace std;
using namespace zich;

// Helpers -> exeption checking and multiplication.
void dim_check_Addition( const Matrix &m1, const Matrix &m2){
    if(m1.rows != m2.rows || m1.colums !=m2.colums){
        throw invalid_argument("Eror: Matrixes Dimentions Must be Equals!.");
    }
}
void dim_check_Mult(const Matrix &m1,const  Matrix &m2){
    if(m1.colums != m2.rows){
        throw invalid_argument("Eror: The first Matrix Column Dimention Must be Equal to the second Matrix Row Dimention!");
    }
}
double row_col_mult(const Matrix& m1,const Matrix& m2, size_t row,size_t col){
    double ans = 0;
    for(int i=0; i<m1.colums; i++){
        ans += m1.mat[row][i] * m2.mat[i][col];
    }
    return ans;
}
double sum_of_matrix(const Matrix &m){
    double ans =0;
    for(size_t i=0; i<m.rows; i++){
        for(size_t j=0; j<m.colums; j++){
            ans += m.mat[i][j];
        }
    }
    return ans;
}
void check_valid_compare(const Matrix &m1, const Matrix &m2){
    if(m1.rows != m2.rows || m1.colums !=m2.colums){
        throw invalid_argument("Matrixes Dimentions Must be Equals!.");
    }
}
void check_valid_string(string str){
    int counter1 =0;
    int counter2 =0;
    bool flag = true;
    for(size_t i=0; i<str.size();i++){
        if(str[i] == ']'){
            counter1 ++;
            flag = true;
        }
        
        if(str[i] == ','){
            counter2 ++;
            if(!flag || str[i+1] != ' '){
                throw invalid_argument("Eror: Bad Input, code 1.");
            }
            flag = false;
        }
    }

    if(counter2 + 1 != counter1){
        throw invalid_argument("Eror: Bad Input, code 2.");
    }
}
void check_for_constructor(const Matrix &m1, const Matrix &m2){
    if(m1.colums < 0 || m2.colums<0 || m1.rows<0 || m2.rows<0){
        throw invalid_argument("Matrixes Dimentions Must be grater then 0 !.");
    }

}


// Constractors -> Empty Matrix, Identity Metrix, Given values Matrix. 
Matrix::Matrix(){
    this->colums = 0;
    this->rows = 0;
    this->mat = NULL;
}
Matrix::Matrix(int Row,int Colums){
    if(Row<=0 || Colums<=0){
         throw invalid_argument("Matrixes Dimentions Must be >= 0 !.");
    }
    this->rows = Row;
    this->colums = Colums;
    this->mat = new double*[(size_t)Row];
    for(int i=0; i<Row;i++){
        this->mat[i] = new double[(size_t)Colums]; 
    }
    size_t x=0;
    for(int i=0; i<Row; i++){
        for(int j=0; j<Colums; j++){
            if(i == j){
                this->mat[i][j] = 1;
            }
            this->mat[i][j] = 0;
        }
    }
}
Matrix::Matrix(std::vector<double> values, int row, int col){
    if(row*col != values.size()){
         throw invalid_argument("Eror: bad values!.");
    }
    this->rows = row;
    this->colums = col;
    this->mat = new double*[(size_t)row];
    for(size_t i=0; i<row;i++){
        this->mat[i] = new double[(size_t)col]; 
    }
    size_t x=0;
    for(size_t i=0; i<row; i++){
        for(size_t j=0; j<col; j++){
            this->mat[i][j] = values[x++];
        }
    }
}
Matrix::~Matrix(){
    for(size_t i=0; i< this->rows; i++){
        delete[] this->mat[i];
    }
    delete[] this->mat;
}

// operator =
Matrix& zich::Matrix::operator=(const Matrix& matrix){
    // Check a=a case.
    if(this == &matrix){
        return *this;
    }
    //Delete the old mat.
    Matrix old_mat = +(*this);
    for(size_t i=0; i<this->rows; i++){
        delete [] this->mat[i];
    }
    delete [] this->mat;
    size_t Rows = (size_t) matrix.rows;
    size_t Cols =  (size_t) matrix.colums;
    this->rows = Rows;
    this->colums = Cols;
    this->mat = new double*[Rows];
    for(size_t i=0;i<Rows; i++){
        this->mat[i] = new double[Cols];
    }
    for(size_t i=0;i<Rows; i++){
        for(size_t j=0; j<Cols; j++){
            this->mat[i][j] = matrix.mat[i][j];
        }
    }
    return *this;
}

// Arithmetic actions between two Matrix (+, +=, -, -=, *. *=).
Matrix zich::operator+(const Matrix &matrix1,const Matrix &matrix2){
    dim_check_Addition(matrix1, matrix2);
    std::vector<double> arr;
    arr.resize((unsigned int)(matrix1.colums * matrix1.rows));

    int row = matrix1.rows;
    int col = matrix1.colums;

    size_t curr =0;
    for(size_t i=0; i<row; i++){
        for(size_t j=0; j<col; j++){
            arr[curr++] = matrix1.mat[i][j] + matrix2.mat[i][j];
        }
    }
    return Matrix(arr,row,col);
}
void zich::Matrix::operator+=( const Matrix &matrix2) const{
    dim_check_Addition(*this, matrix2);
    for(size_t i=0; i<this->rows; i++){
        for(size_t j=0; j<matrix2.colums; j++){
            this->mat[i][j] += matrix2.mat[i][j];
        }
    }
}
Matrix zich::operator-(const Matrix &matrix1, const Matrix &matrix2){
    dim_check_Addition(matrix1, matrix2);
    Matrix tmp(matrix1.rows, matrix1.colums);
    int row = matrix1.rows;
    int col = matrix1.colums;
    for(int i=0; i<row; i++){
        for(int j=0; j<col; j++){
            tmp.mat[i][j] = matrix1.mat[i][j] - matrix2.mat[i][j];
        }
    }
    return tmp;

}
void zich::Matrix::operator-=(const Matrix &matrix2) const{
    dim_check_Addition(*this, matrix2);
    for(size_t i=0; i<this->rows; i++){
        for(size_t j=0; j<matrix2.colums; j++){
            this->mat[i][j] -= matrix2.mat[i][j];
        }
    }
}
Matrix zich::operator*(const Matrix& matrix1,const Matrix& matrix2){
    dim_check_Mult(matrix1,matrix2);
    vector<double> arr;
    arr.resize((unsigned int)(matrix1.rows * matrix2.colums));
    size_t count = 0;
    for(size_t i=0;i<matrix1.rows;i++){
        for(size_t j=0;j<matrix2.colums;j++){
            double num = 0;
            for(size_t k=0; k<matrix1.colums; k++){
                num += matrix1.mat[i][k] * matrix2.mat[k][j];
            }
            arr[count++] = num;
        }
    }
    return Matrix(arr,matrix1.rows,matrix2.colums);
}
void Matrix::operator*=(const Matrix& matrix2){
    dim_check_Mult(*this,matrix2);
    vector<double> arr;
    arr.resize((unsigned int)(this->rows * matrix2.colums));
    size_t index = 0;
    for(size_t i=0; i< this->rows; i++){
        for(size_t j=0; j<matrix2.colums; j++){
            double x=0;
            for(size_t k=0; k<this->colums; k++){
                x += this->mat[i][k] * matrix2.mat[k][j];
            }
            arr[index++] = x;
        }
    }
    index=0;
    for(size_t i=0; i<this->rows; i++){
        for(size_t j=0; j<this->colums; j++){
            this->mat[i][j] = arr[index++];
        }
    }
    this->colums = matrix2.colums;
}

// Increment/decrement operators(++ pre, ++ post, -- pre,--post).
Matrix & zich::Matrix::operator++(){
    for(size_t i=0; i<this->rows; i++){
        for(size_t j=0; j< this->colums; j++){
            this->mat[i][j] += 1;
        }
    }
    return *this;
}
Matrix Matrix::operator++(int){
    Matrix tmp = +(*this);
    ++(*this);
    return tmp;
}
Matrix & zich::Matrix::operator--(){
    for(size_t i=0;i<this->rows; i++){
        for(size_t j=0; j< this->colums; j++){
            this->mat[i][j]--;
        }
    }
    return *this;
}
Matrix zich::Matrix::operator--(int){
    Matrix tmp = +(*this);
    --(*this);
    return tmp;
}


// Arithmetic actions between Matrix and scalar (A*a, a*A, A*=a).
Matrix zich::operator*(const double s ,const Matrix &matrix1){
    std::vector<double> arr;
    arr.resize((unsigned int)(matrix1.colums * matrix1.rows));
    size_t curr = 0;
    int rows = matrix1.rows;
    int col = matrix1.colums;
    for(size_t i=0; i< rows; i++){
        for(size_t j=0; j<col; j++){
            arr[curr++] = matrix1.mat[i][j]*s; 
        }
    }
    return Matrix(arr,rows,col);
}
Matrix zich::operator*(const Matrix &matrix1,double s ){
    return s*matrix1;
}
void Matrix::operator*=(double s ) const{
    // Matrix tmp(this->rows,this->colums);
    for(size_t i=0; i< this->rows; i++){
        for(size_t j=0; j<this->colums; j++){
            this->mat[i][j] *= s; 
        }
    }
    // this->mat = tmp.mat;
}

// Unary Actions.
Matrix zich::Matrix::operator-()const{
    std::vector<double> arr;
    arr.resize((unsigned int)(this->colums * this->rows));
    int rows = this->rows;
    int col = this->colums;
    size_t curr = 0;
    for(size_t i=0; i< rows; i++){
        for(size_t j=0; j<col; j++){
            if(this->mat[i][j] != 0){
                arr[curr++]= this->mat[i][j] * -1; 
            }else{
                arr[curr++]=0;
            }
        }
    }
    return Matrix(arr,rows,col);
}
Matrix zich::Matrix::operator+()const{
    std::vector<double> arr;
    arr.resize((unsigned int)(this->colums * this->rows));
    int rows = this->rows;
    int col = this->colums;
    size_t curr = 0;
    for(size_t i=0; i< rows; i++){
        for(size_t j=0; j<col; j++){
            if(this->mat[i][j] != 0){
                arr[curr++]= this->mat[i][j] * 1; 
            }else{
                arr[curr++]=0;
            }
        }
    }
    return Matrix(arr,rows,col);
}

// Compare Operator between two matrix (<, <=, >, >=, ==, !=).
bool zich::operator<(const Matrix &matrix1, const Matrix &matrix2){
    check_valid_compare(matrix1,matrix2);
    bool ans = false;
    if(sum_of_matrix(matrix1)<sum_of_matrix(matrix2)){
        ans = true;
    }
    return ans;
}
bool zich::operator<=(const Matrix &matrix1, const Matrix &matrix2){
    check_valid_compare(matrix1,matrix2);
    bool ans = false;
    if(sum_of_matrix(matrix1)<=sum_of_matrix(matrix2)){
        ans = true;
    }
    return ans;
}
bool zich::operator>(const Matrix &matrix1, const Matrix &matrix2){
    check_valid_compare(matrix1,matrix2);
    bool ans = false;
    if(sum_of_matrix(matrix1)>sum_of_matrix(matrix2)){
        ans = true;
    }
    return ans;
}
bool zich::operator>=(const Matrix &matrix1, const Matrix &matrix2){
    check_valid_compare(matrix1,matrix2);
    bool ans = false;
    if(sum_of_matrix(matrix1)>=sum_of_matrix(matrix2)){
        ans = true;
    }
    return ans;
}
bool zich::operator==(const Matrix &matrix1, const Matrix &matrix2){
    check_valid_compare(matrix1,matrix2);
    bool ans = true;
    for(size_t i=0; i<matrix1.rows ; i++ ){
        for(size_t j=0; j<matrix1.colums; j++){
            if(matrix1.mat[i][j] != matrix2.mat[i][j]){
                ans = false;
                break;
            }
        }
    }
    return ans;
}
bool zich::operator!=(const Matrix &matrix1, const Matrix &matrix2){
    check_valid_compare(matrix1,matrix2);
    bool ans = false;
    if(sum_of_matrix(matrix1)!=sum_of_matrix(matrix2)){
        ans = true;
    }
    return ans;
}

// Ostream/Istream Actions.
std::ostream& zich::operator<<(std::ostream &o, const Matrix &matrix){
    int row = matrix.rows;
    int col = matrix.colums;
    for(size_t i=0; i<row; i++){
        o << "["; 
        for(size_t j=0; j<col; j++){
            if(j==col-1){
                o << matrix.mat[i][j]; 
            }else{
                o << matrix.mat[i][j] << " "; 
            }   
        }
        o << "]";
        if(i<row-1){
            o << "\n";
        } 
    }
    return o;
}
istream & zich::operator>>(std::istream &o, Matrix &matrix){
    // input given as '[a b c], [d e f], [g h i]', Where that means 3X3 Matrix
    string tmp_str;
    string str;
    while (getline(o,tmp_str)){
        str += tmp_str;   
    }
    check_valid_string(str);
    int col = 0;
    int row = 0;
    for(size_t i=0; i<str.size(); i++){
        if(str[i]==']'){
            row++;
        }
    }
    size_t spaces=0;
    size_t index = 0;
    tmp_str = "";
    for(size_t i=0; i<str.size();i++){
        if(str[i] != '[' && str[i] != ']'){
            tmp_str.push_back(str[i]);
        }
    }
    char ch =0 ;
    // unsigned int tmp = 0;
    while(ch != ','){
        ch=tmp_str[index++];
        if(ch==' '){
            spaces++;
        }
    }
    // tmp = spaces + 1;
    str = "";
    for(size_t i=0; i<tmp_str.size(); i++){
        if(tmp_str[i] != ',')
        {
            str.push_back(tmp_str[i]);
        }
    }
    string delimiter = " ";
    vector<double> arr;
    index = 0;
    size_t size =1;
    arr.resize(size++);
    while( arr.size() > 1) {
        string token = str.substr(0,str.find(delimiter));
        arr[index++]=stod(token);
        arr.resize(size++);
        str.erase(0, str.find(delimiter) + delimiter.length());
    }

    arr[index++]=stod(str);


    matrix.mat = new double*[(size_t)row];

    
    for(size_t i=0; i<row; i++){
        matrix.mat[i] = new double[(size_t)col];
    }
    matrix.colums = col;
    matrix.rows = row;
    index=0;

    for(size_t i = 0; i<col; i++){
        for(size_t j = 0; j<row; j++){
            matrix.mat[i][j]=arr[index++];
        }
    }

    return o;
}