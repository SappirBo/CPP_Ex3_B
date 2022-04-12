#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>

using namespace std;

namespace zich
{

    class Matrix{
        private:
        public:

        double** mat;
        int colums;
        int rows;

        Matrix();
        Matrix(vector<double> values,int Rows, int Colums);
        Matrix(int Row,int Colums);
        ~Matrix();

        // Getters / Setters
        int getCols(){
            return this->colums;
        }
        int getRows(){
            return this->rows;
        }

        // Helpers 
        void printMat();

        // Arithmetic actions on single Matrix.
        Matrix& operator++();
        Matrix operator++(int);
        Matrix& operator--();
        Matrix operator--(int);

        // Arithmetic actions between two Matrix.
        friend Matrix operator+(const Matrix &matrix1,const Matrix &matrix2);
        void operator+=(const Matrix &matrix2) const;
        friend Matrix operator-(const Matrix &matrix1, const Matrix &matrix2);
        void operator-=(const Matrix &matrix2) const;
        friend Matrix operator*(const Matrix &matrix1, const Matrix &matrix2);
        void operator*=(const Matrix &matrix2);

        Matrix& operator=(const Matrix& matrix);
        
        // Arithmetic actions between Matrix and scalar.
        friend Matrix operator*(const double s , const Matrix &matrix1);
        friend Matrix operator*(const Matrix &matrix1, const double s );
        void operator*=(double s) const;
        
        // Compare Operator between two matrix.
        friend bool operator<(const Matrix &matrix1, const Matrix &matrix2);
        friend bool operator<=(const Matrix &matrix1, const Matrix &matrix2);
        friend bool operator>(const Matrix &matrix1, const Matrix &matrix2);
        friend bool operator>=(const Matrix &matrix1, const Matrix &matrix2);
        friend bool operator==(const Matrix &matrix1, const Matrix &matrix2);
        friend bool operator!=(const Matrix &matrix1, const Matrix &matrix2);


        // Onitery Actions
        Matrix operator-() const;
        Matrix operator+() const;

        
        // Ostream/Istream
        friend std::ostream& operator<<(std::ostream &o,const Matrix &matrix);
        friend std::istream & operator>>(std::istream &o, Matrix &matrix);

    };
    
} 
