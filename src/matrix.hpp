/**
 * @file matrix.hpp
 * @author Oct10th
 * @brief The header of this simple mathematical matrix implementation.
 * @details A matrix implementation includes addition, multiplication, inverse, etc.
 * @mainpage Simple Matrix Implementation
*/

#ifndef OTEN_MATRIX_H__
#define OTEN_MATRIX_H__

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <format>

namespace OTen {

class MatrixOperationException: public std::logic_error {
public:
    explicit MatrixOperationException(const std::string &arg): logic_error(arg) {}
    explicit MatrixOperationException(const char *arg): logic_error(arg) {}
    virtual ~MatrixOperationException() noexcept = default;
};

class Matrix {
private:
    unsigned int rowCnt;
    unsigned int colCnt;
    /// @brief the pointer to the data of the matrix class
    double **pData;
    void initialize(int row, int col);
    void freeMemory();
    double get_cofactor(int elemRow, int elemCol) const;
public:
    // Matrix get_cofactor(int elemRow, int elemCol) const;
    Matrix();
    Matrix(int rowCount, int colCount);
    Matrix(int rowCount, int colCount, std::initializer_list<double> list);

    Matrix(const Matrix&);
    Matrix& operator=(const Matrix&);
    Matrix(Matrix &&);
    Matrix& operator=(Matrix&&);

    virtual ~Matrix();

    void print_on_ostream(std::ostream &os) const;
    friend std::ostream& operator<<(std::ostream &os, const Matrix &mat);
    double get_elem(int row, int col) const;
    void set_elem(int row, int col, double value);
    void set_by_init_list(std::initializer_list<double> list);
    Matrix& operator=(std::initializer_list<double> list);

    void add(const Matrix &rhs);
    Matrix operator+(const Matrix &rhs) const;
    void subtract(const Matrix &rhs);
    Matrix operator-(const Matrix &rhs) const;
    void scalar_multiply_by(double num);
    friend Matrix operator*(const double num, const Matrix &mat);
    Matrix operator*(const Matrix &rhs) const;

    double det() const;
    Matrix adj() const;
    Matrix inv() const;
};

} // namespace OTen

#endif