/**
 * @file matrix.cpp
 * @author Oct10th
 * @brief The source code of this simple mathematical matrix implementation.
 * @details A matrix implementation includes addition, subtraction, multiplication, inverse, etc.
 * @mainpage Simple Matrix Implementation
*/

#include "matrix.hpp"


namespace OTen {


/**
 * @brief A private function to allocate and initialize memory.
 * @param row the row count of the matrix that will be created
 * @param col the column count of the matrix that will be created
*/
void Matrix::initialize(int row, int col)
{
    if (row > 0 && col > 0){
        pData = new double* [row];
        for (int i = 0; i < row; i++) {
            pData[i] = new double[col];
        }
    }
    else {
        rowCnt = 0;
        colCnt = 0;
        pData = nullptr;
    }
}

/**
 * @brief A private function to free the memory of the Matrix class.
*/
void Matrix::freeMemory()
{
    if (pData != nullptr){
        for (int i = 0; i < rowCnt; i++) {
            delete [] pData[i];
        }
        delete [] pData;
    }
}

/**
 * @brief A private function to give the cofactor of the (excludeRow, excludeCol) entry of the Matrix class.
 * @param excludeRow the row that the cofactor would exclude
 * @param excludeCol the column that the cofactor would exclude
*/
double Matrix::get_cofactor(const int excludeRow, const int excludeCol) const
{
    Matrix cofactor(rowCnt - 1, colCnt - 1);

    for (int i = 0, curr = 0; i < rowCnt; i++) {
        for (int j = 0; j < colCnt; j++) {
            if (i != excludeRow - 1 && j != excludeCol - 1) {
                cofactor.pData[curr/cofactor.colCnt][curr%cofactor.colCnt] = pData[i][j];
                ++curr;
            }
            else {
                continue;
            }
        }
    }

    return ((excludeRow + excludeCol) % 2 == 0) ? cofactor.det() : -cofactor.det();
}

/**
 * @brief The default constructor of the Matrix class.
 * It gives a Matrix class with no data, and the row count and the column count is both 0.
*/
Matrix::Matrix(): rowCnt(0), colCnt(0)
{
    pData = nullptr;
}

/**
 * @brief The constructor with parameters so that you can get any size of the Matrix class.
 * Each element will be initialized to 0.
 * @param rowCount the row count of the matrix
 * @param colCount the column count of the matrix
*/
Matrix::Matrix(int rowCount, int colCount): rowCnt(rowCount), colCnt(colCount)
{
    initialize(rowCount, colCount);
}

/**
 * @brief The constructor with row count, column count and initializer_list 
 * so that you can get any size of the Matrix class.
 * Additionally, you can set the initial value of the elements in the Matrix class.
 * The order of the initialization is 
 * pData[0][0], [0][1], [0][2], ..., [1][0], [1][1], [1][2], ..., [2][0], [2][1], [2][2], ...
 * @attention If the initializer_list's size is bigger than the total element count of the matrix, 
 * it will throw an std::out_of_range exception.
 * @attention If the initializer_list's size is smaller, 
 * then all the left elements will be initialized to 0.
*/
Matrix::Matrix(int rowCount, int colCount, std::initializer_list<double> list)
{
    rowCnt = rowCount;
    colCnt = colCount;
    initialize(rowCount, colCount);
    if (list.size() <= rowCnt * colCnt) {
        int i = 0;
        for (auto &num: list) {
            pData[i/colCnt][i%colCnt] = num;
            ++i;
        }
    }
    else {
        throw std::out_of_range(
            "OTen::Matrix::Matrix(int,int,initializer_list<double>):\n\tthe size of the init_list is bigger than the total element count of the matrix."
            );
    }
}

/**
 * @brief The move constructor of the Matrix class.
*/
Matrix::Matrix(Matrix &&rhs)
{
    rowCnt = rhs.rowCnt;
    colCnt = rhs.colCnt;
    pData = rhs.pData;

    rhs.rowCnt = rhs.colCnt = 0;
    rhs.pData = nullptr;
}

/**
 * @brief The copy constructor of the Matrix class.
*/
Matrix::Matrix(const Matrix &rhs)
{
    rowCnt = rhs.rowCnt;
    colCnt = rhs.colCnt;
    initialize(rhs.rowCnt, rhs.colCnt);

    for (int i = 0; i < rowCnt; i++) {
        for (int j = 0; j < colCnt; j++) {
            pData[i][j] = rhs.pData[i][j];
        }
    }
}

/**
 * @brief The copy assignment operator of the Matrix class.
*/
Matrix& OTen::Matrix::operator=(const Matrix &rhs)
{
    rowCnt = rhs.rowCnt;
    colCnt = rhs.colCnt;
    initialize(rhs.rowCnt, rhs.colCnt);

    for (int i = 0; i < rowCnt; i++) {
        for (int j = 0; j < colCnt; j++) {
            pData[i][j] = rhs.pData[i][j];
        }
    }
    return *this;
}

/**
 * @brief The move assignment operator of the Matrix class.
*/
Matrix& Matrix::operator=(Matrix &&rhs)
{
    rowCnt = rhs.rowCnt;
    colCnt = rhs.colCnt;
    pData = rhs.pData;

    rhs.rowCnt = rhs.colCnt = 0;
    rhs.pData = nullptr;

    return *this;
}

/**
 * @brief The destructor of the Matrix class.
*/
Matrix::~Matrix()
{
    freeMemory();
}

/**
 * @brief Print the current Matrix class to a ostream.
 * @param os the ostream that will be printed to.
*/
void Matrix::print_on_ostream(std::ostream &os) const
{
    std::ios::fmtflags fmt(os.flags());

    os << "[ " << std::right;
    for (int i = 0; i < rowCnt; ++i) {

        for (int j = 0; j < colCnt; ++j) {
            if (i != 0 && j == 0) os << "  ";
            os << std::format("{:^7}", pData[i][j]) << " ";
        }

        if(i != rowCnt - 1) os << "\n";
        else os << "]" << std::flush;
        
        os.flags(fmt);
    }
}

/**
 * @brief The output operator.
*/
std::ostream& operator<<(std::ostream &os, const Matrix &mat)
{
    mat.print_on_ostream(os);
    return os;
}

/**
 * @brief The function can give you the number of the element which your coordinate locates.
 * @param row the row where the element locates
 * @param col the column where the element locates
 * @attention If the coordinate you give is out of the range of the matrix, 
 * it will throw a std::out_of_range exception.
*/
double Matrix::get_elem(int row, int col) const
{
    if (row > 0 && row <= rowCnt && col > 0 && col <= colCnt) {
        return pData[row-1][col-1];
    }
    else {
        throw std::out_of_range("OTen::Matrix::get_elem() out of range");
    }
}

/**
 * @brief The function can set a single element you want to a value you give.
 * @param row the row where the element locates
 * @param col the column where the element locates
 * @param value the value to set to the element
 * @attention If the coordinate you give is out of the range of the matrix, 
 * it will throw a std::out_of_range exception.
*/
void Matrix::set_elem(int row, int col, double value)
{
    if (row > 0 && row <= rowCnt && col > 0 && col <= colCnt) {
        pData[row-1][col-1] = value;
    }
    else {
        throw std::out_of_range("OTen::Matrix::set_elem() out of range");
    }
}

/**
 * @brief Set the value of the elements in the Matrix class according to the initializer_list.
 * The order of the assignment is 
 * pData[0][0], [0][1], [0][2], ..., [1][0], [1][1], [1][2], ..., [2][0], [2][1], [2][2], ...
 * @return The matrix itself which has been changed.
 * @attention If the initializer_list's size is bigger than the total element count of the matrix, 
 * it will throw a std::out_of_range exception.
 * @attention If the initializer_list's size is smaller, then all the left elements will be set to 0.
*/
Matrix& Matrix::operator=(std::initializer_list<double> list)
{
    this->set_by_init_list(list);
    return *this;
}

/**
 * @brief Set the value of the elements in the Matrix class according to the initializer_list.
 * The order of the assignment is 
 * pData[0][0], [0][1], [0][2], ..., [1][0], [1][1], [1][2], ..., [2][0], [2][1], [2][2], ...
 * @attention If the initializer_list's size is bigger than the total element count of the matrix,
 * it will throw a std::out_of_range exception.
 * @attention If the initializer_list's size is smaller, then all the left elements will be set to 0.
*/
void Matrix::set_by_init_list(std::initializer_list<double> list)
{
    if (list.size() <= rowCnt * colCnt) {
        int i = 0;
        for (auto &num: list) {
            pData[i/colCnt][i%colCnt] = num;
            ++i;
        }
        while (i < rowCnt * colCnt) {
            pData[i/colCnt][i%colCnt] = 0;
            ++i;
        }
    }
    else {
        throw std::out_of_range(
            "OTen::Matrix::set(initializer_list<double>):\n\tthe size of the init_list is bigger than the total element count of the matrix."
            );
    }

}

/**
 * @brief The matrix addition algorithm.
 * @attention The matrix addition requires the two Matrix class to be of the same dimension. 
 * If not, the function will throw a MatrixOperationException.
*/
void Matrix::add(const Matrix &rhs)
{
    if (rowCnt == rhs.rowCnt && colCnt == rhs.colCnt) {

        for (int i = 0; i < rowCnt; i++) {
            for (int j = 0; j < colCnt; j++) {
                pData[i][j] += rhs.pData[i][j];
            }
        }

    }
    else {
        throw MatrixOperationException(
            "OTen::Matrix::add(): addition of two Matrix class of different dimensions"
            );
    }
}

/**
 * @brief The addition operator with the matrix addition algorithm.
 * @attention The matrix addition requires the two Matrix class to be of the same dimension. 
 * If not, the function will throw a MatrixOperationException.
*/
Matrix Matrix::operator+(const Matrix &rhs) const
{
        if (rowCnt == rhs.rowCnt && colCnt == rhs.colCnt) {

        Matrix result(rowCnt, colCnt);

        for (int i = 0; i < rowCnt; i++) {
            for (int j = 0; j < colCnt; j++) {
                result.pData[i][j] = pData[i][j] + rhs.pData[i][j];
            }
        }

        return result;
    }
    else {
        throw MatrixOperationException(
            "OTen::Matrix::operator+(): addition of two Matrix class of different dimensions"
            );
    }

}

/**
 * @brief The matrix subtraction algorithm.
 * @attention The matrix subtraction requires the two Matrix class to be of the same dimension. 
 * If not, the function will throw a MatrixOperationException.
*/
void Matrix::subtract(const Matrix &rhs)
{
    if (rowCnt == rhs.rowCnt && colCnt == rhs.colCnt) {

        for (int i = 0; i < rowCnt; i++) {
            for (int j = 0; j < colCnt; j++) {
                pData[i][j] -= rhs.pData[i][j];
            }
        }

    }
    else {
        throw MatrixOperationException(
            "OTen::Matrix::subtract(): subtraction of two Matrix class of different dimensions"
            );
    }
}

/**
 * @brief The subtraction operator with the matrix subtraction algorithm.
 * @attention The matrix subtraction requires the two Matrix class to be of the same dimension. 
 * If not, the function will throw a MatrixOperationException.
*/
Matrix Matrix::operator-(const Matrix &rhs) const
{
        if (rowCnt == rhs.rowCnt && colCnt == rhs.colCnt) {

        Matrix result(rowCnt, colCnt);

        for (int i = 0; i < rowCnt; i++) {
            for (int j = 0; j < colCnt; j++) {
                result.pData[i][j] = pData[i][j] - rhs.pData[i][j];
            }
        }

        return result;
    }
    else {
        throw MatrixOperationException(
            "OTen::Matrix::operator-(): subtraction of two Matrix class of different dimensions"
            );
    }

}

/**
 * @brief The matrix scalar multiplication algorithm.
 * @param num the number which will multiply the current Matrix class.
 * @attention This function has no return value and it only changes the current Matrix class.
*/
void Matrix::scalar_multiply_by(const double num)
{
    for (int i = 0; i < rowCnt; i++) {
        for (int j = 0; j < colCnt; j++) {
            pData[i][j] *= num;
        }
    }
}

/**
 * @brief The multiplication operator with the matrix scalar multiplication algorithm.
 * @param num the scalar which will multiply the Matrix class mat
 * @param mat the Matrix class which will be multiplied by the scalar num
 * @attention This function returns the result and it doesn't change the Matrix class in the parameter.
*/
Matrix operator*(const double num, const Matrix &mat)
{
    Matrix result(mat);
    result.scalar_multiply_by(num);
    return result;
}

/**
 * @brief The multiplication operator with matrix multiplication algorithm.
 * @param rhs the Matrix class that will be multiplied by the current Matrix class.
 * @return the result Matrix class of the matrix multiplication.
 * @attention This operator doesn't change any data and it only returns the result.
*/
Matrix Matrix::operator*(const Matrix &rhs) const
{
    if (colCnt == rhs.rowCnt) {

        Matrix result(rowCnt, rhs.colCnt);

        for (int i = 0; i < result.rowCnt; ++i) {
            for (int j = 0; j < result.colCnt; ++j) {
                for (int k = 0; k < colCnt; ++k) {
                    result.pData[i][j] += pData[i][k] * rhs.pData[k][j];
                }
            }
        }

        return result;
    }
    else {
        throw MatrixOperationException(
            "OTen::Matrix::matrix_multiply_by():\n\tmultiplication of a m*n matirx and a p*q matrix but n != p"
            );
    }
}

/**
 * @brief The function gives the determinant of the current Matrix class.
 * @attention The determinant only exists in square matrix. 
 * If the current Matrix class is a nonsquare matrix,
 * this function will throw a MatrixOperationException.
*/
double Matrix::det() const
{
    if (rowCnt == colCnt) {

        if (rowCnt == 2) {
            // det = ad - bc
            return pData[0][0] * pData[1][1] - pData[0][1] * pData[1][0];
        }


        double result = 0;

        for (int i = 1; i <= rowCnt; ++i) {
            result += pData[i-1][0] * get_cofactor(i, 1);
        }

        return result;
    }
    else {
        throw MatrixOperationException(
            "OTen::Matrix::det(): evaluation of the determinant of nonsquare matrix"
            );
    }
}

} // namespace OTen
