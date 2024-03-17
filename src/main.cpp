#include "matrix.hpp"
#include <iostream>

int main(void)
{
    OTen::Matrix mat1(2, 3, {2, 3, -4, 7, 8, 9});
    OTen::Matrix mat2(3, 3, {8, 9, 13, 1, 4, 7, -1, 3, 5});
    OTen::Matrix mat3 = mat1 * mat2;
    std::cout << mat3 << std::endl;
}