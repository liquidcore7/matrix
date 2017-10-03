#ifndef MATRIX_LIBRARY_H
#define MATRIX_LIBRARY_H

#include <array>
#include <sstream>

template <size_t N, size_t M, typename T>
//  NxM matrix of elements of type T
class Matrix
{
    std::array< std::array<T, M>, N > _grid;

public:
    Matrix() = default;
    explicit Matrix(std::istream&);
    explicit Matrix(const std::string&);
    explicit Matrix(const T&);

    Matrix(const Matrix&) = default;
    Matrix(Matrix&& r) noexcept : _grid(std::move(r._grid))
    {}

    Matrix& operator=(const Matrix&) = default;
    Matrix& operator=(Matrix&& r) noexcept
    {
        this->_grid = std::move(r._grid);
    }



};

template <size_t N, size_t M, typename T>
Matrix::Matrix(std::istream &is)
{
    for (size_t i = 0; i < N; ++i)
        for (size_t j = 0; j < M; ++j)
            is >> _grid[i][j];
}

Matrix::Matrix(const std::string &s)
{
    std::istringstream iss(s);
    Matrix::Matrix(iss);
}

template <size_t N, size_t M, typename T>
Matrix::Matrix(const T &init)
{
    for (size_t i = 0; i < N; ++i)
        for (size_t j = 0; j < M; ++j)
            _grid[i][j] = init;
}


#endif