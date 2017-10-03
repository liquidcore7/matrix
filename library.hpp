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

    template <typename Fnc>
    // functor should be of the following signature: void functor(T&)
    // example: [] (T& a) {a += 15}
    void apply(const Fnc&);

    Matrix& operator+=(const Matrix&);
    Matrix& operator*=(const T&);

    Matrix operator+(const Matrix&, const Matrix&);
    Matrix operator*(const Matrix&, const T&);

};


/* ********************************************************************************
 * IMPLEMENTATION
 * ***************************************************************************** */

template <size_t N, size_t M, typename T>
Matrix<N, M, T>::Matrix(std::istream &is)
{
    for (size_t i = 0; i < N; ++i)
        for (size_t j = 0; j < M; ++j)
            is >> _grid[i][j];
}

template <size_t N, size_t M, typename T>
Matrix<N, M, T>::Matrix(const std::string &s)
{
    std::istringstream iss(s);
    Matrix<N, M, T>::Matrix(iss);
}

template <size_t N, size_t M, typename T>
Matrix<N, M, T>::Matrix(const T &init)
{
    for(auto &row : _grid)
        for (auto &cell : row)
            cell = init;
}

template <size_t N, size_t M, typename T>
template<typename Fnc>
void Matrix<N, M, T>::apply(const Fnc &functor)
{
    for (auto& row : _grid)
        for (auto& cell : row)
            functor(cell);
}

template <size_t N, size_t M, typename T>
Matrix<N, M, T> &Matrix<N, M, T>::operator+=(const Matrix<N, M, T> &rhs)
{
    for (size_t i = 0; i < N; ++i)
        for (size_t j = 0; j < M; ++j)
            _grid[i][j] += rhs._grid[i][j];
    return *this;
}

template <size_t N, size_t M, typename T>
Matrix<N, M, T> &Matrix<N, M, T>::operator*=(const T &scalar)
{
    this->apply([&scalar] (T& elem) {elem *= scalar;});
    return *this;
}

template <size_t N, size_t M, typename T>
Matrix<N, M, T> operator+(const Matrix<N, M, T>& lhs, const Matrix<N, M, T>& rhs)
{
    Matrix<N, M, T> local;
    for (size_t i = 0; i < N; ++i)
        for (size_t j = 0; j < M; ++j)
            local._grid[i][j] = lhs._grid[i][j] + rhs._grid[i][j];
    return local;
};


template <size_t N, size_t M, typename T>
Matrix<N, M, T> operator*(const Matrix<N, M, T>&lhs, const T& scalar)
{
    Matrix<N, M, T> local;
    for (size_t i = 0; i < N; ++i)
        for (size_t j = 0; j < M; ++j)
            local._grid[i][j] = lhs._grid[i][j] * scalar;
    return local;
};


#endif