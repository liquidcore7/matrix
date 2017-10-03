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



};

template <size_t N, size_t M, typename T>
Matrix<N, M, T> operator+(const Matrix<N, M, T>&, const Matrix<N, M, T>&);


template <size_t N, size_t M, typename T>
Matrix<N, M, T> operator*(const Matrix<N, M, T>&, const T&);





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


#endif