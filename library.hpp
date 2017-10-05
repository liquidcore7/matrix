#ifndef MATRIX_LIBRARY_H
#define MATRIX_LIBRARY_H

#include <array>
#include <sstream>

template <size_t N, size_t M, typename T>
//  NxM matrix of elements of type T
class Matrix
{
protected:
    std::array< std::array<T, M>, N > _grid;

public:
    Matrix() = default;
    explicit Matrix(std::istream&);
    explicit Matrix(const std::string&);
    explicit Matrix(const T&);
    explicit Matrix(const std::array< std::array<T, M>, N >& a)
    : _grid(a) {}
    explicit Matrix(std::array< std::array<T, M>, N >&& a)
    : _grid(std::move(a)) {}

    Matrix(const Matrix&) = default;
    Matrix(Matrix&& r) noexcept : _grid(std::move(r._grid))
    {}

    Matrix& operator=(const Matrix&) = default;
    Matrix& operator=(Matrix&& r) noexcept
    {
        this->_grid = std::move(r._grid);
    }

    inline T get(size_t i, size_t j) const
    { return _grid[i][j]; }

    inline T& get(size_t i, size_t j)
    { return _grid[i][j]; }

    inline decltype(_grid) asArray() const
    { return _grid; }


    Matrix<M, N, T> transpose() const;


    template <typename Fnc>
    // functor should be of the following signature: void functor(T&)
    // example: [] (T& a) {a += 15}
    void apply(const Fnc&);

    Matrix& operator+=(const Matrix&);
    Matrix& operator*=(const T&);

    // requires T to be an integral type
    Matrix& operator!()
    {*this *= -1;}

    // seems too useless, rather delete nor implement
    // Matrix& operator*=(const Matrix<M, M, T>&);

    operator std::string() const;

    template <size_t Nf, size_t Mf, typename Tf>
    // workaround required to avoid 'undefined reference' or '-Wno-non-template-friend', don`t touch
    friend std::ostream& operator<<(std::ostream&, const Matrix<Nf, Mf, Tf>&);

};

template <size_t N, typename T>
struct SqMatrix : public Matrix<N, N, T>
{
    SqMatrix() = default;
    explicit SqMatrix(std::istream& is) : Matrix<N, N, T>(is) {}
    explicit SqMatrix(const std::string& str) : Matrix<N, N, T>(str) {}
    explicit SqMatrix(const T& val) : Matrix<N, N, T>(val) {}
    explicit SqMatrix(const std::array< std::array<T, N>, N >& a)
            : Matrix<N, N, T>(a) {}
    explicit SqMatrix(std::array< std::array<T, N>, N >&& a)
    : Matrix<N, N, T>(std::move(a)) {}
    SqMatrix(const Matrix<N, N, T>& m) : Matrix<N, N, T>(m) {}

    SqMatrix(const SqMatrix&r) = default;
    SqMatrix(SqMatrix&& r) noexcept = default;

    SqMatrix& operator=(const SqMatrix&) = default;
    SqMatrix& operator=(SqMatrix&& r) noexcept = default;


    SqMatrix inverse() const;
    virtual T det() const;
};

template <size_t N, typename T>
struct TriangleMatrix : public SqMatrix<N, T>
{
    TriangleMatrix() : SqMatrix<N, T>() {}
    explicit TriangleMatrix(const std::array< std::array<T, N>, N >& a)
            : SqMatrix<N, T>(a) {}
    explicit TriangleMatrix(std::array< std::array<T, N>, N >&& a)
    : SqMatrix<N, T>(std::move(a)) {}

    TriangleMatrix(const SqMatrix<N, T>& sq);

    T det() const override;

};


namespace {

    template<size_t N, typename T>
    class Det
    {
        SqMatrix<N, T> _data;

    public:

        explicit Det(const SqMatrix<N, T> &m) : _data(m) {}

        explicit Det(SqMatrix<N, T> &&m) : _data(std::move(m)) {}

        T _minor(size_t i, size_t j) const;

        T operator()() const;
    }; // class Det


    template<size_t N, typename T>
    T Det<N, T>::_minor(size_t i, size_t j) const
    {
        std::array<std::array<T, N - 1>, N - 1> newMatrix;
        for (int im = 0, ishift = 0; im < N; ++im) {
            if (im == i) {
                ishift = 1;
                continue;
            }
            for (int jm = 0, jshift = 0; jm < N; ++jm) {
                if (jm == j) {
                    jshift = 1;
                    continue;
                }
                newMatrix[im - ishift][jm - jshift] = _data.get(im, jm);
            }
        }
        return Det<N - 1, T>(std::move(SqMatrix<N - 1, T>(newMatrix))) ();
    }

    template <typename T>
    class Det<2, T>
    {
        SqMatrix<2, T> _data;
    public:
        explicit Det(const SqMatrix<2, T> &m) : _data(m) {}

        explicit Det(SqMatrix<2, T> &&m) : _data(std::move(m)) {}

        T _minor(size_t i, size_t j) const;

        T operator()() const;
    };

    template <typename T>
    T Det<2, T>::operator()() const
    {
        return _data.get(0, 0) * _data.get(1, 1) - _data.get(0, 1) * _data.get(1, 0);
    }

    template <typename T>
    T Det<2, T>::_minor(size_t i, size_t j) const
    {
        return _data.get(
                i == 0 ? 1 : 0,
                j == 0 ? 1 : 0
        );
    };


    template<size_t N, typename T>
    T Det<N, T>::operator()() const {
        T accu = 0;
        int sign = 1;
        for (size_t j = 0; j < N; ++j, sign *= -1)
            accu += sign * _data.get(0, j) * _minor(0, j);
        return accu;
    }

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
Matrix<N, M, T>::operator std::string() const
{
    std::stringstream ss;
    ss << *this;
    return ss.str();
}

template <size_t N, size_t M, typename T>
std::istream& operator>>(std::istream &is, Matrix<N, M, T> &obj)
{
    obj.apply([&is] (T& cell) {is >> cell;});
    return is;
}

template <size_t N, size_t M, typename T>
std::ostream& operator<<(std::ostream &out, const Matrix<N, M, T> &obj)
{
    for (const auto& row : obj._grid)
    {
        for (const auto& cell : row)
            out << cell << ' ';
        out << '\n';
    }
    return out;
}

template <size_t N, size_t M, typename T>
Matrix<M, N, T> Matrix<N, M, T>::transpose() const
{
    Matrix<M, N, T> transposed;
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < M; ++j)
            transposed._grid[j][i] = this->_grid[i][j];
    return transposed;
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
            local.get(i,j) = lhs.get(i,j) * scalar;
    return local;
};

template <size_t N, size_t L, size_t M, typename T>
Matrix<N, M, T> operator*(const Matrix<N, L, T>& lhs, const Matrix<L, M, T>& rhs)
{
    Matrix<N, M, T> local(0);
    for (size_t i = 0; i < N; ++i)
        for (size_t j = 0; j < M; ++j)
        {
            for (size_t k = 0; k < L; ++k)
                local.get(i, j) += lhs.get(i, k) * rhs.get(k, j);
        }

    return local;
};

template <size_t N, typename T>
SqMatrix<N, T> operator*(const SqMatrix<N, T>& l, const SqMatrix<N, T>& r)
{
    SqMatrix<N, T> local(0);
    for (size_t i = 0; i < N; ++i)
        for (size_t j = 0; j < N; ++j)
        {
            for (size_t k = 0; k < N; ++k)
                local.get(i, j) += l.get(i, k) * r.get(k, j);
        }

    return local;
};


template <size_t N, typename T>
T SqMatrix<N, T>::det() const
{
    return Det<N, T>(*this)();
}

template <size_t N, typename T>
SqMatrix<N, T> SqMatrix<N, T>::inverse() const
{
    T det = this->det();
    if (det == 0)
        throw std::logic_error("Inversion on matrix A with detA = 0");
    SqMatrix<N, T> inverted;
    int sign = 1;
    for (size_t i = 0; i < N; ++i)
        for (size_t j = 0; j < N; ++j, sign *= -1)
            inverted.get(i, j) = sign * Det<N, T>(*this)._minor(i, j);
    return (inverted * (1 / this->det()) ).transpose();
};


template <size_t N, typename T>
TriangleMatrix<N, T>::TriangleMatrix(const SqMatrix<N, T> &sq)
{

    auto arr(std::move(sq.asArray()));

    size_t startJ = 0;
    for (size_t k = 0; k < N - 1; ++k, ++startJ)
        for (size_t i = k + 1; i < N; ++i)
        {
            double mult = static_cast<double>(arr[i][startJ]) / arr[k][startJ] * -1.0;
            for (size_t j = startJ; j < N; ++j)
                arr[i][j] += mult * arr[k][j];
        }
    this->_grid = arr;
}

template <size_t N, typename T>
T TriangleMatrix<N, T>::det() const
{
    T accu = 0;
    for (int d = 0; d < N; ++d)
        accu *= this->get(d, d);
    return accu;
}


#endif