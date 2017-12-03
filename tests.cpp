//
// Created by liquidcore7 on 30.11.17.
//

#define BOOST_TEST_MODULE MatrixTest

#include <boost/test/unit_test.hpp>
#include "library.hpp"
#include <array>

//TODO: do something with matrix equality tests
//e.g. cases like BOOST_TEST( (matrix1 == matrix2) == true )

BOOST_AUTO_TEST_CASE(from_string_constructable)
{
    Matrix<2, 3, int> t("1 2 3 4 5 6");
    BOOST_TEST(static_cast<std::string>(t) == "1 2 3 \n4 5 6 \n");
    SqMatrix<2, int> s("2 2 4 5");
    BOOST_TEST(static_cast<std::string>(s) == "2 2 \n4 5 \n");
}

BOOST_AUTO_TEST_CASE(with_value_fillable)
{
    Matrix<3, 2, bool> b(true);
    BOOST_TEST(static_cast<std::string>(b) == "1 1 \n1 1 \n1 1 \n");
    SqMatrix<2, char> c('a');
    BOOST_TEST(static_cast<std::string>(c) == "a a \na a \n");
}

BOOST_AUTO_TEST_CASE(from_array_constructable)
{
    std::array<std::array<double, 5>, 2> mmap =
            {{
                     {{1.1, 1.2, 1.3, 1.4, 1.5}},
                     {{2.1, 2.2, 2.3, 2.4, 2.5}}
             }};
    Matrix<2, 5, double> d(mmap);
    BOOST_TEST(static_cast<std::string>(d) == "1.1 1.2 1.3 1.4 1.5 \n2.1 2.2 2.3 2.4 2.5 \n");

    Matrix<2, 5, double> dm(std::move(mmap));
    BOOST_TEST(static_cast<std::string>(dm) == "1.1 1.2 1.3 1.4 1.5 \n2.1 2.2 2.3 2.4 2.5 \n");

    BOOST_TEST(dm.asArray() == mmap);
}

BOOST_AUTO_TEST_CASE(matrix_copyable)
{
    SqMatrix<3, int> k("1 2 3 4 5 6 7 8 9");
    SqMatrix<3, int> cp = k;
    SqMatrix<3, int> mv(std::move(k));

    BOOST_TEST(static_cast<std::string>(cp) == "1 2 3 \n4 5 6 \n7 8 9 \n");
    BOOST_TEST(static_cast<std::string>(mv) == "1 2 3 \n4 5 6 \n7 8 9 \n");
}

BOOST_AUTO_TEST_CASE(getters_n_setters)
{
    Matrix<4, 3, short> a("1 3 4 6 2 4 6 7 8 1 2 3");

    BOOST_TEST(a.get(1, 2) == 4);
    BOOST_CHECK_THROW(a.get(10, 10), std::out_of_range);

    a.get(0, 0) = 15;
    a.get(0, 1) = 14;
    a.get(0, 2) = 13;

    std::array<short, 3> excepted = {15, 14, 13};

    BOOST_TEST(a.asArray()[0] == excepted);
}


BOOST_AUTO_TEST_CASE(transponation)
{
    SqMatrix<2, int> p("1 2 3 4");
    p = p.transpose();

    BOOST_TEST(static_cast<std::string>(p) == "1 3 \n2 4 \n");

    Matrix<1, 4, float > m("1.1 1.2 1.3 1.4");
    auto m_t = m.transpose();

    BOOST_TEST(static_cast<std::string>(m_t) == "1.1 \n1.2 \n1.3 \n1.4 \n");
}

BOOST_AUTO_TEST_CASE(functors)
{
    SqMatrix<8, int> five(5);
    five.apply([] (int& cell, const size_t i, const size_t j)
               {
                   if (i == j)
                       cell *= 2;
               });

    for (int k = 0; k < 8; ++k)
        BOOST_TEST(five.get(k, k) == 10);

    int autoIncr = 2;
    Matrix<3, 5, int> m(8);
    m.apply([&autoIncr] (int& cell, const size_t i, const size_t j)
            {
               if (i == 1)
                   cell = autoIncr++;
            });

    std::array<int, 5> excepted = {2, 3, 4, 5, 6};

    BOOST_TEST(m.asArray()[1] == excepted);
}

BOOST_AUTO_TEST_CASE(linear_operations)
{
    Matrix<2, 3, int> l(3);
    Matrix<2, 3, int> r(5);

    auto sum = l + r;

    BOOST_TEST((sum == Matrix<2, 3, int>(8)) == true);

    sum *= 3;

    BOOST_TEST((sum == Matrix<2, 3, int> (24)) == true);

    sum += !r; // sum + (-r)

    BOOST_TEST((sum == Matrix<2, 3, int> (19)) == true);

    SqMatrix<3, int> p(2);
    !p += p * 2 + p;

    BOOST_TEST((p == SqMatrix<3, int> (4)) == true);
}

BOOST_AUTO_TEST_CASE(stream_operators)
{
    std::string values = "1 2 3 4";
    std::istringstream ss(values);

    SqMatrix<2, int> fromSS(std::move(ss));

    std::ostringstream out;
    out << fromSS;

    BOOST_TEST(out.str() == static_cast<std::string>(fromSS));
}


BOOST_AUTO_TEST_CASE(multiplication)
{
    Matrix<2, 3, int> ones(1);
    Matrix<3, 2, int> twos(2);
    auto mult = ones * twos;

    BOOST_TEST((mult == SqMatrix<2, int>(6)) == true);
}

//TODO: add tests for det(),





