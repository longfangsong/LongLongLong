#ifndef LONGLONGLONG_LIBRARY_H
#define LONGLONGLONG_LIBRARY_H

#include <deque>
#include <string>
#include <complex>
#include <memory>
#include <iostream>

template<typename T>
std::deque<std::complex<double>> discrete_Fourier_transform(std::deque<T> v);

template<>
std::deque<std::complex<double>> discrete_Fourier_transform<std::complex<double>>(std::deque<std::complex<double>> v);

std::deque<std::complex<double>> inverse_discrete_Fourier_transform(std::deque<std::complex<double>> v);

class LongLongLong {
public:
    LongLongLong();

    LongLongLong(const LongLongLong &other) = default;

    template<typename Iterator>
    LongLongLong(Iterator begin, Iterator end, bool should_reverse = false);

    LongLongLong(const int &n);

    LongLongLong(std::string s);

    friend LongLongLong operator+(const LongLongLong &n1, const LongLongLong &n2);

    friend LongLongLong operator-(const LongLongLong &n1, const LongLongLong &n2);

    friend LongLongLong operator*(const LongLongLong &n1, const LongLongLong &n2);

    friend LongLongLong operator/(const LongLongLong &n1, const LongLongLong &n2);

    friend LongLongLong operator%(const LongLongLong &n1, const LongLongLong &n2);

    LongLongLong &operator+=(const LongLongLong &other);

    LongLongLong &operator-=(const LongLongLong &other);

    LongLongLong &operator*=(const LongLongLong &other);

    LongLongLong &operator/=(const LongLongLong &other);

    LongLongLong &operator%=(const LongLongLong &other);

    friend bool operator==(const LongLongLong &n1, const LongLongLong &n2);

    friend bool operator!=(const LongLongLong &n1, const LongLongLong &n2);

    friend bool operator<(const LongLongLong &n1, const LongLongLong &n2);

    friend bool operator>(const LongLongLong &n1, const LongLongLong &n2);

    friend bool operator<=(const LongLongLong &n1, const LongLongLong &n2);

    friend bool operator>=(const LongLongLong &n1, const LongLongLong &n2);

    LongLongLong operator-() const;

    operator std::string() const;

    friend std::ostream &operator<<(std::ostream &os, const LongLongLong &n);

    friend std::istream &operator>>(std::istream &is, LongLongLong &n);

private:
    LongLongLong(std::deque<int> numbers);

    typedef std::complex<double> cmpl_t;
    std::deque<int> m_Numbers;
    bool m_Negative;
};

#endif