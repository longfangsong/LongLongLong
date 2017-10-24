#include "LongLongLong.h"
#include <algorithm>
#include <iterator>
#include <cctype>

using std::sqrt;
using std::deque;

std::complex<double> unity_root(size_t n, int k) {
    static const auto base = std::exp(std::complex<double>(0, 1) * M_PI);
    return std::pow(base, 2.0 * k / static_cast<double>(n));
}

deque<std::complex<double> > fast_Fourier_transform(std::deque<std::complex<double> > v) {
    if (v.size() == 1)
        return v;
    else {
        std::deque<std::complex<double> > even, odd;
        unsigned sz;
        for (sz = 0; sz < v.size(); ++sz) {
            if (sz % 2 == 1)
                odd.push_back(v[sz]);
            else
                even.push_back(v[sz]);
        }
        auto odd_transformed = fast_Fourier_transform(odd),
                even_transformed = fast_Fourier_transform(even);
        std::deque<std::complex<double> > ans(sz);
        for (int i = 0; i < sz / 2; ++i) {
            ans[i] = even_transformed[i] + unity_root(sz, i) * odd_transformed[i];
            ans[i + sz / 2] = even_transformed[i] - unity_root(sz, i) * odd_transformed[i];
        }
        return ans;
    }
}

deque<std::complex<double> > inverse_fast_Fourier_transform(std::deque<std::complex<double> > v) {
    if (v.size() == 1)
        return v;
    else {
        std::deque<std::complex<double> > even, odd;
        unsigned sz;
        for (sz = 0; sz < v.size(); ++sz) {
            if (sz % 2 == 1)
                odd.push_back(v[sz]);
            else
                even.push_back(v[sz]);
        }
        auto odd_transformed = inverse_fast_Fourier_transform(odd),
                even_transformed = inverse_fast_Fourier_transform(even);
        std::deque<std::complex<double> > ans(sz);
        for (int i = 0; i < sz / 2; ++i) {
            ans[i] = (even_transformed[i] + unity_root(sz, -i) * odd_transformed[i]);
            ans[i + sz / 2] = (even_transformed[i] - unity_root(sz, -i) * odd_transformed[i]);
        }
        return ans;
    }
}

template<>
deque<std::complex<double> >
discrete_Fourier_transform<std::complex<double> >(std::deque<std::complex<double> > v) {
    while ((v.size() - 1) & v.size())
        v.push_back(0);
    deque<std::complex<double> > ans;
    auto fft_result = fast_Fourier_transform(v);
    std::copy_n(fft_result.begin(), v.size(), std::back_inserter(ans));
    return ans;
}

template<typename T>
deque<std::complex<double> > discrete_Fourier_transform(std::deque<T> v) {
    std::deque<std::complex<double> > v_cmp;
    std::transform(v.begin(), v.end(), std::back_inserter(v_cmp), [](T n) {
        return std::complex<double>(n);
    });
    return discrete_Fourier_transform(v_cmp);
}

deque<std::complex<double> > inverse_discrete_Fourier_transform(std::deque<std::complex<double> > v) {
    while ((v.size() - 1) & v.size())
        v.push_back(0);
    deque<std::complex<double> > ans;
    auto fft_result = inverse_fast_Fourier_transform(v);
    std::copy_n(fft_result.begin(), v.size(), std::back_inserter(ans));
    std::for_each(ans.begin(), ans.end(), [&](std::complex<double> &c) {
        c /= v.size();
    });
    return ans;
}

LongLongLong::LongLongLong() : m_Numbers(deque<int>(1, 0)), m_Negative(false) {}

LongLongLong::LongLongLong(std::string s) {
    auto it = s.begin();
    if (*it == '-') {
        m_Negative = true;
        ++it;
    } else if (*it == '+') {
        m_Negative = false;
        ++it;
    } else {
        m_Negative = false;
    }
    for (; it != s.end(); ++it) {
        m_Numbers.push_back(static_cast<int> (*it - '0'));
    }
    std::reverse(m_Numbers.begin(), m_Numbers.end());
}

LongLongLong::LongLongLong(std::deque<int> numbers) : m_Negative(false) {
    m_Numbers = numbers;
}

template<typename Iterator>
LongLongLong::LongLongLong(Iterator begin, Iterator end, bool should_reverse):m_Negative(false) {
    std::copy(begin, end, std::back_inserter(m_Numbers));
    if (should_reverse) {
        std::reverse(m_Numbers.begin(), m_Numbers.end());
    }
}

LongLongLong::LongLongLong(const int &n) {
    std::ostringstream ss;
    ss << n;
    *this = LongLongLong(ss.str());
}

LongLongLong &LongLongLong::operator+=(const LongLongLong &other) {
    if (m_Negative && other.m_Negative) {
        *this = -((-*this) + (-other));
    } else if (m_Negative) {
        *this = other - (-*this);
    } else if (other.m_Negative) {
        *this -= (-other);
    } else if (other.m_Numbers.size() > m_Numbers.size()) {
        *this = other + *this;
    } else {
        auto it1 = m_Numbers.begin();
        auto it2 = other.m_Numbers.begin();
        for (;
                it2 != other.m_Numbers.end(); ++it1, ++it2) {
            *it1 += *it2;
        }
        for (it1 = m_Numbers.begin(); it1 < m_Numbers.end() - 1; ++it1) {
            if (*it1 >= 10) {
                *(it1 + 1) += *it1 / 10;
                *(it1) %= 10;
            }
        }
        if (m_Numbers.back() >= 10) {
            m_Numbers.push_back(m_Numbers.back() / 10);
            *(m_Numbers.rbegin() + 1) %= 10;
        }
    }
    return *this;
}

LongLongLong &LongLongLong::operator-=(const LongLongLong &other) {
    if (other.m_Negative) {
        *this += -(other);
    } else if (*this < other) {
        *this = -(other - *this);
    } else {
        auto it1 = m_Numbers.begin();
        auto it2 = other.m_Numbers.begin();
        for (;
                it2 != other.m_Numbers.end(); ++it1, ++it2) {
            *it1 -= *it2;
        }
        for (auto it = m_Numbers.begin(); it != m_Numbers.end() - 1; ++it) {
            while(*it < 0) {
                *(it+1) -= 1;
                *it += 10;
            }
        }
        while (m_Numbers.back() == 0) {
            m_Numbers.pop_back();
        }
    }
    return *this;
}

LongLongLong &LongLongLong::operator*=(const LongLongLong &other) {
    if (*this == LongLongLong() || other == LongLongLong()) {
        *this = LongLongLong();
    } else if (other.m_Numbers.size() < m_Numbers.size()) {
        *this = other * (*this);
    } else if (other == LongLongLong("10")) {
        m_Numbers.push_front(0);
    } else if (other.m_Numbers.size() == 1) {
        m_Negative = m_Negative ^ other.m_Negative;
        int n = other.m_Numbers[0];
        for_each(m_Numbers.begin(), m_Numbers.end(), [&](int &val) {
            val *= n;
        });
        for (auto it = m_Numbers.begin(); it < m_Numbers.end() - 1; ++it) {
            if (*it > 10) {
                *(it + 1) += *it / 10;
                *(it) %= 10;
            }
        }
        if (m_Numbers.back() > 10) {
            m_Numbers.push_back(m_Numbers.back() / 10);
            *(m_Numbers.rbegin() + 1) %= 10;
        }
    } else {
        m_Negative = m_Negative ^ other.m_Negative;
        auto other_numbers = other.m_Numbers;
        for (size_t i = m_Numbers.size(); i < other_numbers.size(); ++i) {
            m_Numbers.push_back(0);
        }
        auto sz = other_numbers.size() * 2;
        for (size_t i = other_numbers.size(); i < sz; ++i) {
            m_Numbers.push_back(0);
            other_numbers.push_back(0);
        }
        auto ptv1 = discrete_Fourier_transform(m_Numbers),
                ptv2 = discrete_Fourier_transform(other_numbers);

        for (auto it1 = ptv1.begin(), it2 = ptv2.begin();
             it1 != ptv1.end(); ++it1, ++it2) {
            *it1 *= *it2;
        }
        auto coef = inverse_discrete_Fourier_transform(ptv1);
        m_Numbers.clear();
        std::transform(coef.begin(), coef.end(), std::back_inserter(m_Numbers), [](std::complex<double> cmp) {
            return floor(cmp.real() + 0.5);
        });
        for (auto it = m_Numbers.begin(); it < m_Numbers.end() - 1; ++it) {
            if (*it >= 10) {
                *(it + 1) += *it / 10;
                *(it) %= 10;
            }
        }
        if (m_Numbers.back() >= 10) {
            m_Numbers.push_back(m_Numbers.back() / 10);
            *(m_Numbers.rbegin() + 1) %= 10;
        }
        while (m_Numbers.back() == 0) {
            m_Numbers.pop_back();
        }
    }
    return *this;
}

LongLongLong &LongLongLong::operator/=(const LongLongLong &other) {
    if (other.m_Numbers.size() > m_Numbers.size()) {
        *this = LongLongLong();
    } else {
        m_Negative = m_Negative ^ other.m_Negative;
        LongLongLong ans;
        auto left_side = m_Numbers.rbegin();
        auto it = m_Numbers.rbegin() + other.m_Numbers.size();
        auto num = LongLongLong(left_side, it, true);
        while (it < m_Numbers.rend()) {
            if (num > other) {
                int try_divide;
                if (num.m_Numbers.size() == other.m_Numbers.size()) {
                    try_divide = num.m_Numbers.back() / other.m_Numbers.back();
                } else {
                    try_divide = (num.m_Numbers.back() * 10 + *(num.m_Numbers.rbegin() + 1)) / other.m_Numbers.back();
                }
                while (other * try_divide > num) {
                    --try_divide;
                }
                num -= other * try_divide;
                ans.m_Numbers.push_back(try_divide);
            } else if (num == other) {
                ans.m_Numbers.push_back(1);
                num = LongLongLong();
            } else {
                ans.m_Numbers.push_back(0);
            }
            num *= LongLongLong("10");
            num += LongLongLong(*it);
            ++it;
        }
        if (num > other) {
            int try_divide;
            if (num.m_Numbers.size() == other.m_Numbers.size()) {
                try_divide = num.m_Numbers.back() / other.m_Numbers.back();
            } else {
                try_divide = (num.m_Numbers.back() * 10 + *(num.m_Numbers.rbegin() + 1)) / other.m_Numbers.back();
            }
            while (other * try_divide > num) {
                --try_divide;
            }
            num -= other * try_divide;
            ans.m_Numbers.push_back(try_divide);
        } else if (num == other) {
            ans.m_Numbers.push_back(1);
            num = LongLongLong();
        } else {
            ans.m_Numbers.push_back(0);
        }
        std::reverse(ans.m_Numbers.begin(), ans.m_Numbers.end());
        m_Numbers = ans.m_Numbers;
        while (m_Numbers.back() == 0) {
            m_Numbers.pop_back();
        }
    }
    return *this;
}

LongLongLong &LongLongLong::operator%=(const LongLongLong &other) {
    *this -= *this / other * other;
    return *this;
}

LongLongLong operator+(const LongLongLong &n1, const LongLongLong &n2) {
    LongLongLong ans(n1);
    ans += n2;
    return ans;
}

LongLongLong operator-(const LongLongLong &n1, const LongLongLong &n2) {
    LongLongLong ans(n1);
    ans -= n2;
    return ans;
}

LongLongLong operator*(const LongLongLong &n1, const LongLongLong &n2) {
    LongLongLong ans(n1);
    ans *= n2;
    return ans;
}

LongLongLong operator/(const LongLongLong &n1, const LongLongLong &n2) {
    LongLongLong ans(n1);
    ans /= n2;
    return ans;
}

LongLongLong operator%(const LongLongLong &n1, const LongLongLong &n2) {
    LongLongLong ans(n1);
    ans %= n2;
    return ans;
}

LongLongLong LongLongLong::operator-() const {
    static LongLongLong zero;
    LongLongLong ans(*this);
    if (ans != zero)
        ans.m_Negative = !m_Negative;
    return ans;
}

bool operator==(const LongLongLong &n1, const LongLongLong &n2) {
    return n1.m_Negative == n2.m_Negative && n1.m_Numbers.size() == n2.m_Numbers.size()
           && std::equal(n1.m_Numbers.begin(), n1.m_Numbers.end(), n2.m_Numbers.begin());
}

bool operator!=(const LongLongLong &n1, const LongLongLong &n2) {
    return !(n1 == n2);
}

bool operator<(const LongLongLong &n1, const LongLongLong &n2) {
    if (n1.m_Negative && n2.m_Negative)
        return !(-n1 < -n2);
    else if (n1.m_Negative && !n2.m_Negative)
        return true;
    else if (!n1.m_Negative && n2.m_Negative)
        return false;
    else if (n1.m_Numbers.size() < n2.m_Numbers.size())
        return true;
    else if (n1.m_Numbers.size() > n2.m_Numbers.size())
        return false;
    else {
        for (auto it1 = n1.m_Numbers.rbegin(), it2 = n2.m_Numbers.rbegin();
             it1 != n1.m_Numbers.rend(); ++it1, ++it2) {
            if (*it1 < *it2)
                return true;
            else if (*it1 > *it2)
                return false;
        }
    }
    return false;
}

bool operator>(const LongLongLong &n1, const LongLongLong &n2) {
    return !(n1 <= n2);
}

bool operator<=(const LongLongLong &n1, const LongLongLong &n2) {
    return n1 < n2 || n1 == n2;
}

bool operator>=(const LongLongLong &n1, const LongLongLong &n2) {
    return n1 > n2 || n1 == n2;
}

LongLongLong::operator std::string() const {
    std::ostringstream ss;
    if (*this == LongLongLong("0")) {
        ss << '0';
    } else {
        if (m_Negative) {
            ss << '-';
        }
        for (auto it = m_Numbers.rbegin(); it != m_Numbers.rend(); ++it) {
            ss << *it;
        }
    }
    return ss.str();
}

std::ostream &operator<<(std::ostream &os, const LongLongLong &n) {
    os << std::string(n);
    return os;
}

std::istream &operator>>(std::istream &is, LongLongLong &n) {
    std::string s;
    std::istreambuf_iterator<char> ch(is);
    for (; isspace(*ch); ++ch);
    n.m_Negative = false;
    if (*ch == '-') {
        n.m_Negative = true;
        ++ch;
    } else if (*ch == '+') {
        ++ch;
    }
    n.m_Numbers.clear();
    for (; isdigit(*ch); ++ch) {
        n.m_Numbers.push_back(static_cast<int> (*ch - '0'));
    }
    std::reverse(n.m_Numbers.begin(), n.m_Numbers.end());
    return is;
}

