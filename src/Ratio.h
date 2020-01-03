// Copyright Alexander Bolz 2019
// Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at https://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <cassert>
#include <cstdint>
#include <climits>
#include <type_traits>

#ifndef RATIO_ASSERT
#define RATIO_ASSERT(X) assert(X)
#endif

//==================================================================================================
// Rational
//  Compile-time rational arithmetic
//==================================================================================================

namespace sc {

namespace impl
{
    constexpr intmax_t Abs(intmax_t x) noexcept {
        RATIO_ASSERT(x != INTMAX_MIN);
        return x < 0 ? -x : x;
    }

    constexpr intmax_t Sgn(intmax_t x) noexcept {
        return x < 0 ? -1 : 1;
    }

    constexpr intmax_t Gcd(intmax_t a, intmax_t b) noexcept {
        intmax_t x = Abs(a);
        intmax_t y = Abs(b);
        while (y > 0) {
            const auto r = x % y;
            x = y;
            y = r;
        }
        return x;
    }

    constexpr intmax_t Lcm(intmax_t a, intmax_t b) noexcept {
        return (a / Gcd(a, b)) * b;
    }

    constexpr bool CheckMul64(intmax_t x, intmax_t y) noexcept {
        return Abs(x) <= INTMAX_MAX / (y == 0 ? 1 : Abs(y));
    }

    constexpr bool CheckAdd64(intmax_t x, intmax_t y) noexcept {
        return Sgn(x) != Sgn(y) || (Abs(x) <= INTMAX_MAX - Abs(y));
    }

#if 0
    constexpr intmax_t Power(intmax_t x, intmax_t n) noexcept {
        RATIO_ASSERT(x >= 1);
        RATIO_ASSERT(n >= 0);

        intmax_t p = 1;
        if (x > 1) {
            for ( ; n > 0; --n) {
                RATIO_ASSERT(p <= INTMAX_MAX / x);
                p *= x;
            }
        }

        return p;
    }

    // Returns x^n > lower (without overflow).
    constexpr bool IsPowerGreaterThan(intmax_t x, intmax_t n, intmax_t lower) noexcept {
        RATIO_ASSERT(x >= 1);
        RATIO_ASSERT(n >= 0);
        RATIO_ASSERT(lower >= 0);

        const intmax_t lim = INTMAX_MAX / x;
        const intmax_t max = lim < lower ? lim : lower; // = min(lim, lower)

        intmax_t p = 1;
        for ( ; n > 0; --n) {
            if (p > max) // p*x will overflow, or p > lower
                return true;
            p *= x;
        }

        return p > lower;
    }

    // Computes the n-th root y of x,
    // i.e. returns the largest y, such that y^n <= x
    constexpr intmax_t Root(intmax_t x, intmax_t n) noexcept {
        RATIO_ASSERT(x >= 0);
        RATIO_ASSERT(n >= 1);

        if (x <= 1 || n <= 1)
            return x;

        intmax_t lo = 1;
        intmax_t hi = 1 + x / n;
        // Bernoulli  ==>  x^(1/n) <= 1 + (x - 1)/n < 1 + x/n
        // Since n >= 2, hi will not overflow here.

        for (;;) {
            const intmax_t y = lo + (hi - lo) / 2;
            if (y == lo)                           // hi - lo <= 1
                return y;
            else if (IsPowerGreaterThan(y, n, x))  // y^n > x
                hi = y;
            else                                   // y^n <= x
                lo = y;
        }
    }
#endif
}

template <intmax_t Num, intmax_t Den>
struct Rational
{
    static_assert(Den != 0,
        "invalid operation");
    static_assert(-INTMAX_MAX <= Num,
        "invalid argument");
    static_assert(-INTMAX_MAX <= Den,
        "invalid argument");

    static constexpr intmax_t num = impl::Abs(Num) / impl::Gcd(Num, Den) * (impl::Sgn(Num) * impl::Sgn(Den));
    static constexpr intmax_t den = impl::Abs(Den) / impl::Gcd(Num, Den);
    using type = Rational<num, den>;

    template <typename ValueT>
    [[nodiscard]] constexpr auto operator()(ValueT x) noexcept
    {
        using ResultT = std::common_type_t<std::remove_cv_t<ValueT>, intmax_t>;

        if constexpr (den == 1)
        {
            if constexpr (num == 1)
                return static_cast<ResultT>(x);
            else
                return static_cast<ResultT>(x) * static_cast<ResultT>(num);
        }
        else
        {
            if constexpr (num == 1)
                return static_cast<ResultT>(x) / static_cast<ResultT>(den);
            else
                return static_cast<ResultT>(x) * static_cast<ResultT>(num) / static_cast<ResultT>(den);
        }
    }

    template <intmax_t Num2, intmax_t Den2>
    [[nodiscard]] constexpr friend auto operator+(Rational, Rational<Num2, Den2>) noexcept
    {
        constexpr intmax_t N1 = num;
        constexpr intmax_t D1 = den;
        constexpr intmax_t N2 = Rational<Num2, Den2>::num;
        constexpr intmax_t D2 = Rational<Num2, Den2>::den;

        constexpr intmax_t Gx = impl::Gcd(D1, D2);

        static_assert(impl::CheckMul64(N1, D2 / Gx),
            "integer arithmetic overflow");
        static_assert(impl::CheckMul64(N2, D1 / Gx),
            "integer arithmetic overflow");
        static_assert(impl::CheckAdd64(N1 * (D2 / Gx), N2 / (D1 / Gx)),
            "integer arithmetic overflow");
        static_assert(impl::CheckMul64(D1, D2 / Gx),
            "integer arithmetic overflow");

        constexpr intmax_t Nr = N1 * (D2 / Gx) + N2 * (D1 / Gx);
        constexpr intmax_t Dr = D1 * (D2 / Gx);

        return typename Rational<Nr, Dr>::type{}; // Nr/Dr needs to be reduced again!
    }

    template <intmax_t Num2, intmax_t Den2>
    [[nodiscard]] constexpr friend auto operator-(Rational lhs, Rational<Num2, Den2>) noexcept
    {
        using RHS = Rational<-Num2, Den2>; // No need to reduce again...
        return lhs * RHS{};
    }

    template <intmax_t Num2, intmax_t Den2>
    [[nodiscard]] constexpr friend auto operator*(Rational, Rational<Num2, Den2>) noexcept
    {
        constexpr intmax_t N1 = num;
        constexpr intmax_t D1 = den;
        constexpr intmax_t N2 = Rational<Num2, Den2>::num;
        constexpr intmax_t D2 = Rational<Num2, Den2>::den;

        constexpr intmax_t Gx = impl::Gcd(N1, D2);
        constexpr intmax_t Gy = impl::Gcd(N2, D1);

        static_assert(impl::CheckMul64(N1 / Gx, N2 / Gy),
            "integer arithmetic overflow");
        static_assert(impl::CheckMul64(D1 / Gy, D2 / Gx),
            "integer arithmetic overflow");

        constexpr intmax_t Nr = (N1 / Gx) * (N2 / Gy);
        constexpr intmax_t Dr = (D1 / Gy) * (D2 / Gx);

        return Rational<Nr, Dr>{};
    }

    template <intmax_t Num2, intmax_t Den2>
    [[nodiscard]] constexpr friend auto operator/(Rational lhs, Rational<Num2, Den2>) noexcept
    {
        using RHS = typename Rational<Den2, Num2>::type; // Reduce!
        return lhs * RHS{};
    }

    template <intmax_t Num2, intmax_t Den2>
    [[nodiscard]] constexpr friend bool operator==(Rational, Rational<Num2, Den2>) noexcept
    {
        constexpr intmax_t N1 = num;
        constexpr intmax_t D1 = den;
        constexpr intmax_t N2 = Rational<Num2, Den2>::num;
        constexpr intmax_t D2 = Rational<Num2, Den2>::den;

        return N1 == N2 && D1 == D2;
    }

    template <intmax_t Num2, intmax_t Den2>
    [[nodiscard]] constexpr friend bool operator!=(Rational lhs, Rational<Num2, Den2> rhs) noexcept
    {
        return !(lhs == rhs);
    }

    //template <intmax_t Num2, intmax_t Den2>
    //[[nodiscard]] constexpr friend bool operator<(Rational, Rational<Num2, Den2>) noexcept
    //{
    //}

    //template <intmax_t Num2, intmax_t Den2>
    //[[nodiscard]] constexpr friend bool operator>(Rational, Rational<Num2, Den2>) noexcept
    //{
    //}

    //template <intmax_t Num2, intmax_t Den2>
    //[[nodiscard]] constexpr friend bool operator<=(Rational, Rational<Num2, Den2>) noexcept
    //{
    //}

    //template <intmax_t Num2, intmax_t Den2>
    //[[nodiscard]] constexpr friend bool operator>=(Rational, Rational<Num2, Den2>) noexcept
    //{
    //}
};

template <intmax_t Num, intmax_t Den = 1>
    using Ratio = typename Rational<Num, Den>::type;

template <typename R1, typename R2>
    using AddRatios = decltype(R1{} + R2{});

template <typename R1, typename R2>
    using SubRatios = decltype(R1{} - R2{});

template <typename R1, typename R2>
    using MulRatios = decltype(R1{} * R2{});

template <typename R1, typename R2>
    using DivRatios = decltype(R1{} / R2{});

template <typename R1, typename R2>
    using RatioDivides = std::bool_constant< DivRatios<R2, R1>::den == 1 >;

// CommonRatio = GCD(R1, R2) = (GCD(N1, N2), LCM(D1, D2))
template <typename R1, typename R2>
    using CommonRatio = Rational< impl::Gcd(R1::num, R2::num), impl::Lcm(R1::den, R2::den) >;

} // namespace sc

#if 0
// ratio_less

template <class _R1, class _R2, bool _Odd = false,
          intmax_t _Q1 = _R1::num / _R1::den, intmax_t _M1 = _R1::num % _R1::den,
          intmax_t _Q2 = _R2::num / _R2::den, intmax_t _M2 = _R2::num % _R2::den>
struct __ratio_less1
{
    static const bool value = _Odd ? _Q2 < _Q1 : _Q1 < _Q2;
};

template <class _R1, class _R2, bool _Odd, intmax_t _Qp>
struct __ratio_less1<_R1, _R2, _Odd, _Qp, 0, _Qp, 0>
{
    static const bool value = false;
};

template <class _R1, class _R2, bool _Odd, intmax_t _Qp, intmax_t _M2>
struct __ratio_less1<_R1, _R2, _Odd, _Qp, 0, _Qp, _M2>
{
    static const bool value = !_Odd;
};

template <class _R1, class _R2, bool _Odd, intmax_t _Qp, intmax_t _M1>
struct __ratio_less1<_R1, _R2, _Odd, _Qp, _M1, _Qp, 0>
{
    static const bool value = _Odd;
};

template <class _R1, class _R2, bool _Odd, intmax_t _Qp, intmax_t _M1,
                                                        intmax_t _M2>
struct __ratio_less1<_R1, _R2, _Odd, _Qp, _M1, _Qp, _M2>
{
    static const bool value = __ratio_less1<ratio<_R1::den, _M1>,
                                            ratio<_R2::den, _M2>, !_Odd>::value;
};

template <class _R1, class _R2, intmax_t _S1 = __static_sign<_R1::num>::value,
                                intmax_t _S2 = __static_sign<_R2::num>::value>
struct __ratio_less
{
    static const bool value = _S1 < _S2;
};

template <class _R1, class _R2>
struct __ratio_less<_R1, _R2, 1LL, 1LL>
{
    static const bool value = __ratio_less1<_R1, _R2>::value;
};

template <class _R1, class _R2>
struct __ratio_less<_R1, _R2, -1LL, -1LL>
{
    static const bool value = __ratio_less1<ratio<-_R2::num, _R2::den>, ratio<-_R1::num, _R1::den> >::value;
};
#endif
