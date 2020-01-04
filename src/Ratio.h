// Copyright Alexander Bolz 2019
// Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at https://www.boost.org/LICENSE_1_0.txt)

#pragma once

#define RATIO_COMPARE_LESS() 0
#define RATIO_ROOT() 0

#include <cassert>
#include <cstdint>
#include <climits>
#include <ratio>
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
    constexpr int64_t Abs(int64_t x) noexcept {
        RATIO_ASSERT(x != INT64_MIN);
        return x < 0 ? -x : x;
    }

    // Returns: x >= 0 ? 1 : -1
    constexpr int64_t Sgn(int64_t x) noexcept {
        return x < 0 ? -1 : 1;
    }

    constexpr int64_t Gcd(int64_t a, int64_t b) noexcept {
        int64_t x = Abs(a);
        int64_t y = Abs(b);
        while (y > 0) {
            const auto r = x % y;
            x = y;
            y = r;
        }
        return x;
    }

    constexpr int64_t Lcm(int64_t a, int64_t b) noexcept {
        return (a / Gcd(a, b)) * b;
    }

    constexpr bool CheckMul64(int64_t x, int64_t y) noexcept {
        return Abs(x) <= INT64_MAX / (y == 0 ? 1 : Abs(y));
    }

    constexpr bool CheckAdd64(int64_t x, int64_t y) noexcept {
        return Sgn(x) != Sgn(y) || (Abs(x) <= INT64_MAX - Abs(y));
    }

#if RATIO_COMPARE_LESS()
    struct Uint64x2 {
        uint64_t hi;
        uint64_t lo;

        constexpr friend bool operator==(Uint64x2 lhs, Uint64x2 rhs) noexcept {
            return lhs.hi == rhs.hi && lhs.lo == rhs.lo;
        }

        constexpr friend bool operator<(Uint64x2 lhs, Uint64x2 rhs) noexcept {
            return lhs.hi < rhs.hi || (lhs.hi == rhs.hi && lhs.lo < rhs.lo);
        }
    };

    constexpr uint32_t Lo32(uint64_t x) noexcept {
        return static_cast<uint32_t>(x & 0xFFFFFFFFu);
    }

    constexpr uint32_t Hi32(uint64_t x) noexcept {
        return static_cast<uint32_t>(x >> 32);
    }

    constexpr uint64_t Load64(uint32_t hi, uint32_t lo) noexcept {
        return uint64_t{hi} << 32 | lo;
    }

    constexpr Uint64x2 Mul64x64(uint64_t a, uint64_t b) noexcept {
        const uint64_t b00 = uint64_t{Lo32(a)} * Lo32(b);
        const uint64_t b01 = uint64_t{Lo32(a)} * Hi32(b);
        const uint64_t b10 = uint64_t{Hi32(a)} * Lo32(b);
        const uint64_t b11 = uint64_t{Hi32(a)} * Hi32(b);

        const uint64_t mid1 = b10 + Hi32(b00);
        const uint64_t mid2 = b01 + Lo32(mid1);

        const uint64_t hi = b11 + Hi32(mid1) + Hi32(mid2);
        const uint64_t lo = Load64(Lo32(mid2), Lo32(b00));
        return {hi, lo};
    }

    constexpr Uint64x2 Mul64x64(int64_t a, int64_t b) noexcept {
        RATIO_ASSERT(a >= 0);
        RATIO_ASSERT(b >= 0);
        return Mul64x64(static_cast<uint64_t>(a), static_cast<uint64_t>(b));
    }
#endif

#if RATIO_ROOT()
    // Computes y^n.
    // Does not check for overflow.
    constexpr int64_t Power(int64_t y, int64_t n) noexcept {
        RATIO_ASSERT(y >= 1);
        RATIO_ASSERT(n >= 0);

        int64_t p = 1;
        //if (y > 1) {
            for ( ; n > 0; --n) {
                RATIO_ASSERT(p <= INT64_MAX / y);
                p *= y;
            }
        //}

        return p;
    }

    // Returns y^n <=> x (without overflow).
    constexpr int ComparePower(int64_t y, int64_t n, int64_t x) noexcept {
        RATIO_ASSERT(y >= 1);
        RATIO_ASSERT(n >= 0);
        RATIO_ASSERT(x >= 0);

        const int64_t lim = INT64_MAX / y;
        const int64_t max = lim < x ? lim : x; // = min(lim, x)

        int64_t p = 1;
        for ( ; n > 0; --n) {
            if (p > max) // p*y will overflow, or p > x
                return +1;
            p *= y;
        }

        return p < x ? -1 : (x < p);
    }

    //struct RootResult
    //{
    //    int64_t value;
    //    bool is_exact;
    //};

    // Computes the n-th root y of x,
    // i.e. returns the largest y, such that y^n <= x
    constexpr int64_t Root(int64_t x, int64_t n) noexcept {
        RATIO_ASSERT(x >= 0);
        RATIO_ASSERT(n >= 1);

        if (x <= 1 || n <= 1)
            return x;

        int64_t lo = 1;
        int64_t hi = 1 + x / n;
        // Bernoulli  ==>  x^(1/n) <= 1 + (x - 1)/n < 1 + x/n
        // Since n >= 2, hi will not overflow here.

        for (;;) {
            const auto y = lo + (hi - lo) / 2;
            if (y == lo)        // hi - lo <= 1
                return y;
            const auto cmp = ComparePower(y, n, x);
            if (cmp < 0)        // y^n < x
                lo = y;
            else if (0 < cmp)   // y^n > x
                hi = y;
            else                // y^n = x
                return y;
        }
    }
#endif
}

template <int64_t Num, int64_t Den>
struct Rational
{
    static_assert(-INT64_MAX <= Num,
        "invalid argument");
    static_assert(-INT64_MAX <= Den,
        "invalid argument");
    static_assert(Den != 0,
        "invalid argument");

    static constexpr int64_t num = impl::Abs(Num) / impl::Gcd(Num, Den) * (impl::Sgn(Num) * impl::Sgn(Den));
    static constexpr int64_t den = impl::Abs(Den) / impl::Gcd(Num, Den);
    using type = Rational<num, den>;

    template <typename ValueT>
    [[nodiscard]] constexpr auto operator()(ValueT x) noexcept
    {
        using ResultT = std::common_type_t<std::remove_cv_t<ValueT>, int64_t>;

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

    template <int64_t Num2, int64_t Den2>
    [[nodiscard]] constexpr friend auto operator+(Rational, Rational<Num2, Den2>) noexcept
    {
        constexpr int64_t N1 = num;
        constexpr int64_t D1 = den;
        constexpr int64_t N2 = Rational<Num2, Den2>::num;
        constexpr int64_t D2 = Rational<Num2, Den2>::den;

        constexpr int64_t Gd = impl::Gcd(D1, D2);

        static_assert(impl::CheckMul64(N1, D2 / Gd),
            "integer arithmetic overflow");
        static_assert(impl::CheckMul64(N2, D1 / Gd),
            "integer arithmetic overflow");
        static_assert(impl::CheckAdd64(N1 * (D2 / Gd), N2 / (D1 / Gd)),
            "integer arithmetic overflow");
        static_assert(impl::CheckMul64(D1, D2 / Gd),
            "integer arithmetic overflow");

        constexpr int64_t Nr = N1 * (D2 / Gd) + N2 * (D1 / Gd);
        constexpr int64_t Dr = D1 * (D2 / Gd);

        return typename Rational<Nr, Dr>::type{}; // Nr/Dr needs to be reduced again!
    }

    template <int64_t Num2, int64_t Den2>
    [[nodiscard]] constexpr friend auto operator-(Rational lhs, Rational<Num2, Den2>) noexcept
    {
        using RHS = Rational<-Num2, Den2>; // No need to reduce again...
        return lhs * RHS{};
    }

    template <int64_t Num2, int64_t Den2>
    [[nodiscard]] constexpr friend auto operator*(Rational, Rational<Num2, Den2>) noexcept
    {
        constexpr int64_t N1 = num;
        constexpr int64_t D1 = den;
        constexpr int64_t N2 = Rational<Num2, Den2>::num;
        constexpr int64_t D2 = Rational<Num2, Den2>::den;

        constexpr int64_t Gx = impl::Gcd(N1, D2);
        constexpr int64_t Gy = impl::Gcd(N2, D1);

        static_assert(impl::CheckMul64(N1 / Gx, N2 / Gy),
            "integer arithmetic overflow");
        static_assert(impl::CheckMul64(D1 / Gy, D2 / Gx),
            "integer arithmetic overflow");

        constexpr int64_t Nr = (N1 / Gx) * (N2 / Gy);
        constexpr int64_t Dr = (D1 / Gy) * (D2 / Gx);

        return Rational<Nr, Dr>{};
    }

    template <int64_t Num2, int64_t Den2>
    [[nodiscard]] constexpr friend auto operator/(Rational lhs, Rational<Num2, Den2>) noexcept
    {
        using RHS = typename Rational<Den2, Num2>::type; // Reduce!
        return lhs * RHS{};
    }

    template <int64_t Num2, int64_t Den2>
    [[nodiscard]] constexpr friend bool operator==(Rational, Rational<Num2, Den2>) noexcept
    {
        constexpr int64_t N1 = num;
        constexpr int64_t D1 = den;
        constexpr int64_t N2 = Rational<Num2, Den2>::num;
        constexpr int64_t D2 = Rational<Num2, Den2>::den;

        return N1 == N2 && D1 == D2;
    }

    template <int64_t Num2, int64_t Den2>
    [[nodiscard]] constexpr friend bool operator!=(Rational lhs, Rational<Num2, Den2> rhs) noexcept
    {
        return !(lhs == rhs);
    }

#if RATIO_COMPARE_LESS()
    template <int64_t Num2, int64_t Den2>
    [[nodiscard]] constexpr friend bool operator<(Rational, Rational<Num2, Den2>) noexcept
    {
        constexpr int64_t N1 = num;
        constexpr int64_t D1 = den;
        constexpr int64_t N2 = Rational<Num2, Den2>::num;
        constexpr int64_t D2 = Rational<Num2, Den2>::den;

        if constexpr (impl::Sgn(N1) != impl::Sgn(N2))
            return N1 < N2;

        if constexpr (N1 >= 0)
            return impl::Mul64x64( N1, D2) < impl::Mul64x64( N2, D1);
        else
            return impl::Mul64x64(-N2, D1) < impl::Mul64x64(-N1, D2);
    }

    template <int64_t Num2, int64_t Den2>
    [[nodiscard]] constexpr friend bool operator>(Rational lhs, Rational<Num2, Den2> rhs) noexcept
    {
        return rhs < lhs;
    }

    template <int64_t Num2, int64_t Den2>
    [[nodiscard]] constexpr friend bool operator<=(Rational lhs, Rational<Num2, Den2> rhs) noexcept
    {
        return !(rhs < lhs);
    }

    template <int64_t Num2, int64_t Den2>
    [[nodiscard]] constexpr friend bool operator>=(Rational lhs, Rational<Num2, Den2> rhs) noexcept
    {
        return !(lhs < rhs);
    }
#endif
};

template <int64_t Num, int64_t Den = 1>
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
