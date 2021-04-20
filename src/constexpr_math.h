// Copyright 2021 Alexander Bolz
//
// Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at https://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <cstdint>

#if defined(_MSC_VER) && !defined(__clang__) && _MSC_VER >= 1925
#	define UNITS_HAS_CONSTEXPR_MATH() 1
#elif defined(__has_builtin)
#	if __has_builtin(__builtin_bit_cast) && __has_builtin(__builtin_is_constant_evaluated)
#		define UNITS_HAS_CONSTEXPR_MATH() 1
#	else
#		define UNITS_HAS_CONSTEXPR_MATH() 0
#	endif
#else
#	define UNITS_HAS_CONSTEXPR_MATH() 0
#endif

#if UNITS_HAS_CONSTEXPR_MATH()
namespace uom::constexpr_math {

//==================================================================================================
//
//==================================================================================================

#if 0

using Float = float;
using Bits  = uint32_t;

inline constexpr int32_t SignificandBits = 24; // == p (includes the hidden bit)
inline constexpr int32_t ExponentBits    = 8;
inline constexpr int32_t ExponentBias    = 0x7F;
inline constexpr Bits    SignificandMask = 0x007FFFFF;
inline constexpr Bits    ExponentMask    = 0x7F800000;
inline constexpr Bits    SignMask        = 0x80000000;
inline constexpr Bits    HiddenBit       = 0x00800000;

#else

using Float = double;
using Bits  = uint64_t;

inline constexpr int32_t SignificandBits = 53; // == p (includes the hidden bit)
inline constexpr int32_t ExponentBits    = 11;
inline constexpr int32_t ExponentBias    = 0x3FF;
inline constexpr Bits    SignificandMask = 0x000FFFFFFFFFFFFF;
inline constexpr Bits    ExponentMask    = 0x7FF0000000000000;
inline constexpr Bits    SignMask        = 0x8000000000000000;
inline constexpr Bits    HiddenBit       = 0x0010000000000000;

#endif

constexpr Bits  bitcast(Float x) noexcept { return __builtin_bit_cast(Bits,  x); }
constexpr Float bitcast(Bits  x) noexcept { return __builtin_bit_cast(Float, x); }

constexpr int32_t exponent_from_bits(Bits bits) noexcept {
    return static_cast<int32_t>((bits & ExponentMask) >> (SignificandBits - 1)) - ExponentBias;
}

//==================================================================================================
// floor
//==================================================================================================

constexpr Float ieee754_floor(Float x) noexcept
{
    auto bits = bitcast(x);

    // x = 1.f * 2^e2, where f has at most (p-1) bits.
    //  (The subnormal case is handled as |x| < 1 below.)

    const auto e2 = exponent_from_bits(bits);
    if (e2 >= SignificandBits - 1)
    {
        // |x| >= 2^p

        // x is integral (or Inf or NaN)
        return x;
    }
    else if (e2 >= 0)
    {
        // |x| >= 1

        const auto fraction_mask = SignificandMask >> e2;
        if ((bits & fraction_mask) == 0)
        {
            // x is integral
            return x;
        }

        if ((bits & SignMask) != 0) // (x < 0)
        {
            bits += (fraction_mask + 1); // -1
        }

        return bitcast(bits & ~fraction_mask);
    }
    else
    {
        // |x| < 1

        if ((bits & SignMask) == 0) // (x >= +0)
            return 0;

        if (bits == SignMask)
            return bitcast(SignMask);

        return -1;
    }
}

//==================================================================================================
// ceil
//==================================================================================================

constexpr Float ieee754_ceil(Float x) noexcept
{
    auto bits = bitcast(x);

    // x = 1.f * 2^e2, where f has at most (p-1) bits.
    //  (The subnormal case is handled as |x| < 1 below.)

    const auto e2 = exponent_from_bits(bits);
    if (e2 >= SignificandBits - 1)
    {
        // |x| >= 2^p

        // x is integral (or Inf or NaN)
        return x;
    }
    else if (e2 >= 0)
    {
        // |x| >= 1

        const auto fraction_mask = SignificandMask >> e2;
        if ((bits & fraction_mask) == 0)
        {
            // x is integral
            return x;
        }

        if ((bits & SignMask) == 0) // (x > 0)
        {
            bits += (fraction_mask + 1); // +1
        }

        return bitcast(bits & ~fraction_mask);
    }
    else
    {
        // |x| < 1

        if ((bits & SignMask) != 0) // (x <= -0)
            return bitcast(SignMask);

        if (bits == 0)
            return 0;

        return +1;
    }
}

//==================================================================================================
// trunc
//==================================================================================================

constexpr Float ieee754_trunc(Float x) noexcept
{
    const auto bits = bitcast(x);

    // x = 1.f * 2^e2, where f has at most (p-1) bits.
    //  (The subnormal case is handled as |x| < 1 below.)

    const auto e2 = exponent_from_bits(bits);
    if (e2 >= SignificandBits - 1)
    {
        // |x| >= 2^p

        // x is integral (or Inf or NaN)
        return x;
    }
    else if (e2 >= 0)
    {
        // |x| >= 1

        const auto fraction_mask = SignificandMask >> e2;
        return bitcast(bits & ~fraction_mask);
    }
    else
    {
        // |x| < 1

        return bitcast(bits & SignMask);
    }
}

//==================================================================================================
// round
//==================================================================================================

constexpr Float ieee754_round(Float x) noexcept
{
    auto bits = bitcast(x);

    // x = 1.f * 2^e2, where f has at most (p-1) bits.
    //  (The subnormal case is handled as |x| < 1 below.)

    const auto e2 = exponent_from_bits(bits);
    if (e2 >= SignificandBits - 1)
    {
        // |x| >= 2^p

        // x is integral (or Inf or NaN)
        return x;
    }
    else if (e2 >= 0)
    {
        // |x| >= 1

        const auto fraction_mask = SignificandMask >> e2;
        if ((bits & fraction_mask) == 0)
        {
            // x is integral
            return x;
        }

        bits += (fraction_mask + 1) / 2; // +- 1/2

        return bitcast(bits & ~fraction_mask);
    }
    else
    {
        // |x| < 1

        if (2*x >= +1)
            return +1;
        if (2*x <= -1)
            return -1;

        return bitcast(bits & SignMask);
    }
}

//==================================================================================================
// sqrt
//==================================================================================================

// From:
// https://github.com/freebsd/freebsd-src/blob/de1aa3dab23c06fec962a14da3e7b4755c5880cf/lib/msun/src/e_sqrt.c

/*
 * ====================================================
 * Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
 *
 * Developed at SunSoft, a Sun Microsystems, Inc. business.
 * Permission to use, copy, modify, and distribute this
 * software is freely granted, provided that this notice
 * is preserved.
 * ====================================================
 */

/*
Method:
Bit by bit method using integer arithmetic. (Slow, but portable)

1. Normalization

    Scale x to y in [1,4) with even powers of 2:
    Find an integer k such that 1 <= (y = x * 2^(2k)) < 4, then

        sqrt(x) = 2^k * sqrt(y)

2. Bit by bit computation

    Let q_i = sqrt(y) truncated to i bit after binary point (q_0 = 1),

        s_i = 2 * q_i,
        y_i = 2^(i+1) * (y - q_i^2)                                     (1)

    To compute q_{i+1} from q_i, one checks whether

        (q_i + 2^-(i+1))^2 <= y.                                        (2)

    If (2) is false, then
        q_{i+1} = q_i;
    otherwise
        q_{i+1} = q_i + 2^-(i+1).

    With some algebric manipulation, it is not difficult to see that (2) is
    equivalent to

        s_i + 2^-(i+1) <= y_i                                           (3)

    The advantage of (3) is that s_i and y_i can be computed by the
    following recurrence formula:

        if (3) is false

            s_{i+1} = s_i
            y_{i+1} = y_i;                                              (4)

        otherwise

            s_{i+1} = s_i + 2^-i,
            y_{i+1} = y_i - s_i - 2^-(i+1)                              (5)

    One may easily use induction to prove (4) and (5).

    Note. Since the left hand side of (3) contains only i+2 bits, it does
    not necessary to do a full (53-bit) comparison in (3).

3. Final rounding

    After generating the 53 bits result, we compute one more bit. Together
    with the remainder, we can decide whether the result is exact, bigger
    than 1/2 ulp, or less than 1/2 ulp (it will never equal to 1/2 ulp).

    The rounding mode can be detected by checking whether huge + tiny is
    equal to huge, and whether huge - tiny is equal to huge for some
    floating point number "huge" and "tiny".
*/

constexpr Float ieee754_sqrt(Float x) noexcept
{
    const Bits bits = bitcast(x);

    // sqrt( nan) =  nan
    // sqrt(+inf) = +inf
    // sqrt(-inf) =  snan
    if ((bits & ExponentMask) == ExponentMask)
    {
        return x * x + x;
    }

    // sqrt(+-0) = +-0
    if ((bits & ~SignMask) == 0)
    {
        return x;
    }

    // sqrt(x < 0) = snan
    if ((bits & SignMask) != 0)
    {
        return (x - x) / (x - x);
    }

    // x = m2 * 2^e2,
    // where 2^(p-1) <= m2 < 2^p

    auto m2 = bits & SignificandMask;
    auto e2 = static_cast<int32_t>(bits >> (SignificandBits - 1));
    if (e2 == 0)
    {
        e2 = 1 - ExponentBias;
        while ((m2 & HiddenBit) == 0)
        {
            m2 <<= 1;
            e2  -= 1;
        }
    }
    else
    {
        m2 |= HiddenBit;
        e2 -= ExponentBias;
    }

    // If e2 is odd, make it even.
    if (e2 % 2 != 0)
    {
        m2 <<= 1;
        e2  -= 1;
    }

    m2 <<= 1; // Generate one more bit to correctly round the result

    Bits q = 0;              // q = sqrt(x)
    Bits s = 0;
    Bits r = HiddenBit * 2;  // r = moving bit from right to left
    while (r != 0)
    {
        const auto t = s | r;
        if (t <= m2)
        {
            m2 -= t;
            q |= r;
            s |= r << 1;
        }
        m2 <<= 1;
        r  >>= 1;
    }

    if (m2 != 0)
    {
        // Round to nearest-even.
        q += (q & 1);
    }

    q >>= 1; // Discard rounding bit

    return bitcast(q + (static_cast<Bits>(e2 / 2 - 1 + ExponentBias) << (SignificandBits - 1)));
}

#if 0
constexpr Float kSqrt00 = ieee754_sqrt( 0.0);
constexpr Float kSqrt01 = ieee754_sqrt( 1.0);
constexpr Float kSqrt02 = ieee754_sqrt( 2.0);
constexpr Float kSqrt03 = ieee754_sqrt( 3.0);
constexpr Float kSqrt04 = ieee754_sqrt(-1.0);
constexpr Float kSqrt05 = ieee754_sqrt(-0.0);
constexpr Float kSqrt06 = ieee754_sqrt(-1.0e+300 * 1.0e+300);
constexpr Float kSqrt07 = ieee754_sqrt( 1.0e+300 * 1.0e+300);
#endif

//==================================================================================================
// cbrt
//==================================================================================================

//==================================================================================================
// pow
//==================================================================================================

} // namespace uom::constexpr_math
#endif // ^^^ UNITS_HAS_CONSTEXPR_MATH() ^^^
