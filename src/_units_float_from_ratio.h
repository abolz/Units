// Copyright 2021 Alexander Bolz
//
// Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at https://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <cassert>
#include <cstdint>

#if defined(__SIZEOF_INT128__)
#define UNITS_BUILTIN_UINT128() 1
#else
#define UNITS_BUILTIN_UINT128() 0
#endif

namespace uom::impl {

//==================================================================================================
//
//==================================================================================================

#if UNITS_BUILTIN_UINT128()

__extension__ using uint128_t = unsigned __int128;

constexpr uint128_t ShiftLeft(uint128_t value, int count) noexcept
{
    assert(count > 0);
    assert(count < 64);
    assert(value == (value << count) >> count);
    return value << count;
}

struct DivModResult {
    uint128_t q = 0;
    uint128_t r = 0;
};

constexpr DivModResult DivMod(uint128_t num, uint128_t den) noexcept
{
    const uint128_t q = num / den;
    const uint128_t r = num % den;
    assert(q >= 0);
    assert(q <= UINT64_MAX);
    return {q, r};
}

#else

//
// TODO:
// Fallback-implementation (for MSVC...)
//
// Operations:
//     Compare
//     ShiftLeft (limited)
//     DivMod
//

#endif

//==================================================================================================
//
//==================================================================================================

// Converts num/den to the closest IEEE floating-point number x = f * 2^e.
// PRE: num > 0
// PRE: den > 0
constexpr double F64FromPositiveRatio(int64_t num_, int64_t den_) noexcept
{
    assert(num_ >= 1);
    assert(den_ >= 1);

#if UNITS_BUILTIN_UINT128()

    constexpr int P             = 53;
    constexpr int ExponentBits  = 11;
    constexpr int ExponentBias  = (1 << (ExponentBits - 1)) - 1 + (P - 1);

    constexpr uint64_t TwoP            = uint64_t{1} << P;  // = 2^p
    constexpr uint64_t HiddenBit       = TwoP / 2;          // = 2^(p-1)
    constexpr uint64_t SignificandMask = HiddenBit - 1;     // = 2^(p-1) - 1
    constexpr uint64_t MaxIeeeExponent = (uint64_t{1} << ExponentBits) - 2;

#if 0
    if (num_ <= TwoP && den_ <= TwoP)
        return static_cast<double>(num_) / static_cast<double>(den_);
#endif

    uint128_t num = static_cast<uint64_t>(num_);
    uint128_t den = static_cast<uint64_t>(den_);

    // 1)
    // Scale into [2^(p-1), 2^p).
    int e = 0;
    while (num >= ShiftLeft(den, P)) // num/den >= 2^p
    {
        den = ShiftLeft(den, 1);
        ++e;
    }
    while (num < ShiftLeft(den, P - 1)) // num/den < 2^(p-1)
    {
        num = ShiftLeft(num, 1);
        --e;
    }

    // 2)
    // Divide...
    const auto [q, r] = DivMod(num, den);
    assert(q <= UINT64_MAX);

    uint64_t f = static_cast<uint64_t>(q);
    assert(f >= 0);
    assert(f <= TwoP - 1);

    // 3)
    // Round to nearest-even.
    const bool isExact = (r == 0);
    if (!isExact)
    {
        const uint128_t r2 = ShiftLeft(r, 1);
        if (r2 > den || (r2 == den && (f % 2) != 0))
        {
            ++f;
            if (f == TwoP) // Overflow. Move a trailing zero into the exponent.
            {
                f /= 2;
                ++e;
            }
        }
    }

    assert(f >= 1);
    assert(f <= TwoP - 1);
    assert(e >= 1 - ExponentBias);

    const uint64_t ieee_significand = f & SignificandMask;
    const uint64_t ieee_exponent = static_cast<uint64_t>(e + ExponentBias);
    assert(ieee_exponent <= MaxIeeeExponent);

    return __builtin_bit_cast(double, ieee_significand | ieee_exponent << (P - 1));

#else

    return static_cast<double>(num_) / static_cast<double>(den_);

#endif
}

// Converts num/den to the closest IEEE floating-point number x = f * 2^e.
// PRE: den > 0
constexpr double F64FromRatio(int64_t num, int64_t den) noexcept
{
    assert(num != INT64_MIN);
    assert(den >= 1);

    if (num == 0)
        return 0;

    if (num > 0)
        return  F64FromPositiveRatio( num, den);
    else
        return -F64FromPositiveRatio(-num, den);
}

} // namespace uom::impl
