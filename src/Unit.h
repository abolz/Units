// Copyright 2021 Alexander Bolz
//
// Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at https://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <cassert>
#include <ratio>
#include <type_traits>

#ifndef UNITS_ASSERT
#define UNITS_ASSERT(X) assert(X)
#endif

namespace uom {

//==================================================================================================
//
//==================================================================================================

template <typename R>
inline constexpr bool IsReducedRatio = std::ratio_equal_v<R, typename R::type>;

template <typename T>
inline constexpr bool IsRatio = false;

#if 1
template <int64_t Num, int64_t Den>
inline constexpr bool IsRatio<std::ratio<Num, Den>> = IsReducedRatio<std::ratio<Num, Den>>;
#else
template <int64_t Num, int64_t Den>
inline constexpr bool IsRatio<std::ratio<Num, Den>> = true;
#endif

//==================================================================================================
// Dimension
//==================================================================================================

template <int64_t Num = 1, int64_t Den = 1>
using Dimension = typename std::ratio<Num, Den>::type;

template <typename D1, typename D2>
using MulDimensions = typename std::ratio_multiply<D1, D2>::type;

template <typename D1, typename D2>
using DivDimensions = typename std::ratio_divide<D1, D2>::type;

//==================================================================================================
// Kind
//==================================================================================================

template <typename D, typename Tag>
struct Kind
{
    static_assert(IsRatio<D>,
        "dimension must be a (reduced) std::ratio");

    using type      = Kind;
    using dimension = D;
    using tag       = Tag;
};

template <typename T>
inline constexpr bool IsKind = false;

template <typename D, typename Tag>
inline constexpr bool IsKind<Kind<D, Tag>> = true;

namespace kinds
{
    struct Simple {};

    template <typename ...Factors>
    struct Complex {};

    template <typename Tag, int64_t Exponent>
    struct Factor
    {
        using type = Factor;
        using tag  = Tag;
        static constexpr int64_t exponent = Exponent;
    };

    // Dimensionless [1]
//  using One           = Kind<Dimension<1>, Simple>;
    using Dimensionless = Kind<Dimension<1>, Simple>;
}

namespace kinds::impl
{
    constexpr uint64_t FNV1a(const char* str) noexcept
    {
        uint64_t hash = 14695981039346656037u;
        for ( ; *str != '\0'; ++str)
            hash = (hash ^ static_cast<unsigned char>(*str)) * 1099511628211u;

        return hash;
    }

    template <typename T>
    constexpr uint64_t TypeId() noexcept
    {
#if defined(_MSC_VER) && !defined(__clang__)
        return impl::FNV1a(__FUNCSIG__);
#else
        return impl::FNV1a(__PRETTY_FUNCTION__);
#endif
    }

    template <typename T, typename ...Ts>
    constexpr auto Head(Complex<T, Ts...>) noexcept
    {
        return T{};
    }

    template <typename T, typename ...Ts>
    constexpr auto Tail(Complex<T, Ts...>) noexcept
    {
        return Complex<Ts...>{};
    }

    template <typename ...Ts, typename ...Us>
    constexpr auto Concat(Complex<Ts...>, Complex<Us...>) noexcept
    {
        return Complex<Ts..., Us...>{};
    }

    template <typename ...Fs1, typename ...Fs2>
    constexpr auto Merge([[maybe_unused]] Complex<Fs1...> lhs, [[maybe_unused]] Complex<Fs2...> rhs) noexcept
    {
        if constexpr (sizeof...(Fs1) == 0)
        {
            return rhs;
        }
        else if constexpr (sizeof...(Fs2) == 0)
        {
            return lhs;
        }
        else
        {
            using H1 = decltype(impl::Head(lhs));
            using H2 = decltype(impl::Head(rhs));

            using T1 = typename H1::tag;
            using T2 = typename H2::tag;

            constexpr uint64_t h1 = impl::TypeId<T1>();
            constexpr uint64_t h2 = impl::TypeId<T2>();
            if constexpr (h1 < h2)
            {
                return impl::Concat(Complex<H1>{}, impl::Merge(impl::Tail(lhs), rhs));
            }
            else if constexpr (h1 > h2)
            {
                return impl::Concat(Complex<H2>{}, impl::Merge(lhs, impl::Tail(rhs)));
            }
            else
            {
                static_assert(std::is_same_v<T1, T2>,
                    "collision detected - please try changing the tag name");

                constexpr int64_t e1 = H1::exponent;
                constexpr int64_t e2 = H2::exponent;
                constexpr int64_t e = (e1 + e2 == 0) ? 0 : (std::is_same_v<Simple, T1> ? 1 : (e1 + e2));

                if constexpr (e != 0)
                {
                    using F = Factor<T1, e>;
                    return impl::Concat(Complex<F>{}, impl::Merge(impl::Tail(lhs), impl::Tail(rhs)));
                }
                else
                {
                    return impl::Merge(impl::Tail(lhs), impl::Tail(rhs));
                }
            }
        }
    }

    template <typename Tag, int64_t Exponent>
    constexpr auto Rcp(Factor<Tag, Exponent>) noexcept
    {
        return Factor<Tag, -Exponent>{};
    }

    template <typename ...Fs>
    constexpr auto Rcp(Complex<Fs...>) noexcept
    {
        return Complex<decltype(impl::Rcp(Fs{}))...>{};
    }

    template <typename T>
    struct Wrap { using type = Complex<Factor<T, 1>>; };

    template <typename ...Fs>
    struct Wrap<Complex<Fs...>> { using type = Complex<Fs...>; };

    template <typename T>
    struct Unwrap { using type = T; };

    template <>
    struct Unwrap<Complex<>> { using type = Simple; };

    template <typename T>
    struct Unwrap<Complex<Factor<T, 1>>> { using type = T; };

    template <typename T1, typename T2>
    struct MulTags
    {
        using W1 = typename Wrap<T1>::type;
        using W2 = typename Wrap<T2>::type;
        using WR = decltype( impl::Merge(W1{}, W2{}) );

        using type = typename Unwrap<WR>::type;
    };

    template <typename T1, typename T2>
    struct DivTags
    {
        using W1 = typename Wrap<T1>::type;
        using W2 = typename Wrap<T2>::type;
        using WR = decltype( impl::Merge(W1{}, impl::Rcp(W2{})) );

        using type = typename Unwrap<WR>::type;
    };

#if 1
    // (reduce compile-time???)
    template <>
    struct MulTags<Simple, Simple> { using type = Simple; };

   // (reduce compile-time???)
    template <>
    struct DivTags<Simple, Simple> { using type = Simple; };
#endif

} // namespace kinds::impl

template <typename T1, typename T2>
using MulTags = typename kinds::impl::MulTags<T1, T2>::type;

template <typename T1, typename T2>
using DivTags = typename kinds::impl::DivTags<T1, T2>::type;

template <typename K1, typename K2>
using MulKinds = Kind< MulDimensions<typename K1::dimension, typename K2::dimension>,
                       MulTags<typename K1::tag, typename K2::tag> >;

template <typename K1, typename K2>
using DivKinds = Kind< DivDimensions<typename K1::dimension, typename K2::dimension>,
                       DivTags<typename K1::tag, typename K2::tag> >;

//==================================================================================================
// Conversion
//==================================================================================================

template <int64_t Num, int64_t Den = 1>
using Ratio = typename std::ratio<Num, Den>::type;

// (R::num / R::den) * PI^PiExp
template <typename R, int64_t PiExp = 0>
struct Conversion final
{
    static_assert(IsRatio<R>,
        "R must be a (reduced) std::ratio");

    using type = Conversion;
    using ratio = typename R::type;
    static constexpr int64_t num = ratio::num;
    static constexpr int64_t den = ratio::den;
    static constexpr int64_t exp = PiExp;

    // All integers <= 2^53 are exactly representable as 'double'
    static constexpr int64_t Two53 = 9007199254740992; // == 2^53

    static_assert(num > 0,
        "invalid argument");
    static_assert(num <= Two53,
        "invalid argument");
    static_assert(den > 0,
        "invalid argument");
    static_assert(den <= Two53,
        "invalid argument");
    static_assert(exp >= -4,
        "argument out of range (sorry, not implemented...)");
    static_assert(exp <= 4,
        "argument out of range (sorry, not implemented...)");

    // Returns: (x * num / den) * pi^exp
    [[nodiscard]] constexpr double operator()(double x) const noexcept
    {
        return applyPi(applyRatio(x));
    }

    // Returns: (x * num / den)
    [[nodiscard]] static constexpr double applyRatio(double x) noexcept
    {
//      static_assert(__builtin_is_constant_evaluated());
        static_assert(num >= 1);
        static_assert(den >= 1);
        static_assert(den == 1 || num % den != 0);

        if constexpr (num == 1 && den == 1)
            return x;

        if constexpr (num == 1)
            return x / den;

        if constexpr (den == 1)
            return x * num;

        if constexpr (num >= den)
            return x * (static_cast<double>(num) / static_cast<double>(den));
        else
            return x / (static_cast<double>(den) / static_cast<double>(num));
    }

    // Returns (x * pi^exp)
    [[nodiscard]] static constexpr double applyPi(double x) noexcept
    {
        constexpr double Powers[] = {
            1.0,                // pi^0
            3.141592653589793,  // pi^1
            9.869604401089358,  // pi^2
            31.00627668029982,  // pi^3
            97.40909103400244,  // pi^4
        };

        if constexpr (exp == 0)
            return x;
        else if constexpr (exp > 0)
            return x * Powers[ exp];
        else
            return x / Powers[-exp];
    }
};

template <typename T>
inline constexpr bool IsConversion = false;

template <typename R, int64_t E>
inline constexpr bool IsConversion<Conversion<R, E>> = true; // IsRatio<R>;

template <typename C1, typename C2 /* = C1 */>
using MulConversions = Conversion<typename std::ratio_multiply<typename C1::ratio, typename C2::ratio>::type, C1::exp + C2::exp>;

template <typename C1, typename C2>
using DivConversions = Conversion<typename std::ratio_divide<typename C1::ratio, typename C2::ratio>::type, C1::exp - C2::exp>;

namespace impl
{
#if 0
    template <typename C1>
    using IsIntegralConversion = std::bool_constant<(C1::den == 1 && C1::exp >= 0)>;
#else
    template <typename C1>
    using IsIntegralConversion = std::bool_constant<(C1::den == 1 && C1::exp == 0)>;
#endif

    template <typename C1, typename C2> // (C1 | C2)?
    using ConversionDivides = IsIntegralConversion<DivConversions<C2, C1>>;

} // namespace impl

//==================================================================================================
// Unit
//==================================================================================================

template <typename C, typename K>
struct Unit final
{
    static_assert(IsConversion<C>,
        "C must be a Conversion");
    static_assert(IsKind<K>,
        "K must be a Kind");

    using type       = Unit;
    using conversion = C;
    using kind       = K;
    using dimension  = typename kind::dimension;
    using tag        = typename kind::tag;
};

template <typename T>
inline constexpr bool IsUnit = false;

template <typename C, typename K>
inline constexpr bool IsUnit<Unit<C, K>> = true;

template <typename U1, typename U2>
using MulUnits = typename Unit< MulConversions<typename U1::conversion, typename U2::conversion>,
                                MulKinds<typename U1::kind, typename U2::kind> >::type;

template <typename U1, typename U2>
using DivUnits = typename Unit< DivConversions<typename U1::conversion, typename U2::conversion>,
                                DivKinds<typename U1::kind, typename U2::kind> >::type;

namespace units
{
//  using One               = Unit<Conversion<Ratio<1>>, kinds::One>;
    using Dimensionless     = Unit<Conversion<Ratio<1>>, kinds::Dimensionless>;
}

template <typename U>
inline constexpr bool IsDimensionlessUnit = false;

template <>
inline constexpr bool IsDimensionlessUnit<units::Dimensionless> = true;

//==================================================================================================
// Quantity (value + compile-time unit)
//==================================================================================================

template <typename U>
class Quantity final
{
    // Represents a linear transformation y(x),
    //      y = Cx
    // where C is a fixed scaling factor.

    static_assert(IsUnit<U>,
        "U must be a Unit");

public:
    using type        = Quantity;
    using scalar_type = double;
    using unit        = U;
    using conversion  = typename U::conversion;
    using kind        = typename U::kind;
    using dimension   = typename U::dimension;
    using tag         = typename U::tag;
    using zero        = Ratio<0>;

    using simplified_type
        = Quantity<Unit<conversion, Kind<dimension, kinds::Simple>>>;

private:
    scalar_type _count = 0;

private:
    // asymmetric
    template <typename C1, typename C2, typename T = void>
    using EnableImplicitConversion = std::enable_if_t<impl::ConversionDivides<C1, C2>::value, T>;

    // symmetric
    template <typename K1, typename K2, typename T = void>
    using EnableExplicitConversion
        = std::enable_if_t<std::is_same_v<typename K1::dimension, typename K2::dimension>, T>;

public:
    constexpr Quantity() noexcept = default;
    constexpr Quantity(const Quantity&) noexcept = default;
    constexpr Quantity& operator=(const Quantity&) noexcept = default;

    constexpr explicit Quantity(scalar_type c) noexcept
        : _count(c)
    {
    }

    template <typename C2, EnableImplicitConversion<conversion, C2, int> = 0>
    constexpr Quantity(Quantity<Unit<C2, kind>> q) noexcept
        : _count(DivConversions<C2, conversion>{}(q.count_internal()))
    {
    }

    template <typename C2, typename K2, EnableExplicitConversion<kind, K2, int> = 0>
    constexpr explicit Quantity(Quantity<Unit<C2, K2>> q) noexcept
        : _count(DivConversions<C2, conversion>{}(q.count_internal()))
    {
    }

    [[nodiscard]] constexpr scalar_type count_internal() const noexcept
    {
        return _count;
    }

    [[nodiscard]] constexpr simplified_type value() const noexcept
    {
        return simplified_type(_count);
    }

    template <typename T>
    [[nodiscard]] constexpr bool is_convertible_to() const noexcept;

    template <typename T>
    [[nodiscard]] constexpr T convert_to() const noexcept;

    template <typename T>
    [[nodiscard]] constexpr scalar_type count_as() const noexcept;

    [[nodiscard]] constexpr friend Quantity operator+(Quantity q) noexcept
    {
        return q;
    }

    [[nodiscard]] constexpr friend Quantity operator-(Quantity q) noexcept
    {
        return Quantity(-q.count_internal());
    }

    [[nodiscard]] constexpr friend Quantity operator+(Quantity lhs, Quantity rhs) noexcept
    {
        return Quantity(lhs.count_internal() + rhs.count_internal());
    }

    template <typename _ = unit, std::enable_if_t<IsDimensionlessUnit<_>, int> = 0>
    [[nodiscard]] constexpr friend Quantity operator+(Quantity lhs, scalar_type rhs) noexcept
    {
        return Quantity(lhs.count_internal() + rhs);
    }

    template <typename _ = unit, std::enable_if_t<IsDimensionlessUnit<_>, int> = 0>
    [[nodiscard]] constexpr friend Quantity operator+(scalar_type lhs, Quantity rhs) noexcept
    {
        return Quantity(lhs + rhs.count_internal());
    }

    [[nodiscard]] constexpr friend Quantity operator-(Quantity lhs, Quantity rhs) noexcept
    {
        return Quantity(lhs.count_internal() - rhs.count_internal());
    }

    template <typename _ = unit, std::enable_if_t<IsDimensionlessUnit<_>, int> = 0>
    [[nodiscard]] constexpr friend Quantity operator-(Quantity lhs, scalar_type rhs) noexcept
    {
        return Quantity(lhs.count_internal() - rhs);
    }

    template <typename _ = unit, std::enable_if_t<IsDimensionlessUnit<_>, int> = 0>
    [[nodiscard]] constexpr friend Quantity operator-(scalar_type lhs, Quantity rhs) noexcept
    {
        return Quantity(lhs - rhs.count_internal());
    }

    template <typename U2>
    [[nodiscard]] constexpr friend auto operator*(Quantity lhs, Quantity<U2> rhs) noexcept
    {
        return Quantity<MulUnits<unit, U2>>(lhs.count_internal() * rhs.count_internal());
    }

    template <typename U2>
    [[nodiscard]] constexpr friend auto operator/(Quantity lhs, Quantity<U2> rhs) noexcept
    {
        return Quantity<DivUnits<unit, U2>>(lhs.count_internal() / rhs.count_internal());
    }

    [[nodiscard]] constexpr friend Quantity operator*(Quantity lhs, scalar_type rhs) noexcept
    {
        return Quantity(lhs.count_internal() * rhs);
    }

    [[nodiscard]] constexpr friend Quantity operator/(Quantity lhs, scalar_type rhs) noexcept
    {
        return Quantity(lhs.count_internal() / rhs);
    }

    [[nodiscard]] constexpr friend Quantity operator*(scalar_type lhs, Quantity rhs) noexcept
    {
        return Quantity(lhs * rhs.count_internal());
    }

    [[nodiscard]] constexpr friend auto operator/(scalar_type lhs, Quantity rhs) noexcept
    {
        return Quantity<units::Dimensionless>(lhs) / rhs;
    }

    constexpr friend Quantity& operator+=(Quantity& lhs, Quantity rhs) noexcept
    {
        lhs._count += rhs.count_internal();
        return lhs;
    }

    template <typename _ = unit, std::enable_if_t<IsDimensionlessUnit<_>, int> = 0>
    constexpr friend Quantity& operator+=(Quantity& lhs, scalar_type rhs) noexcept
    {
        lhs._count += rhs;
        return lhs;
    }

    constexpr friend Quantity& operator-=(Quantity& lhs, Quantity rhs) noexcept
    {
        lhs._count -= rhs.count_internal();
        return lhs;
    }

    template <typename _ = unit, std::enable_if_t<IsDimensionlessUnit<_>, int> = 0>
    constexpr friend Quantity& operator-=(Quantity& lhs, scalar_type rhs) noexcept
    {
        lhs._count -= rhs;
        return lhs;
    }

    constexpr friend Quantity& operator*=(Quantity& lhs, scalar_type rhs) noexcept
    {
        lhs._count *= rhs;
        return lhs;
    }

    constexpr friend Quantity& operator/=(Quantity& lhs, scalar_type rhs) noexcept
    {
        lhs._count /= rhs;
        return lhs;
    }

    template <typename _ = unit, std::enable_if_t<IsDimensionlessUnit<_>, int> = 0>
    constexpr friend Quantity& operator*=(Quantity& lhs, Quantity rhs) noexcept
    {
        lhs._count *= rhs.count_internal();
        return lhs;
    }

    template <typename _ = unit, std::enable_if_t<IsDimensionlessUnit<_>, int> = 0>
    constexpr friend Quantity& operator/=(Quantity& lhs, Quantity rhs) noexcept
    {
        lhs._count /= rhs.count_internal();
        return lhs;
    }

    [[nodiscard]] constexpr friend bool operator==(Quantity lhs, Quantity rhs) noexcept
    {
        return lhs.count_internal() == rhs.count_internal();
    }

    [[nodiscard]] constexpr friend bool operator!=(Quantity lhs, Quantity rhs) noexcept
    {
        return lhs.count_internal() != rhs.count_internal();
    }

    [[nodiscard]] constexpr friend bool operator<(Quantity lhs, Quantity rhs) noexcept
    {
        return lhs.count_internal() < rhs.count_internal();
    }

    [[nodiscard]] constexpr friend bool operator>(Quantity lhs, Quantity rhs) noexcept
    {
        return lhs.count_internal() > rhs.count_internal();
    }

    [[nodiscard]] constexpr friend bool operator<=(Quantity lhs, Quantity rhs) noexcept
    {
        return lhs.count_internal() <= rhs.count_internal();
    }

    [[nodiscard]] constexpr friend bool operator>=(Quantity lhs, Quantity rhs) noexcept
    {
        return lhs.count_internal() >= rhs.count_internal();
    }
};

template <typename T>
inline constexpr bool IsQuantity = false;

template <typename U>
inline constexpr bool IsQuantity<Quantity<U>> = true;

template <typename Conv, typename Q>
using ScaledQuantity
    = Quantity<Unit<MulConversions<Conv, typename Q::conversion>, typename Q::kind>>;

template <typename Q, typename Tag>
using TaggedQuantity // a.k.a. Change-Kind
    = Quantity<Unit<typename Q::conversion, Kind<typename Q::dimension, Tag>>>;

//==================================================================================================
// Absolute
//==================================================================================================

template <typename Q, typename Z = Ratio<0>>
class Absolute final
{
    // Represents an affine transformation y(x),
    //      y = C(x + Z) = Cx + Z'
    // where C is a fixed (rational) scaling factor, and Z is a fixed (rational) offset.

    static_assert(IsQuantity<Q>,
        "Q must be a Quantity");

public:
    using type          = Absolute;
    using relative_type = Q;
    using scalar_type   = typename relative_type::scalar_type;
    using unit          = typename relative_type::unit;
    using conversion    = typename relative_type::conversion;
    using kind          = typename relative_type::kind;
    using dimension     = typename relative_type::dimension;
    using tag           = typename relative_type::tag;
    using zero          = Z;

private:
    scalar_type _count;

public:
    constexpr Absolute() noexcept = default;
    constexpr Absolute(const Absolute&) noexcept = default;
    constexpr Absolute& operator=(const Absolute&) noexcept = default;

    constexpr explicit Absolute(scalar_type c) noexcept
        : _count(c)
    {
    }

    [[nodiscard]] constexpr scalar_type count_internal() const noexcept
    {
        return _count;
    }

    template <typename T>
    [[nodiscard]] constexpr bool is_convertible_to() const noexcept;

    template <typename T>
    [[nodiscard]] constexpr T convert_to() const noexcept;

    template <typename T>
    [[nodiscard]] constexpr scalar_type count_as() const noexcept;

    [[nodiscard]] constexpr friend Absolute operator+(Absolute lhs, relative_type rhs) noexcept
    {
        return Absolute(lhs.count_internal() + rhs.count_internal());
    }

    [[nodiscard]] constexpr friend Absolute operator+(relative_type lhs, Absolute rhs) noexcept
    {
        return Absolute(lhs.count_internal() + rhs.count_internal());
    }

    [[nodiscard]] constexpr friend relative_type operator-(Absolute lhs, Absolute rhs) noexcept
    {
        return relative_type(lhs.count_internal() - rhs.count_internal());
    }

    [[nodiscard]] constexpr friend Absolute operator-(Absolute lhs, relative_type rhs) noexcept
    {
        return Absolute(lhs.count_internal() - rhs.count_internal());
    }

    friend Absolute& operator+=(Absolute& lhs, relative_type rhs) noexcept
    {
        lhs._count += rhs.count_internal();
        return lhs;
    }

    friend Absolute& operator-=(Absolute& lhs, relative_type rhs) noexcept
    {
        lhs._count -= rhs.count_internal();
        return lhs;
    }

    [[nodiscard]] constexpr friend bool operator==(Absolute lhs, Absolute rhs) noexcept
    {
        return lhs.count_internal() == rhs.count_internal();
    }

    [[nodiscard]] constexpr friend bool operator!=(Absolute lhs, Absolute rhs) noexcept
    {
        return lhs.count_internal() != rhs.count_internal();
    }

    [[nodiscard]] constexpr friend bool operator<(Absolute lhs, Absolute rhs) noexcept
    {
        return lhs.count_internal() < rhs.count_internal();
    }

    [[nodiscard]] constexpr friend bool operator>(Absolute lhs, Absolute rhs) noexcept
    {
        return lhs.count_internal() > rhs.count_internal();
    }

    [[nodiscard]] constexpr friend bool operator<=(Absolute lhs, Absolute rhs) noexcept
    {
        return lhs.count_internal() <= rhs.count_internal();
    }

    [[nodiscard]] constexpr friend bool operator>=(Absolute lhs, Absolute rhs) noexcept
    {
        return lhs.count_internal() >= rhs.count_internal();
    }
};

template <typename T>
inline constexpr bool IsAbsolute = false;

template <typename Q, typename Z>
inline constexpr bool IsAbsolute<Absolute<Q, Z>> = true;

//==================================================================================================
// convert_to / count_as
//==================================================================================================

namespace impl
{
    template <typename T>
    constexpr double AsDouble() noexcept
    {
        if constexpr (IsRatio<T>)
        {
            static_assert(T::den != 0);
            return static_cast<double>(T::num) / static_cast<double>(T::den);
        }
        else
        {
            // return static_cast<double>(T::value);
            return static_cast<double>(T{}());
        }
    }

    enum class Direction {
        forward,
        backward,
    };

    template <Direction Dir, typename C1, typename Z1, typename C2, typename Z2>
    constexpr double Convert(double x) noexcept
    {
        static_assert(C1::exp == C2::exp,
            "sorry, not supported (yet)");

        // Forward:
        //  convert from C1(x + Z1) to C2(y + Z2)
        //
        //  y = (C1 / C2)( x + Z1 ) - Z2
        //    = (C1 / C2)( x + Z1 - Z2 * (C2 / C1) )
        //    = (C1 / C2)( x ) + (Z1 * (C1 / C2) - Z2)
        //    = a x + b
        //
        // Backward:
        //  convert from C2(y + Z2) to C1(x + Z1)
        //
        //  x = (y - b) / a

        if constexpr (IsRatio<Z1> && IsRatio<Z2>)
        {
            using R1 = std::ratio_divide<typename C1::ratio, typename C2::ratio>;
            using R2 = std::ratio_subtract<std::ratio_multiply<Z1, R1>, Z2>;

            constexpr double a = AsDouble<R1>();
            constexpr double b = AsDouble<R2>();

            if constexpr (Dir == Direction::forward)
                return a * x + b;
            else
                return (x - b) / a;
        }
        else
        {
            // XXX:
            // Use same expression as above?!

            constexpr double a  = AsDouble<std::ratio_divide<typename C1::ratio, typename C2::ratio>>();
            constexpr double z1 = AsDouble<Z1>();
            constexpr double z2 = AsDouble<Z2>();

            if constexpr (Dir == Direction::forward)
                return (x + z1) * a - z2;
            else
                return (x + z2) / a - z1;
        }
    }
}

//--------------------------------------------------------------------------------------------------
//
//--------------------------------------------------------------------------------------------------

template <typename Target, typename SourceUnit>
constexpr Target convert_to(Quantity<SourceUnit> q) noexcept
{
    using Source = Quantity<SourceUnit>;

    static_assert(q.is_convertible_to<Target>(),
        "invalid conversion between unrelated quantities");

    using CS = typename Source::conversion;
    using CT = typename Target::conversion;

    if constexpr (IsQuantity<Target>)
    {
        // (backward)
        return Target(DivConversions<CS, CT>{}(q.count_internal()));
    }
    else
    {
        using ZS = typename Source::zero;
        using ZT = typename Target::zero;
        return Target(impl::Convert<impl::Direction::forward, CS, ZS, CT, ZT>(q.count_internal()));
    }
}

template <typename T, typename U>
constexpr double count_as(Quantity<U> q) noexcept
{
    return uom::convert_to<T>(q).count_internal();
}

template <typename U>
template <typename T>
constexpr bool Quantity<U>::is_convertible_to() const noexcept
{
    using Target = T;
    using Source = Quantity<U>;

    if constexpr (IsQuantity<Target> || IsAbsolute<Target>)
        return std::is_same_v<typename Target::dimension, typename Source::dimension>;
    else
        return false;
}

template <typename U>
template <typename T>
constexpr T Quantity<U>::convert_to() const noexcept
{
    return uom::convert_to<T>(*this);
}

template <typename U>
template <typename T>
constexpr typename Quantity<U>::scalar_type Quantity<U>::count_as() const noexcept
{
    return uom::convert_to<T>(*this).count_internal();
}

//--------------------------------------------------------------------------------------------------
//
//--------------------------------------------------------------------------------------------------

template <typename Target, typename SourceQuantity, typename SourceZero>
constexpr Target convert_to(Absolute<SourceQuantity, SourceZero> a) noexcept
{
    using Source = Absolute<SourceQuantity, SourceZero>;

    static_assert(a.is_convertible_to<Target>(),
        "invalid conversion between unrelated quantities");

    using CS = typename Source::conversion;
    using ZS = typename Source::zero;
    using CT = typename Target::conversion;
    using ZT = typename Target::zero;

    if constexpr (IsQuantity<Target>)
    {
        return Target(impl::Convert<impl::Direction::backward, CT, ZT, CS, ZS>(a.count_internal()));
    }
    else
    {
        return Target(impl::Convert<impl::Direction::forward, CS, ZS, CT, ZT>(a.count_internal()));
    }
}

template <typename T, typename Q, typename Z>
constexpr double count_as(Absolute<Q, Z> q) noexcept
{
    return uom::convert_to<T>(q).count_internal();
}

template <typename Q, typename Z>
template <typename T>
constexpr bool Absolute<Q, Z>::is_convertible_to() const noexcept
{
    using Target = T;
    using Source = Absolute<Q, Z>;

    if constexpr (IsQuantity<Target> || IsAbsolute<Target>)
        return std::is_same_v<typename Target::dimension, typename Source::dimension>;
    else
        return false;
}

template <typename Q, typename Z>
template <typename T>
constexpr T Absolute<Q, Z>::convert_to() const noexcept
{
    return uom::convert_to<T>(*this);
}

template <typename Q, typename Z>
template <typename T>
constexpr typename Absolute<Q, Z>::scalar_type Absolute<Q, Z>::count_as() const noexcept
{
    return uom::convert_to<T>(*this).count_internal();
}

//==================================================================================================
// SI
//==================================================================================================

namespace kinds
{
    // Some prime numbers:
    // 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, ...
    //                         ~~  ^^  ^^  ^^  ^^

    using Length            = Kind< Dimension< 2>, Simple >;
    using Mass              = Kind< Dimension< 3>, Simple >;
    using Time              = Kind< Dimension< 5>, Simple >;
    using ElectricCurrent   = Kind< Dimension< 7>, Simple >;
    using Temperature       = Kind< Dimension<11>, Simple >;
    using AmountOfSubstance = Kind< Dimension<13>, Simple >;
    using LuminousIntensity = Kind< Dimension<17>, Simple >;
    using PlaneAngle        = Kind< Dimension<19>, Simple >;
    using Bit               = Kind< Dimension<23>, Simple >;
    using Entity            = Kind< Dimension<29>, Simple >;
    using Event             = Kind< Dimension<31>, Simple >;
    using Cycle             = Kind< Dimension<37>, Simple >;
}

namespace units
{
    using Meter             = Unit<Conversion<Ratio<1>>, kinds::Length>;
    using Second            = Unit<Conversion<Ratio<1>>, kinds::Time>;
    using Ampere            = Unit<Conversion<Ratio<1>>, kinds::ElectricCurrent>;
    using Kilogram          = Unit<Conversion<Ratio<1>>, kinds::Mass>;
    using Kelvin            = Unit<Conversion<Ratio<1>>, kinds::Temperature>;
    using Mole              = Unit<Conversion<Ratio<1>>, kinds::AmountOfSubstance>;
    using Candela           = Unit<Conversion<Ratio<1>>, kinds::LuminousIntensity>;
    using Radian            = Unit<Conversion<Ratio<1>>, kinds::PlaneAngle>;
    using Bit               = Unit<Conversion<Ratio<1>>, kinds::Bit>;
    using Entity            = Unit<Conversion<Ratio<1>>, kinds::Entity>;
    using Event             = Unit<Conversion<Ratio<1>>, kinds::Event>;
    using Cycle             = Unit<Conversion<Ratio<1>>, kinds::Cycle>;
}

namespace tags
{
    class ApparentPower;
    class AreaPerLength;
    class Radians;
    class ReactivePower;
    class SolidAngle;
}

//------------------------------------------------------------------------------
// One

using Dimensionless     = Quantity<units::Dimensionless>;
using Percent           = ScaledQuantity<Conversion<Ratio<1, 100>>, Dimensionless>;
using Permill           = ScaledQuantity<Conversion<Ratio<1, 1000>>, Dimensionless>;

using Entities          = Quantity<units::Entity>;
using Events            = Quantity<units::Event>;
using Cycles            = Quantity<units::Cycle>;

//------------------------------------------------------------------------------
// Amount of substance

using Moles             = Quantity<units::Mole>;

//------------------------------------------------------------------------------
// Length

using Meters            = Quantity<units::Meter>;
using Millimeters       = ScaledQuantity<Conversion<Ratio<1, 1000>>, Meters>;
using Centimeters       = ScaledQuantity<Conversion<Ratio<1, 100>>, Meters>;
using Decimeters        = ScaledQuantity<Conversion<Ratio<1, 10>>, Meters>;
using Kilometers        = ScaledQuantity<Conversion<Ratio<1000>>, Meters>;

using Inches            = ScaledQuantity<Conversion<Ratio<254, 100>>, Centimeters>; // (international)
using Feet              = ScaledQuantity<Conversion<Ratio<12>>, Inches>;            // (international)
using Yards             = ScaledQuantity<Conversion<Ratio<3>>, Feet>;               // (international)
using Miles             = ScaledQuantity<Conversion<Ratio<1760>>, Yards>;           // (international)

//------------------------------------------------------------------------------
// Area

using SquareMillimeters = decltype(Millimeters{} * Millimeters{});
using SquareCentimeters = decltype(Centimeters{} * Centimeters{});
using SquareDecimeters  = decltype(Decimeters{} * Decimeters{});
using SquareMeters      = decltype(Meters{} * Meters{});
using SquareKilometers  = decltype(Kilometers{} * Kilometers{});

//------------------------------------------------------------------------------
// Volume

using CubicMillimeters  = decltype(SquareMillimeters{} * Millimeters{});
using CubicCentimeters  = decltype(SquareCentimeters{} * Centimeters{});
using CubicDecimeters   = decltype(SquareDecimeters{} * Decimeters{});
using CubicMeters       = decltype(SquareMeters{} * Meters{});

//------------------------------------------------------------------------------
// Plane angle

using Radians           = Quantity<units::Radian>;
using Degrees           = ScaledQuantity<Conversion<Ratio<1, 180>, /* pi^ */ 1>, Radians>;
using Gons              = ScaledQuantity<Conversion<Ratio<1, 200>, /* pi^ */ 1>, Radians>;
using Revolutions       = ScaledQuantity<Conversion<Ratio<2,   1>, /* pi^ */ 1>, Radians>;

using ArcMinutes        = ScaledQuantity<Conversion<Ratio<1, 60>>, Degrees>;
using ArcSeconds        = ScaledQuantity<Conversion<Ratio<1, 60>>, ArcMinutes>;

using ReciprocalRadians = decltype(Dimensionless{} / Radians{});

//------------------------------------------------------------------------------
// Solid angle

#if 1
using Steradians        = decltype(Radians{} * Radians{});
using SquareDegrees     = decltype(Degrees{} * Degrees{});
#else
using Steradians        = TaggedQuantity<decltype(Radians{} * Radians{}), tags::SolidAngle>;
using SquareDegrees     = TaggedQuantity<decltype(Degrees{} * Degrees{}), tags::SolidAngle>;
#endif

//------------------------------------------------------------------------------
// Mass

using Kilograms         = Quantity<units::Kilogram>;
using Grams             = ScaledQuantity<Conversion<Ratio<1, 1000>>, Kilograms>;
using Milligrams        = ScaledQuantity<Conversion<Ratio<1, 1000>>, Grams>;
using Tons              = ScaledQuantity<Conversion<Ratio<1000>>, Kilograms>;

//------------------------------------------------------------------------------
// Density

using KilogramsPerCubicMeter = decltype(Kilograms{} / CubicMeters{});
using TonsPerCubicMeters     = decltype(Tons{}      / CubicMeters{});

using SquareCentimetersPerMeter = TaggedQuantity<decltype(SquareCentimeters{} / Meters{}), tags::AreaPerLength>;
using SquareMetersPerMeter      = TaggedQuantity<decltype(SquareMeters{}      / Meters{}), tags::AreaPerLength>;

//------------------------------------------------------------------------------
// Time

using Seconds           = Quantity<units::Second>;
using Milliseconds      = ScaledQuantity<Conversion<Ratio<1, 1000>>, Seconds>;
using Minutes           = ScaledQuantity<Conversion<Ratio<60>>, Seconds>;
using Hours             = ScaledQuantity<Conversion<Ratio<60>>, Minutes>;
using Days              = ScaledQuantity<Conversion<Ratio<24>>, Hours>;
using Weeks             = ScaledQuantity<Conversion<Ratio<7>>, Days>;
using Years             = ScaledQuantity<Conversion<Ratio<146097, 400>>, Days>;
using Months            = ScaledQuantity<Conversion<Ratio<1, 12>>, Years>;

//------------------------------------------------------------------------------
// Frequency

using Hertz             = decltype(Dimensionless{} / Seconds{});
using Kilohertz         = ScaledQuantity<Conversion<Ratio<1000>>, Hertz>;
using Megahertz         = ScaledQuantity<Conversion<Ratio<1000>>, Kilohertz>;
using Gigahertz         = ScaledQuantity<Conversion<Ratio<1000>>, Megahertz>;

//------------------------------------------------------------------------------
// Speed or velocity

using MetersPerSecond   = decltype(Meters{} / Seconds{});
using KilometersPerHour = decltype(Kilometers{} / Hours{});
using MilesPerHour      = decltype(Miles{} / Hours{});

//------------------------------------------------------------------------------
// Flow (volume)

using CubicMetersPerSecond = decltype(CubicMeters() / Seconds());

//------------------------------------------------------------------------------
// Acceleration

using MetersPerSecondSquared = decltype(Meters{} / (Seconds{} * Seconds{}));

//--------------------------------------------------------------------------
// Force

using Newtons           = decltype(Kilograms{} * Meters{} / (Seconds{} * Seconds{}));
using Kilonewtons       = ScaledQuantity<Conversion<Ratio<1000>>, Newtons>;

//------------------------------------------------------------------------------
// Pressure or mechanical stress

using Pascals           = decltype(Newtons{} / SquareMeters{});

//--------------------------------------------------------------------------
// Torque or moment of force

using NewtonMeters      = decltype(Newtons{} * Meters{} / Radians{});

//--------------------------------------------------------------------------
// Energy

using Joules            = decltype(Newtons{} * Meters{});
using Kilojoules        = ScaledQuantity<Conversion<Ratio<1000>>, Joules>;

//------------------------------------------------------------------------------
// Power or heat flow rate

using Watts             = decltype(Joules{} / Seconds{});
using Kilowatts         = ScaledQuantity<Conversion<Ratio<1000>>, Watts>;

//using Vars              = TaggedQuantity<Watts, tags::ReactivePower>;
//using Kilovars          = ScaledQuantity<Conversion<Ratio<1000>>, Vars>;

//using VoltAmperes       = TaggedQuantity<Watts, tags::ApparentPower>;
//using KiloVoltAmperes   = ScaledQuantity<Conversion<Ratio<1000>>, VoltAmperes>;

//------------------------------------------------------------------------------
// Action

//------------------------------------------------------------------------------
// Dynamic viscosity

using PascalSeconds     = decltype(Pascals{} * Seconds{});

//------------------------------------------------------------------------------
// Kinematic viscosity

using SquareMetersPerSeconds = decltype(SquareMeters{} / Seconds{});

//------------------------------------------------------------------------------
// Electric current

using Amperes           = Quantity<units::Ampere>;

//------------------------------------------------------------------------------
// Electric charge

using Coulombs          = decltype(Amperes{} * Seconds{});

//------------------------------------------------------------------------------
// Electric dipole

//------------------------------------------------------------------------------
// Electromotive force, electric potential difference

using Volts             = decltype(Watts{} / Amperes{});

//------------------------------------------------------------------------------
// Electrical resistance

using Ohms              = decltype(Volts{} / Amperes{});

//------------------------------------------------------------------------------
// Capacitance

using Farads            = decltype(Coulombs{} / Volts{});

//------------------------------------------------------------------------------
// Magnetic flux

using Webers            = decltype(Volts{} * Seconds{});

//------------------------------------------------------------------------------
// Magnetic flux density

using Teslas            = decltype(Webers{} / SquareMeters{});

//------------------------------------------------------------------------------
// Inductance

using Henrys            = decltype(Webers{} / Amperes{});

//------------------------------------------------------------------------------
// Temperature

using Kelvin            = Quantity<units::Kelvin>;
using Rankine           = ScaledQuantity<Conversion<Ratio<5, 9>>, Kelvin>;
using Reaumurs          = ScaledQuantity<Conversion<Ratio<5, 4>>, Kelvin>;

using Millikelvin       = ScaledQuantity<Conversion<Ratio<1, 1000>>, Kelvin>;

//using Celsius           = ScaledQuantity<Conversion<Ratio<1>>, Kelvin>;
//using Fahrenheit        = ScaledQuantity<Conversion<Ratio<5, 9>>, Kelvin>;

using DegKelvin         = Absolute<Kelvin>;
using DegCelsius        = Absolute<Kelvin, Ratio<27315, 100>>;
using DegRankine        = Absolute<Rankine>;
using DegFahrenheit     = Absolute<Rankine, Ratio<45967, 100>>;
using DegReaumur        = Absolute<Reaumurs, Ratio<21852, 100>>;

//--------------------------------------------------------------------------
// Information entropy

using Bits              = Quantity<units::Bit>;
using Nibbles           = ScaledQuantity<Conversion<Ratio<4>>, Bits>;
using Bytes             = ScaledQuantity<Conversion<Ratio<8>>, Bits>;
using Kilobytes         = ScaledQuantity<Conversion<Ratio<1000>>, Bytes>;
using Megabytes         = ScaledQuantity<Conversion<Ratio<1000>>, Kilobytes>;
using Gigabytes         = ScaledQuantity<Conversion<Ratio<1000>>, Megabytes>;

//------------------------------------------------------------------------------
// Luminous intensity

using Candelas          = Quantity<units::Candela>;

//------------------------------------------------------------------------------
// Luminance

using Nits              = decltype(Candelas{} / SquareMeters{});

//------------------------------------------------------------------------------
// Luminous flux

using Lumens            = decltype(Candelas{} * Steradians{});

//------------------------------------------------------------------------------
// Illuminance

using Luxs              = decltype(Lumens{} / SquareMeters{});

//==================================================================================================
//
//==================================================================================================

} // namespace uom
