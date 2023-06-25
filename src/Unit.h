// Copyright 2021 Alexander Bolz
//
// Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at https://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <cassert>
#include <limits>
#include <ratio>
#include <type_traits>

#ifndef UNITS_ASSERT
#define UNITS_ASSERT(X) assert(X)
#endif

namespace uom {

// Internal representation.
// This is NOT an option (yet).
#if 0 // defined(UNITS_SCALAR_F32)
using Scalar = float;
#define UNITS_SCALAR_F32() 1
#define UNITS_SCALAR_F64() 0
#else
using Scalar = double;
#define UNITS_SCALAR_F32() 0
#define UNITS_SCALAR_F64() 1
#endif

//==================================================================================================
//
//==================================================================================================

template <typename T>
inline constexpr bool IsRatio = false;

template <int64_t Num, int64_t Den>
inline constexpr bool IsRatio<std::ratio<Num, Den>> = true;

//==================================================================================================
// Dimension
//==================================================================================================

template <int64_t Num = 1, int64_t Den = 1>
using Dimension = typename std::ratio<Num, Den>::type;

namespace impl
{
    template <typename D1, typename D2>
    using MulDimensions = typename std::ratio_multiply<D1, D2>::type;

    template <typename D1, typename D2>
    using DivDimensions = typename std::ratio_divide<D1, D2>::type;

} // namespace impl

//==================================================================================================
// Kind
//==================================================================================================

template <typename D, typename Tag>
struct Kind
{
    static_assert(IsRatio<D>,
        "dimension must be a (reduced) std::ratio");

    using dimension = typename D::type;
    using tag       = Tag;
    using type      = Kind<dimension, Tag>;
};

template <typename T>
inline constexpr bool IsKind = false;

template <typename D, typename Tag>
inline constexpr bool IsKind<Kind<D, Tag>> = true;

namespace kinds
{
    using Untagged = void;

    template <typename ...Factors>
    struct Product {};

    template <typename T, int64_t Exponent>
    struct Factor
    {
        using type = Factor;
        using tag  = T;
        static constexpr int64_t exponent = Exponent;
    };

    // Dimensionless [1]
    using Dimensionless = Kind<Dimension<1>, Untagged>;
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

    template <typename T1, typename T2>
    constexpr int CompareTags() noexcept
    {
        constexpr uint64_t h1 = impl::TypeId<T1>();
        constexpr uint64_t h2 = impl::TypeId<T2>();

        if constexpr (h1 < h2)
        {
            return -1;
        }
        else if constexpr (h1 > h2)
        {
            return +1;
        }
        else
        {
            static_assert(std::is_same_v<T1, T2>,
                "collision detected - please try changing the tag name");
            return 0;
        }
    }

    template <typename T, typename ...Ts>
    constexpr auto Head(Product<T, Ts...>) noexcept
    {
        return T{};
    }

    template <typename T, typename ...Ts>
    constexpr auto Tail(Product<T, Ts...>) noexcept
    {
        return Product<Ts...>{};
    }

    template <typename ...Ts, typename ...Us>
    constexpr auto Concat(Product<Ts...>, Product<Us...>) noexcept
    {
        return Product<Ts..., Us...>{};
    }

    template <typename ...Fs1, typename ...Fs2>
    constexpr auto Merge([[maybe_unused]] Product<Fs1...> lhs, [[maybe_unused]] Product<Fs2...> rhs) noexcept
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

            constexpr int cmp = impl::CompareTags<T1, T2>();
            if constexpr (0 > cmp)
            {
                return impl::Concat(Product<H1>{}, impl::Merge(impl::Tail(lhs), rhs));
            }
            else if constexpr (cmp > 0)
            {
                return impl::Concat(Product<H2>{}, impl::Merge(lhs, impl::Tail(rhs)));
            }
            else
            {
                static_assert(std::is_same_v<T1, T2>);

                constexpr int64_t e1 = H1::exponent;
                constexpr int64_t e2 = H2::exponent;
                constexpr int64_t e = (e1 + e2 == 0) ? 0 : (std::is_same_v<Untagged, T1> ? 1 : (e1 + e2));

                if constexpr (e != 0)
                {
                    using F = Factor<T1, e>;
                    return impl::Concat(Product<F>{}, impl::Merge(impl::Tail(lhs), impl::Tail(rhs)));
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
    constexpr auto Rcp(Product<Fs...>) noexcept
    {
        return Product<decltype(impl::Rcp(Fs{}))...>{};
    }

    template <typename T>
    struct WrapFactor { using type = Factor<T, 1>; };

    template <typename T, int64_t Exponent>
    struct WrapFactor<Factor<T, Exponent>> { using type = Factor<T, Exponent>; };

    template <typename T>
    struct Wrap { using type = Product<Factor<T, 1>>; };

    template <typename ...Fs>
    struct Wrap<Product<Fs...>> { using type = Product<typename WrapFactor<Fs>::type...>; };

    template <typename T>
    struct Unwrap { using type = T; };

    template <typename T>
    struct Unwrap<Factor<T, 1>> { using type = T; };

    template <>
    struct Unwrap<Product<>> { using type = Untagged; };

    template <typename T>
    struct Unwrap<Product<Factor<T, 1>>> { using type = T; };

    template <typename ...Fs>
    struct Unwrap<Product<Fs...>> { using type = Product<typename Unwrap<Fs>::type...>; };

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
    struct MulTags<Untagged, Untagged> { using type = Untagged; };

    // (reduce compile-time???)
    template <>
    struct DivTags<Untagged, Untagged> { using type = Untagged; };
#endif

} // namespace kinds::impl

namespace impl
{
    template <typename T1, typename T2>
    using MulTags = typename kinds::impl::MulTags<T1, T2>::type;

    template <typename T1, typename T2>
    using DivTags = typename kinds::impl::DivTags<T1, T2>::type;

    template <typename K1, typename K2>
    using MulKinds = typename Kind< typename MulDimensions<typename K1::dimension, typename K2::dimension>::type,
                                    MulTags<typename K1::tag, typename K2::tag> >::type;

    template <typename K1, typename K2>
    using DivKinds = typename Kind< typename DivDimensions<typename K1::dimension, typename K2::dimension>::type,
                                    DivTags<typename K1::tag, typename K2::tag> >::type;

} // namespace impl

//==================================================================================================
//
//==================================================================================================

namespace impl
{
    constexpr Scalar ScalarFromRatio(int64_t num, int64_t den)
    {
        // UNITS_ASSERT(num != INT64_MIN);
        // UNITS_ASSERT(den >= 1);

        UNITS_ASSERT(den != 0);

        return static_cast<Scalar>(static_cast<double>(num) / static_cast<double>(den));
    }

} // namespace impl

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

    using ratio = typename R::type;
    using type = Conversion<ratio, PiExp>;

    // Using 64-bit integers, the conversion factor has a range of approx.
    //  [1.08 10^-19, 9.22 10^18]

    static constexpr int64_t num = ratio::num;
    static constexpr int64_t den = ratio::den;
    static constexpr int64_t exp = PiExp;

#if 0
    // All integers <= 2^24 are exactly representable as 'float'.
    // All integers <= 2^53 are exactly representable as 'double'.
    static constexpr int64_t Max = int64_t{1} << std::numeric_limits<double>::digits;
#else
    static constexpr int64_t Max = INT64_MAX;
#endif

    static_assert(num > 0,
        "invalid argument - numerator must be positive");
    static_assert(num <= Max,
        "invalid argument - numerator too large");
    static_assert(den > 0,
        "invalid argument - denominator must be positive");
    static_assert(den <= Max,
        "invalid argument - denominator too large");

    // Returns: (x * num / den) * pi^exp
    [[nodiscard]] constexpr Scalar operator()(Scalar x) const noexcept
    {
        return applyPi(applyRatio(x));
    }

    // Returns (x * num / den)
    [[nodiscard]] static constexpr Scalar applyRatio(Scalar x) noexcept
    {
        static_assert(num >= 1); // conversion factors are always positive
        static_assert(den >= 1);
        static_assert(den == 1 || num % den != 0);

        if constexpr (num == 1 && den == 1)
        {
            return x;
        }
        else if constexpr (num > den)
        {
            constexpr Scalar scale = impl::ScalarFromRatio(num, den);
            return x * scale;
        }
        else
        {
            static_assert(num != den);

            constexpr Scalar scale = impl::ScalarFromRatio(den, num);
            return x / scale;
        }
    }

    // Returns (x * pi^exp)
    [[nodiscard]] static constexpr Scalar applyPi(Scalar x) noexcept
    {
        static_assert(-4 <= exp && exp <= 4,
            "argument out of range (sorry, not implemented...)");

        constexpr double Powers[] = {
            1.0,                // pi^0
            3.141592653589793,  // pi^1 ==  884279719003555 / 2^48 (IEEE double)
            9.869604401089358,  // pi^2 == 2778046668940015 / 2^48 (IEEE double)
            31.00627668029982,  // pi^3 == 2181872751617887 / 2^46 (IEEE double)
            97.40909103400244,  // pi^4 == 3427277703775251 / 2^45 (IEEE double)
        };

        if constexpr (exp == 0)
            return x;
        else if constexpr (exp > 0)
            return x * static_cast<Scalar>(Powers[ exp]);
        else
            return x / static_cast<Scalar>(Powers[-exp]);
    }
};

template <typename T>
inline constexpr bool IsConversion = false;

template <typename R, int64_t E>
inline constexpr bool IsConversion<Conversion<R, E>> = true; // IsRatio<R>;

namespace impl
{
    template <typename C1, typename C2 /* = C1 */>
    using MulConversions
        = typename Conversion<typename std::ratio_multiply<typename C1::ratio, typename C2::ratio>::type, C1::exp + C2::exp>::type;

    template <typename C1, typename C2>
    using DivConversions
        = typename Conversion<typename std::ratio_divide<typename C1::ratio, typename C2::ratio>::type, C1::exp - C2::exp>::type;

    template <typename C1>
    using IsIntegralConversion = std::bool_constant<(C1::den == 1 && C1::exp == 0)>;

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

    using conversion = typename C::type;
    using kind       = typename K::type;
    using dimension  = typename kind::dimension;
    using tag        = typename kind::tag;
    using type       = Unit<conversion, kind>;
};

template <typename T>
inline constexpr bool IsUnit = false;

template <typename C, typename K>
inline constexpr bool IsUnit<Unit<C, K>> = true;

namespace impl
{
    template <typename U1, typename U2>
    using MulUnits = typename Unit< typename MulConversions<typename U1::conversion, typename U2::conversion>::type,
                                    typename MulKinds<typename U1::kind, typename U2::kind>::type >::type;

    template <typename U1, typename U2>
    using DivUnits = typename Unit< typename DivConversions<typename U1::conversion, typename U2::conversion>::type,
                                    typename DivKinds<typename U1::kind, typename U2::kind>::type >::type;

} // namespace impl

namespace units
{
    using Dimensionless = Unit<Conversion<Ratio<1>>, kinds::Dimensionless>;
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
    using scalar_type = Scalar;
    using unit        = typename U::type;
    using conversion  = typename unit::conversion;
    using kind        = typename unit::kind;
    using dimension   = typename unit::dimension;
    using tag         = typename unit::tag;
    using zero        = Ratio<0>;
    using type        = Quantity<unit>;

    using untagged_type
        = typename Quantity<Unit<conversion, Kind<dimension, kinds::Untagged>>>::type;

    // static constexpr Unit unit = {};

private:
    scalar_type _count = 0;

private:
    // asymmetric
    template <typename UnitTo, typename UnitFrom>
    using IsImplicitlyConvertible
        = std::conjunction<
            // The dimensions must match!
            std::is_same<typename UnitTo::dimension, typename UnitFrom::dimension>,
            // The conversion factor must be an integer!
            impl::ConversionDivides<typename UnitTo::conversion, typename UnitFrom::conversion>,
            // And either
            std::disjunction<
              // the tags must be the same,
              std::is_same<typename UnitFrom::tag, typename UnitTo::tag>,
              // or the right hand side must be a "simple" type (untagged).
              std::is_same<typename UnitFrom::tag, kinds::Untagged>
            >
          >;

    // asymmetric
    template <typename UnitTo, typename UnitFrom>
    using IsExplicitlyConvertible
        = std::conjunction<
            // The dimensions must match!
            std::is_same<typename UnitTo::dimension, typename UnitFrom::dimension>
          >;

    template <typename UnitType, typename T = void>
    using EnableIfDimensionless
        = std::enable_if_t<
            IsDimensionlessUnit<UnitType>,
            T
          >;

    template <typename UnitTo, typename UnitFrom, typename T = void>
    using EnableImplicitConversion
        = std::enable_if_t<
            IsImplicitlyConvertible<UnitTo, UnitFrom>::value,
            T
          >;

    template <typename UnitTo, typename UnitFrom, typename T = void>
    using EnableExplicitConversion
        = std::enable_if_t<
            std::conjunction_v<
              std::negation<IsImplicitlyConvertible<UnitTo, UnitFrom>>,
              IsExplicitlyConvertible<UnitTo, UnitFrom>
            >,
            T
          >;

public:
    constexpr Quantity() noexcept = default;
    constexpr Quantity(const Quantity&) noexcept = default;
    constexpr Quantity& operator=(const Quantity&) noexcept = default;

    constexpr explicit Quantity(scalar_type c) noexcept
        : _count(c)
    {
    }

    template <typename U2, EnableImplicitConversion<unit, U2, int> = 0>
    constexpr Quantity(Quantity<U2> q) noexcept
        : _count(impl::DivConversions<typename U2::conversion, conversion>{}(q._count_internal()))
    {
    }

    template <typename U2, EnableExplicitConversion<unit, U2, int> = 0>
    constexpr explicit Quantity(Quantity<U2> q) noexcept
        : _count(impl::DivConversions<typename U2::conversion, conversion>{}(q._count_internal()))
    {
    }

    template <typename _unit = unit, EnableIfDimensionless<_unit, int> = 0>
    [[nodiscard]] constexpr explicit operator Scalar() const noexcept
    {
        return _count;
    }

    // These are not the droids you're looking for!
    // Use `count_as` to convert this `Quantity` into a scalar.
    [[nodiscard]] constexpr scalar_type _count_internal() const noexcept
    {
        return _count;
    }

    [[nodiscard]] constexpr untagged_type untagged() const noexcept
    {
        return untagged_type(_count);
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
        return Quantity(-q._count_internal());
    }

    [[nodiscard]] constexpr friend Quantity operator+(Quantity lhs, Quantity rhs) noexcept
    {
        return Quantity(lhs._count_internal() + rhs._count_internal());
    }

    [[nodiscard]] constexpr friend Quantity operator-(Quantity lhs, Quantity rhs) noexcept
    {
        return Quantity(lhs._count_internal() - rhs._count_internal());
    }

    template <typename _unit = unit, EnableIfDimensionless<_unit, int> = 0>
    [[nodiscard]] constexpr friend Quantity operator+(Quantity lhs, scalar_type rhs) noexcept
    {
        return Quantity(lhs._count_internal() + rhs);
    }

    template <typename _unit = unit, EnableIfDimensionless<_unit, int> = 0>
    [[nodiscard]] constexpr friend Quantity operator-(Quantity lhs, scalar_type rhs) noexcept
    {
        return Quantity(lhs._count_internal() - rhs);
    }

    template <typename _unit = unit, EnableIfDimensionless<_unit, int> = 0>
    [[nodiscard]] constexpr friend Quantity operator+(scalar_type lhs, Quantity rhs) noexcept
    {
        return Quantity(lhs + rhs._count_internal());
    }

    template <typename _unit = unit, EnableIfDimensionless<_unit, int> = 0>
    [[nodiscard]] constexpr friend Quantity operator-(scalar_type lhs, Quantity rhs) noexcept
    {
        return Quantity(lhs - rhs._count_internal());
    }

    template <typename U2>
    [[nodiscard]] constexpr friend auto operator*(Quantity lhs, Quantity<U2> rhs) noexcept
    {
        using R = impl::MulUnits<unit, U2>;
        return Quantity<Unit<typename R::conversion, typename R::kind>>(lhs._count_internal() * rhs._count_internal());
    }

    template <typename U2>
    [[nodiscard]] constexpr friend auto operator/(Quantity lhs, Quantity<U2> rhs) noexcept
    {
        using R = impl::DivUnits<unit, U2>;
        return Quantity<Unit<typename R::conversion, typename R::kind>>(lhs._count_internal() / rhs._count_internal());
    }

    [[nodiscard]] constexpr friend Quantity operator*(Quantity lhs, scalar_type rhs) noexcept
    {
        return Quantity(lhs._count_internal() * rhs);
    }

    [[nodiscard]] constexpr friend Quantity operator/(Quantity lhs, scalar_type rhs) noexcept
    {
        return Quantity(lhs._count_internal() / rhs);
    }

    [[nodiscard]] constexpr friend Quantity operator*(scalar_type lhs, Quantity rhs) noexcept
    {
        return Quantity(lhs * rhs._count_internal());
    }

    [[nodiscard]] constexpr friend auto operator/(scalar_type lhs, Quantity rhs) noexcept
    {
        using R = impl::DivUnits<units::Dimensionless, unit>;
        return Quantity<Unit<typename R::conversion, typename R::kind>>(lhs / rhs._count_internal());
    }

    constexpr friend Quantity& operator+=(Quantity& lhs, Quantity rhs) noexcept
    {
        lhs._count += rhs._count_internal();
        return lhs;
    }

    constexpr friend Quantity& operator-=(Quantity& lhs, Quantity rhs) noexcept
    {
        lhs._count -= rhs._count_internal();
        return lhs;
    }

    template <typename _unit = unit, EnableIfDimensionless<_unit, int> = 0>
    constexpr friend Quantity& operator+=(Quantity& lhs, scalar_type rhs) noexcept
    {
        lhs._count += rhs;
        return lhs;
    }

    template <typename _unit = unit, EnableIfDimensionless<_unit, int> = 0>
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

    [[nodiscard]] constexpr friend bool operator==(Quantity lhs, Quantity rhs) noexcept
    {
        return lhs._count_internal() == rhs._count_internal();
    }

    [[nodiscard]] constexpr friend bool operator!=(Quantity lhs, Quantity rhs) noexcept
    {
        return lhs._count_internal() != rhs._count_internal();
    }

    [[nodiscard]] constexpr friend bool operator<(Quantity lhs, Quantity rhs) noexcept
    {
        return lhs._count_internal() < rhs._count_internal();
    }

    [[nodiscard]] constexpr friend bool operator>(Quantity lhs, Quantity rhs) noexcept
    {
        return lhs._count_internal() > rhs._count_internal();
    }

    [[nodiscard]] constexpr friend bool operator<=(Quantity lhs, Quantity rhs) noexcept
    {
        return lhs._count_internal() <= rhs._count_internal();
    }

    [[nodiscard]] constexpr friend bool operator>=(Quantity lhs, Quantity rhs) noexcept
    {
        return lhs._count_internal() >= rhs._count_internal();
    }
};

template <typename T>
inline constexpr bool IsQuantity = false;

template <typename U>
inline constexpr bool IsQuantity<Quantity<U>> = true;

template <typename Conv, typename Q>
using ScaledQuantity
    = typename Quantity<
        typename Unit<
          typename impl::MulConversions<Conv, typename Q::conversion>::type,
          typename Q::kind
        >::type
      >::type;

template <typename Q, typename Tag>
using TaggedQuantity // a.k.a. Change-Kind
    = typename Quantity<
        typename Unit<
          typename Q::conversion,
          typename Kind<typename Q::dimension, Tag>::type
        >::type
      >::type;

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
    using relative_type = typename Q::type;
    using scalar_type   = typename relative_type::scalar_type;
    using unit          = typename relative_type::unit;
    using conversion    = typename relative_type::conversion;
    using kind          = typename relative_type::kind;
    using dimension     = typename relative_type::dimension;
    using tag           = typename relative_type::tag;
    using zero          = Z;
    using type          = Absolute<relative_type, Z>;

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

    // These are not the droids you're looking for!
    // Use `count_as` to convert this `Absolute` into a scalar.
    [[nodiscard]] constexpr scalar_type _count_internal() const noexcept
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
        return Absolute(lhs._count_internal() + rhs._count_internal());
    }

    [[nodiscard]] constexpr friend Absolute operator+(relative_type lhs, Absolute rhs) noexcept
    {
        return Absolute(lhs._count_internal() + rhs._count_internal());
    }

    [[nodiscard]] constexpr friend relative_type operator-(Absolute lhs, Absolute rhs) noexcept
    {
        return relative_type(lhs._count_internal() - rhs._count_internal());
    }

    [[nodiscard]] constexpr friend Absolute operator-(Absolute lhs, relative_type rhs) noexcept
    {
        return Absolute(lhs._count_internal() - rhs._count_internal());
    }

    friend Absolute& operator+=(Absolute& lhs, relative_type rhs) noexcept
    {
        lhs._count += rhs._count_internal();
        return lhs;
    }

    friend Absolute& operator-=(Absolute& lhs, relative_type rhs) noexcept
    {
        lhs._count -= rhs._count_internal();
        return lhs;
    }

    [[nodiscard]] constexpr friend bool operator==(Absolute lhs, Absolute rhs) noexcept
    {
        return lhs._count_internal() == rhs._count_internal();
    }

    [[nodiscard]] constexpr friend bool operator!=(Absolute lhs, Absolute rhs) noexcept
    {
        return lhs._count_internal() != rhs._count_internal();
    }

    [[nodiscard]] constexpr friend bool operator<(Absolute lhs, Absolute rhs) noexcept
    {
        return lhs._count_internal() < rhs._count_internal();
    }

    [[nodiscard]] constexpr friend bool operator>(Absolute lhs, Absolute rhs) noexcept
    {
        return lhs._count_internal() > rhs._count_internal();
    }

    [[nodiscard]] constexpr friend bool operator<=(Absolute lhs, Absolute rhs) noexcept
    {
        return lhs._count_internal() <= rhs._count_internal();
    }

    [[nodiscard]] constexpr friend bool operator>=(Absolute lhs, Absolute rhs) noexcept
    {
        return lhs._count_internal() >= rhs._count_internal();
    }
};

template <typename T>
inline constexpr bool IsAbsolute = false;

template <typename Q, typename Z>
inline constexpr bool IsAbsolute<Absolute<Q, Z>> = true;

//==================================================================================================
// convert_to / count
//==================================================================================================

namespace impl
{
    template <typename T>
    constexpr Scalar AsScalar() noexcept
    {
        // static_assert(IsRatio<T>,
        //     "T must be a std::ratio");

        if constexpr (IsRatio<T>)
        {
            constexpr Scalar value = impl::ScalarFromRatio(T::num, T::den);
            return value;
        }
        // else
        // {
        //     constexpr Scalar value = static_cast<Scalar>(T::value);
        //     return value;
        // }
        // else
        // {
        //     constexpr Scalar value = static_cast<Scalar>(T{});
        //     return value;
        // }
    }

    enum class Direction {
        forward,
        backward,
    };

    template <Direction Dir, typename C1, typename Z1, typename C2, typename Z2>
    constexpr Scalar Convert(Scalar x) noexcept
    {
        static_assert(IsConversion<C1>,
            "C1 must be a conversion");
        static_assert(IsConversion<C2>,
            "C2 must be a conversion");
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

            constexpr Scalar a = ScalarFromRatio(R1::num, R1::den);
            constexpr Scalar b = ScalarFromRatio(R2::num, R2::den);

            if constexpr (Dir == Direction::forward)
                return a * x + b;
            else
                return (x - b) / a;
        }
        else
        {
            // XXX:
            // Use same expression as above?!

            using CR = std::ratio_divide<typename C1::ratio, typename C2::ratio>;

            constexpr Scalar a  = impl::AsScalar<CR>();
            constexpr Scalar z1 = impl::AsScalar<Z1>();
            constexpr Scalar z2 = impl::AsScalar<Z2>();

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
        return Target(impl::DivConversions<CS, CT>{}(q._count_internal()));
    }
    else
    {
        using ZS = typename Source::zero;
        using ZT = typename Target::zero;
        return Target(impl::Convert<impl::Direction::forward, CS, ZS, CT, ZT>(q._count_internal()));
    }
}

template <typename T, typename U>
constexpr Scalar count_as(Quantity<U> q) noexcept
{
    return uom::convert_to<T>(q)._count_internal();
}

template <typename U>
template <typename T>
constexpr bool Quantity<U>::is_convertible_to() const noexcept
{
    using Target = T;
    using Source = Quantity<U>;

    //
    // FIXME:
    // This is too lenient...
    //
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
    return uom::convert_to<T>(*this)._count_internal();
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
        return Target(impl::Convert<impl::Direction::backward, CT, ZT, CS, ZS>(a._count_internal()));
    }
    else
    {
        return Target(impl::Convert<impl::Direction::forward, CS, ZS, CT, ZT>(a._count_internal()));
    }
}

template <typename T, typename Q, typename Z>
constexpr Scalar count_as(Absolute<Q, Z> q) noexcept
{
    return uom::convert_to<T>(q)._count_internal();
}

template <typename Q, typename Z>
template <typename T>
constexpr bool Absolute<Q, Z>::is_convertible_to() const noexcept
{
    using Target = T;
    using Source = Absolute<Q, Z>;

    //
    // FIXME:
    // This is too lenient...
    //
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
    return uom::convert_to<T>(*this)._count_internal();
}

//==================================================================================================
// SI
//==================================================================================================

namespace tags
{
    // struct Length;
    // struct Mass;
    // struct Time;
    // struct ElectricCurrent;
    // struct Temperature;
    // struct AmountOfSubstance;
    // struct LuminousIntensity;
    // struct PlaneAngle;
    // struct Bit;
    // struct Entity;
    // struct Event;
    // struct Cycle;

    struct ApparentPower;
    struct AreaPerLength;
    struct Radians;
    struct ReactivePower;
    struct SolidAngle;
}

namespace kinds
{
    // Some prime numbers:
    // 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, ...
    //                         ~~  ^^  ^^  ^^  ^^

    using Length            = Kind< Dimension<   2,    1>, Untagged >;
    using Mass              = Kind< Dimension<   3,    1>, Untagged >;
    using Time              = Kind< Dimension<   5,    1>, Untagged >;
    using ElectricCurrent   = Kind< Dimension<   7,    1>, Untagged >;
    using Temperature       = Kind< Dimension<  11,    1>, Untagged >;
    using AmountOfSubstance = Kind< Dimension<  13,    1>, Untagged >;
    using LuminousIntensity = Kind< Dimension<  17,    1>, Untagged >;
    using PlaneAngle        = Kind< Dimension<  19,    1>, Untagged >; // ~~
    using Bit               = Kind< Dimension<  23,    1>, Untagged >; // ^^
    using Entity            = Kind< Dimension<  29,    1>, Untagged >; // ^^
    using Event             = Kind< Dimension<  31,    1>, Untagged >; // ^^
    using Cycle             = Kind< Dimension<  37,    1>, Untagged >; // ^^

#if 1 // #ifdef __INTELLISENSE__
    // EXPERIMENT:
    // These typedef's are only used to generate slightly prettier error messages and IntelliSense hints.
    // They are otherwise unused and useless.

    using Area              = Kind< Dimension<   4,    1>, Untagged >;
    using Volume            = Kind< Dimension<   6,    1>, Untagged >;
    using SolidAngle        = Kind< Dimension< 361,    1>, Untagged >;
    using Density           = Kind< Dimension<   3,    8>, Untagged >;
    using Frequency         = Kind< Dimension<   1,    5>, Untagged >;
    using Velocity          = Kind< Dimension<   2,    5>, Untagged >;
    using Acceleration      = Kind< Dimension<   2,   25>, Untagged >;
    using Force             = Kind< Dimension<   6,   25>, Untagged >;
    using Pressure          = Kind< Dimension<   6,  100>, Untagged >;
    using Energy            = Kind< Dimension<  12,   25>, Untagged >;
    using Power             = Kind< Dimension<  12,  125>, Untagged >;
    using ElectricCharge    = Kind< Dimension<  35,    1>, Untagged >;
    using Luminance         = Kind< Dimension<  17,    4>, Untagged >;
    using LuminousFlux      = Kind< Dimension<6137,    1>, Untagged >;
    using Illuminance       = Kind< Dimension<6137,    4>, Untagged >;
#endif
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

//------------------------------------------------------------------------------
// Dimensionless

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

using AstronomicalUnits = ScaledQuantity<Conversion<Ratio<149597870700>>, Meters>;
using Parsecs           = ScaledQuantity<Conversion<Ratio<648000>, /*pi^*/ -1>, AstronomicalUnits>;

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

using KilogramsPerCubicMeter
                        = decltype(Kilograms{} / CubicMeters{});
using TonsPerCubicMeters
                        = decltype(Tons{}      / CubicMeters{});

using SquareCentimetersPerMeter
                        = TaggedQuantity<decltype(SquareCentimeters{} / Meters{}), tags::AreaPerLength>;
using SquareMetersPerMeter
                        = TaggedQuantity<decltype(SquareMeters{}      / Meters{}), tags::AreaPerLength>;

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

using CubicMetersPerSecond
                        = decltype(CubicMeters() / Seconds());

//------------------------------------------------------------------------------
// Acceleration

using MetersPerSecondSquared
                        = decltype(Meters{} / (Seconds{} * Seconds{}));

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

using SquareMetersPerSeconds
                        = decltype(SquareMeters{} / Seconds{});

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
