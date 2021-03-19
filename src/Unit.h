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
// Gcd/Lcm
//==================================================================================================

namespace impl
{
    // <numeric> std::gcd ?
    constexpr int64_t Gcd(int64_t x, int64_t y) noexcept
    {
        UNITS_ASSERT(x >= 0);
        UNITS_ASSERT(y >= 1);

        while (y > 0) {
            const auto r = x % y;
            x = y;
            y = r;
        }
        return x;
    }

    // <numeric> std::lcm ?
    constexpr int64_t Lcm(int64_t x, int64_t y) noexcept
    {
        return (x / Gcd(x, y)) * y;
    }

} // namespace impl

//==================================================================================================
// Kind
//==================================================================================================

template <typename Tag, typename D>
struct Kind
{
    using type      = Kind;
    using tag       = Tag;
    using dimension = D;
};

namespace kinds
{
    template <int64_t Num = 1, int64_t Den = 1>
    using Dimension = typename std::ratio<Num, Den>::type;

    template <typename D1, typename D2>
    using MulDimensions = typename std::ratio_multiply<D1, D2>::type;

    template <typename D1, typename D2>
    using DivDimensions = typename std::ratio_divide<D1, D2>::type;

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

    namespace impl
    {
        constexpr uint64_t FNV1a(const char* str) noexcept
        {
            uint64_t hash = 14695981039346656037u;
            for ( ; *str != '\0'; ++str)
            {
                hash = (hash ^ *str) * 1099511628211u;
            }

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
        constexpr T Head(Complex<T, Ts...>) noexcept
        {
            return {};
        }

        template <typename T, typename ...Ts>
        constexpr Complex<Ts...> Tail(Complex<T, Ts...>) noexcept
        {
            return {};
        }

        template <typename ...Ts, typename ...Us>
        constexpr Complex<Ts..., Us...> operator+(Complex<Ts...>, Complex<Us...>) noexcept
        {
            return {};
        }

        enum class MergeType
        {
            mul,
            div,
        };

        template <MergeType MT, typename ...Fs1, typename ...Fs2>
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
                    return Complex<H1>{} + impl::Merge<MT>(impl::Tail(lhs), rhs);
                }
                else if constexpr (h1 > h2)
                {
                    return Complex<H2>{} + impl::Merge<MT>(lhs, impl::Tail(rhs));
                }
                else
                {
                    static_assert(std::is_same_v<T1, T2>,
                        "collision detected - please try changing the tag name");

                    constexpr int64_t e1 = H1::exponent;
                    constexpr int64_t e2 = MT == MergeType::mul ? H2::exponent : -H2::exponent;

                    constexpr int64_t e = (e1 + e2 == 0) ? 0 : (std::is_same_v<Simple, T1> ? 1 : (e1 + e2));
                    if constexpr (e != 0)
                    {
                        using F = Factor<T1, e>;
                        return Complex<F>{} + impl::Merge<MT>(impl::Tail(lhs), impl::Tail(rhs));
                    }
                    else
                    {
                        return impl::Merge<MT>(impl::Tail(lhs), impl::Tail(rhs));
                    }
                }
            }
        }

        template <typename T>
        struct Wrap1 { using type = Complex<Factor<T, 1>>; };

        template <typename ...Fs>
        struct Wrap1<Complex<Fs...>> { using type = Complex<Fs...>; };

        template <typename T>
        struct Unwrap1 { using type = T; };

        template <>
        struct Unwrap1<Complex<>> { using type = Simple; };

        template <typename T>
        struct Unwrap1<Complex<Factor<T, 1>>> { using type = T; };

    } // namespace impl

    template <typename T>
    using Wrap = typename impl::Wrap1<T>::type;

    template <typename T>
    using Unwrap = typename impl::Unwrap1<T>::type;

    template <typename T1, typename T2>
    struct MulTags1
    {
        using type = Unwrap< decltype( impl::Merge<impl::MergeType::mul>(Wrap<T1>{}, Wrap<T2>{}) ) >;
    };

    template <typename T1, typename T2>
    struct DivTags1
    {
        using type = Unwrap< decltype( impl::Merge<impl::MergeType::div>(Wrap<T1>{}, Wrap<T2>{}) ) >;
    };

#if 1
    // (reduce compile-time???)
    template <>
    struct MulTags1<Simple, Simple> { using type = Simple; };

    // (reduce compile-time???)
    template <>
    struct DivTags1<Simple, Simple> { using type = Simple; };
#endif

    template <typename T1, typename T2>
    using MulTags = typename kinds::MulTags1<T1, T2>::type;

    template <typename T1, typename T2>
    using DivTags = typename kinds::DivTags1<T1, T2>::type;

} // namespace kinds

template <typename K1, typename K2>
using MulKinds = Kind< kinds::MulTags<typename K1::tag, typename K2::tag>,
                       kinds::MulDimensions<typename K1::dimension, typename K2::dimension> >;

template <typename K1, typename K2>
using DivKinds = Kind< kinds::DivTags<typename K1::tag, typename K2::tag>,
                       kinds::DivDimensions<typename K1::dimension, typename K2::dimension> >;

namespace kinds
{
    //----------------------------------------------------------------------
    // Base kinds

    // Some prime numbers:
    // 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, ...
    //                             ^^  ^^  ^^

    namespace dim
    {
        using One               = Dimension< 1>; // 1
        using Length            = Dimension< 2>; // Meter m
        using Mass              = Dimension< 3>; // Kilogram kg
        using Time              = Dimension< 5>; // Second s
        using ElectricCurrent   = Dimension< 7>; // Ampere A
        using Temperature       = Dimension<11>; // Kelvin K
        using AmountOfSubstance = Dimension<13>; // Mole mol
        using LuminousIntensity = Dimension<17>; // Candela cd
        using PlaneAngle        = Dimension<19>; // Radian rad

//      using Entity            = Dimension<23>;
//      using Event             = Dimension<29>;
        using Bit               = Dimension<31>;

    } // namespace dim

    using One               = Kind< Simple, dim::One               >; // 1
    using Length            = Kind< Simple, dim::Length            >; // Meter m
    using Mass              = Kind< Simple, dim::Mass              >; // Kilogram kg
    using Time              = Kind< Simple, dim::Time              >; // Second s
    using ElectricCurrent   = Kind< Simple, dim::ElectricCurrent   >; // Ampere A
    using Temperature       = Kind< Simple, dim::Temperature       >; // Kelvin K
    using AmountOfSubstance = Kind< Simple, dim::AmountOfSubstance >; // Mole mol
    using LuminousIntensity = Kind< Simple, dim::LuminousIntensity >; // Candela cd
    using PlaneAngle        = Kind< Simple, dim::PlaneAngle        >; // Radian rad

    using Bit               = Kind< Simple, dim::Bit >;

    //----------------------------------------------------------------------
    // Derived kinds

    // m^2/m
    struct AreaPerLength {};

    // m^3/m^2
    struct VolumePerArea {};

    // Steradian sr = rad^2
    struct SolidAngle {};

} // namespace kinds

//==================================================================================================
// Conversion
//==================================================================================================

template <int64_t Num, int64_t Den = 1>
using Ratio = typename std::ratio<Num, Den>::type;

// (R::num / R::den) * PI^PiExp
template <typename R, int64_t PiExp = 0>
struct Conversion final
{
    using ratio = typename R::type;
    static constexpr int64_t num = ratio::num;
    static constexpr int64_t den = ratio::den;
    static constexpr int64_t exp = PiExp;

    // All integers <= 2^53 are exactly representable as 'double'
    static constexpr int64_t Two53 = 9007199254740992; // == 2^53

    static_assert(num >= -Two53,
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

    // Returns x * num / den
    [[nodiscard]] static constexpr double applyRatio(double x) noexcept
    {
        if constexpr (den == 1)
        {
            if constexpr (num == 1)
                return x;
            else if constexpr (num == -1)
                return -x;
            else
                return x * num;
        }
        else
        {
            if constexpr (num == 1)
                return x / den;
            else if constexpr (num == -1)
                return -x / den;
            else
                return x * (static_cast<double>(num) / static_cast<double>(den));
        }
    }

    // Returns x * pi^exp
    [[nodiscard]] static constexpr double applyPi(double x) noexcept
    {
        constexpr double Powers[] = {
            0.010265982254684336, // pi^-4
            0.03225153443319949,  // pi^-3
            0.10132118364233778,  // pi^-2
            0.3183098861837907,   // pi^-1
            1,                    // pi^ 0
            3.141592653589793,    // pi^ 1
            9.869604401089358,    // pi^ 2
            31.00627668029982,    // pi^ 3
            97.40909103400244,    // pi^ 4
        };

        if constexpr (exp == 0)
            return x;
        else
            return x * Powers[exp + 4];
    }
};

template <typename C1, typename C2 /* = C1 */>
using MulConversions = Conversion<std::ratio_multiply<typename C1::ratio, typename C2::ratio>, C1::exp + C2::exp>;

template <typename C1, typename C2>
using DivConversions = Conversion<std::ratio_divide<typename C1::ratio, typename C2::ratio>, C1::exp - C2::exp>;

namespace impl
{
    template <typename C1>
    using IsIntegralConversion = std::bool_constant<(C1::den == 1 && C1::exp == 0)>;

    template <typename C1, typename C2> // (C1 | C2)?
    using ConversionDivides = IsIntegralConversion<DivConversions<C2, C1>>;

    template <typename C1, typename C2>
    struct CommonConversionImpl
    {
        // SFINAE:
        // missing 'type' !!!
    };

    template <typename R1, typename R2, int64_t CommonPiExp>
    struct CommonConversionImpl<Conversion<R1, CommonPiExp>, Conversion<R2, CommonPiExp>>
    {
        // = gcd(N1/D1, N2/D2)
        //   = gcd(N1, N2) / lcm(D1, D2)
        //
        static constexpr int64_t num = impl::Gcd(R1::num, R2::num);
        static constexpr int64_t den = impl::Lcm(R1::den, R2::den);

        using type = Conversion<Ratio<num, den>, CommonPiExp>;
    };

    template <typename C1, typename C2>
    using CommonConversion = typename CommonConversionImpl<C1, C2>::type;

} // namespace impl

//==================================================================================================
// Unit
//==================================================================================================

template <typename C, typename K>
struct Unit final
{
    using type       = Unit;
    using conversion = C;
    using kind       = K;
    using dimension  = typename kind::dimension;
};

template <typename U1, typename U2>
using MulUnits = typename Unit< MulConversions<typename U1::conversion, typename U2::conversion>,
                                MulKinds<typename U1::kind, typename U2::kind> >::type;

template <typename U1, typename U2>
using DivUnits = typename Unit< DivConversions<typename U1::conversion, typename U2::conversion>,
                                DivKinds<typename U1::kind, typename U2::kind> >::type;

namespace impl
{
    template <typename T>
    inline constexpr bool IsUnit = false;

    template <typename C, typename K>
    inline constexpr bool IsUnit<Unit<C, K>> = true;
    //
    // TODO:
    //  IsConversion<C>
    //  IsKind<K>
    //

} // namespace impl

//--------------------------------------------------------------------------------------------------
// Typedefs
//--------------------------------------------------------------------------------------------------

namespace units
{
    //--------------------------------------------------------------------------

    using One               = Unit<Conversion<Ratio<1>>, kinds::One>;
    using Dimensionless     = Unit<Conversion<Ratio<1>>, kinds::One>;

    //--------------------------------------------------------------------------
    // SI base units

    using Meter             = Unit<Conversion<Ratio<1>>, kinds::Length>;
    using Second            = Unit<Conversion<Ratio<1>>, kinds::Time>;
    using Kilogram          = Unit<Conversion<Ratio<1>>, kinds::Mass>;
    using Kelvin            = Unit<Conversion<Ratio<1>>, kinds::Temperature>;
    using Mole              = Unit<Conversion<Ratio<1>>, kinds::AmountOfSubstance>;
    using Candela           = Unit<Conversion<Ratio<1>>, kinds::LuminousIntensity>;
    using Radian            = Unit<Conversion<Ratio<1>>, kinds::PlaneAngle>;

    //--------------------------------------------------------------------------

    using Bit               = Unit<Conversion<Ratio<1>>, kinds::Bit>;

} // namespace units

//==================================================================================================
// Quantity (value + compile-time unit)
//==================================================================================================

template <typename U>
class Quantity;

template <typename OldTag, typename NewTag>
struct Retag
{
};

namespace impl
{
    template <typename T>
    inline constexpr bool IsQuantity = false;

    template <typename U>
    inline constexpr bool IsQuantity<Quantity<U>> = true;

    template <typename T>
    struct Untag1
    {
        using type = T;
    };

    template <typename OldTag, typename NewTag>
    struct Untag1<Retag<OldTag, NewTag>>
    {
        using type = OldTag;
    };

    template <typename T>
    using Untag = typename Untag1<T>::type;

} // namespace impl

template <typename Q, typename Tag>
using Tagged // a.k.a. Change-Kind
    = Quantity<Unit<typename Q::conversion, Kind<Retag<typename Q::tag, Tag>, typename Q::dimension>>>;

template <typename Q>
using Untagged // a.k.a. Change-Kind
    = Quantity<Unit<typename Q::conversion, Kind<impl::Untag<typename Q::tag>, typename Q::dimension>>>;

//--------------------------------------------------------------------------------------------------
//
//--------------------------------------------------------------------------------------------------

template <typename U>
class Quantity final
{
    static_assert(impl::IsUnit<U>,
        "'Quantity' can only be used with 'Unit's");

public:
    using type        = Quantity;
    using scalar_type = double;
    using unit        = U;
    using conversion  = typename U::conversion;
    using kind        = typename U::kind;
    using dimension   = typename U::dimension;
    using tag         = typename kind::tag;

    using untagged_type = Untagged<Quantity>;
    using simplified_type = Tagged<Quantity, kinds::Simple>;

private:
    scalar_type _count = 0;

private:
    template <typename C2>
    using CommonQuantity = Quantity<Unit<impl::CommonConversion<conversion, C2>, kind>>;

    // asymmetric
    template <typename C1, typename C2, typename T = void>
    using EnableImplicitConversion = std::enable_if_t<impl::ConversionDivides<C1, C2>::value, T>;

    // symmetric
    template <typename K1, typename K2, typename T = void>
    using EnableExplicitConversion = std::enable_if_t<std::is_same_v<typename K1::dimension, typename K2::dimension>, T>; // (ratio_equal?)

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

    [[nodiscard]] constexpr untagged_type untag() const noexcept
    {
        return untagged_type(_count);
    }

    [[nodiscard]] constexpr simplified_type simplify() const noexcept
    {
        return simplified_type(_count);
    }

    [[nodiscard]] constexpr double count_internal() const noexcept
    {
        return _count;
    }

    template <typename Q, EnableExplicitConversion<typename Q::kind, kind, int> = 0>
    [[nodiscard]] constexpr double count() const noexcept
    {
        if constexpr (std::is_same_v<Q, Quantity>)
            return count_internal();
        else
            return Q(*this).count_internal();
    }

    //--------------------------------------------------------------------------
    // Arithmetic

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

    [[nodiscard]] constexpr friend Quantity operator-(Quantity lhs, Quantity rhs) noexcept
    {
        return Quantity(lhs.count_internal() - rhs.count_internal());
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
        return Quantity<DivUnits<units::One, unit>>(lhs / rhs.count_internal());
    }

    constexpr friend Quantity& operator+=(Quantity& lhs, Quantity rhs) noexcept
    {
        lhs._count += rhs.count_internal();
        return lhs;
    }

    constexpr friend Quantity& operator-=(Quantity& lhs, Quantity rhs) noexcept
    {
        lhs._count -= rhs.count_internal();
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

    //--------------------------------------------------------------------------
    // Comparisons

    template <typename C2, typename Q = CommonQuantity<C2>>
    [[nodiscard]] constexpr friend int compare(Quantity lhs, Quantity<Unit<C2, kind>> rhs) noexcept
    {
        const auto x = Q(lhs).count_internal();
        const auto y = Q(rhs).count_internal();
        if (x < y)
            return -1;
        if (x > y)
            return +1;
        return 0;
    }

    template <typename C2, typename Q = CommonQuantity<C2>>
    [[nodiscard]] constexpr friend bool operator==(Quantity lhs, Quantity<Unit<C2, kind>> rhs) noexcept
    {
        return Q(lhs).count_internal() == Q(rhs).count_internal();
    }

    template <typename C2, typename Q = CommonQuantity<C2>>
    [[nodiscard]] constexpr friend bool operator!=(Quantity lhs, Quantity<Unit<C2, kind>> rhs) noexcept
    {
        return Q(lhs).count_internal() != Q(rhs).count_internal();
    }

    template <typename C2, typename Q = CommonQuantity<C2>>
    [[nodiscard]] constexpr friend bool operator<(Quantity lhs, Quantity<Unit<C2, kind>> rhs) noexcept
    {
        return Q(lhs).count_internal() < Q(rhs).count_internal();
    }

    template <typename C2, typename Q = CommonQuantity<C2>>
    [[nodiscard]] constexpr friend bool operator>(Quantity lhs, Quantity<Unit<C2, kind>> rhs) noexcept
    {
        return Q(lhs).count_internal() > Q(rhs).count_internal();
    }

    template <typename C2, typename Q = CommonQuantity<C2>>
    [[nodiscard]] constexpr friend bool operator<=(Quantity lhs, Quantity<Unit<C2, kind>> rhs) noexcept
    {
        return Q(lhs).count_internal() <= Q(rhs).count_internal();
    }

    template <typename C2, typename Q = CommonQuantity<C2>>
    [[nodiscard]] constexpr friend bool operator>=(Quantity lhs, Quantity<Unit<C2, kind>> rhs) noexcept
    {
        return Q(lhs).count_internal() >= Q(rhs).count_internal();
    }
};

template <typename U>
[[nodiscard]] constexpr auto simplify(Quantity<U> q) noexcept -> decltype(q.simplify())
{
    return q.simplify();
}

template <typename Q, typename U>
[[nodiscard]] constexpr auto count(Quantity<U> q) noexcept -> decltype(q.template count<Q>())
{
    return q.template count<Q>();
}

//--------------------------------------------------------------------------------------------------
// Typedefs
//--------------------------------------------------------------------------------------------------

template <typename Conv, typename Q>
using ScaledQuantity
    = Quantity<Unit<MulConversions<Conv, typename Q::conversion>, typename Q::kind>>;

//------------------------------------------------------------------------------
// One

using Dimensionless     = Quantity<units::Dimensionless>;
using Percent           = ScaledQuantity<Conversion<Ratio<1, 100>>, Dimensionless>;
using Permill           = ScaledQuantity<Conversion<Ratio<1, 1000>>, Dimensionless>;

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
// Mass

using Kilograms         = Quantity<units::Kilogram>;
using Grams             = ScaledQuantity<Conversion<Ratio<1, 1000>>, Kilograms>;
using Milligrams        = ScaledQuantity<Conversion<Ratio<1, 1000>>, Grams>;
using Tons              = ScaledQuantity<Conversion<Ratio<1000>>, Kilograms>;

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
// Temperature

using Kelvin            = Quantity<units::Kelvin>;

//------------------------------------------------------------------------------
// Amount of substance

using Moles             = Quantity<units::Mole>;

//------------------------------------------------------------------------------
// Luminous intensity

using Candelas          = Quantity<units::Candela>;

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
// Area

using SquareMillimeters = decltype(Millimeters{} * Millimeters{});
using SquareCentimeters = decltype(Centimeters{} * Centimeters{});
using SquareDecimeters  = decltype(Decimeters{} * Decimeters{});
using SquareMeters      = decltype(Meters{} * Meters{});
using SquareKilometers  = decltype(Kilometers{} * Kilometers{});

//--------------------------------------------------------------------------
// Area per Length

using SquareCentimetersPerMeter
    = Tagged<decltype(SquareCentimeters{} / Meters{}), kinds::AreaPerLength>;

using SquareMetersPerMeter
    = Tagged<decltype(SquareMeters{}      / Meters{}), kinds::AreaPerLength>;


//------------------------------------------------------------------------------
// Volume

using CubicMillimeters  = decltype(SquareMillimeters{} * Millimeters{});
using CubicCentimeters  = decltype(SquareCentimeters{} * Centimeters{});
using CubicDecimeters   = decltype(SquareDecimeters{} * Decimeters{});
using CubicMeters       = decltype(SquareMeters{} * Meters{});

//------------------------------------------------------------------------------
// Frequency

using Hertz             = decltype(Dimensionless{} / Seconds{});
using Kilohertz         = ScaledQuantity<Conversion<Ratio<1000>>, Hertz>;
using Megahertz         = ScaledQuantity<Conversion<Ratio<1000>>, Kilohertz>;
using Gigahertz         = ScaledQuantity<Conversion<Ratio<1000>>, Megahertz>;

//------------------------------------------------------------------------------
// Velocity

using MetersPerSecond   = decltype(Meters{} / Seconds{});

using KilometersPerHour = decltype(Kilometers{} / Hours{});

//------------------------------------------------------------------------------
// Solid angle

using Steradians        = Tagged<decltype(Radians{} * Radians{}), kinds::SolidAngle>;
using SquareDegrees     = Tagged<decltype(Degrees{} * Degrees{}), kinds::SolidAngle>;

//--------------------------------------------------------------------------
// Force

using Newtons           = decltype(Kilograms{} * Meters{} / (Seconds{} * Seconds{}));
using Kilonewtons       = ScaledQuantity<Conversion<Ratio<1000>>, Newtons>;

//--------------------------------------------------------------------------
// Energy

using Joules            = decltype(Newtons{} * Meters{});
using Kilojoules        = ScaledQuantity<Conversion<Ratio<1000>>, Joules>;

//--------------------------------------------------------------------------
// Torque

using NewtonMeters      = decltype(Newtons{} * Meters{} / Radians{});

//--------------------------------------------------------------------------
// Power

using Watts             = decltype(Joules{} / Seconds{});
using Kilowatts         = ScaledQuantity<Conversion<Ratio<1000>>, Watts>;

//using Vars              = Tagged<Watts, kinds::ElectricPower>;
//using Kilovars          = ScaledQuantity<Conversion<Ratio<1000>>, Vars>;

//--------------------------------------------------------------------------
// Data

using Bits              = Quantity<units::Bit>;
using Nibbles           = ScaledQuantity<Conversion<Ratio<4>>, Bits>;
using Bytes             = ScaledQuantity<Conversion<Ratio<8>>, Bits>;
using Kilobytes         = ScaledQuantity<Conversion<Ratio<1000>>, Bytes>;
using Megabytes         = ScaledQuantity<Conversion<Ratio<1000>>, Kilobytes>;
using Gigabytes         = ScaledQuantity<Conversion<Ratio<1000>>, Megabytes>;

//==================================================================================================
//
//==================================================================================================

#if 0
[[nodiscard]] constexpr Dimensionless operator+(Dimensionless lhs, double rhs) noexcept
{
    return Dimensionless(lhs.count_internal() + rhs);
}

[[nodiscard]] constexpr Dimensionless operator-(Dimensionless lhs, double rhs) noexcept
{
    return Dimensionless(lhs.count_internal() - rhs);
}

[[nodiscard]] constexpr Dimensionless operator+(double lhs, Dimensionless rhs) noexcept
{
    return Dimensionless(lhs + rhs.count_internal());
}

[[nodiscard]] constexpr Dimensionless operator-(double lhs, Dimensionless rhs) noexcept
{
    return Dimensionless(lhs - rhs.count_internal());
}
#endif

//==================================================================================================
// QuantityPoint
//==================================================================================================

template <typename DifferenceType>
class QuantityPoint final
{
public:
    using difference_type = DifferenceType;
    using scalar_type     = typename DifferenceType::scalar_type;
    using unit            = typename DifferenceType::unit;
    using conversion      = typename DifferenceType::conversion;
    using kind            = typename DifferenceType::kind;
    using dimension       = typename DifferenceType::dimension;

private:
    difference_type _value;

public:
    constexpr QuantityPoint() noexcept = default;

    constexpr explicit QuantityPoint(difference_type value) noexcept
        : _value(value)
    {
    }

    constexpr explicit QuantityPoint(scalar_type value) noexcept
        : _value(value)
    {
    }

    [[nodiscard]] constexpr difference_type value() const noexcept
    {
        return _value;
    }

    [[nodiscard]] constexpr double count_internal() const noexcept
    {
        return _value.count_internal();
    }

    template <typename Q>
    [[nodiscard]] constexpr auto count() const noexcept -> decltype(_value.template count<Q>())
    {
        return _value.template count<Q>();
    }

    [[nodiscard]] constexpr friend QuantityPoint operator+(QuantityPoint lhs, difference_type rhs) noexcept
    {
        return QuantityPoint(lhs._value + rhs);
    }

    [[nodiscard]] constexpr friend QuantityPoint operator+(difference_type lhs, QuantityPoint rhs) noexcept
    {
        return QuantityPoint(lhs + rhs._value);
    }

    [[nodiscard]] constexpr friend difference_type operator-(QuantityPoint lhs, QuantityPoint rhs) noexcept
    {
        return lhs._value - rhs._value;
    }

    [[nodiscard]] constexpr friend QuantityPoint operator-(QuantityPoint lhs, difference_type rhs) noexcept
    {
        return QuantityPoint(lhs._value - rhs);
    }

    friend QuantityPoint& operator+=(QuantityPoint& lhs, difference_type rhs) noexcept
    {
        lhs._value += rhs;
        return lhs;
    }

    friend QuantityPoint& operator-=(QuantityPoint& lhs, difference_type rhs) noexcept
    {
        lhs._value -= rhs;
        return lhs;
    }

    [[nodiscard]] constexpr friend bool operator==(QuantityPoint lhs, QuantityPoint rhs) noexcept
    {
        return lhs._value == rhs._value;
    }

    [[nodiscard]] constexpr friend bool operator!=(QuantityPoint lhs, QuantityPoint rhs) noexcept
    {
        return lhs._value != rhs._value;
    }

    [[nodiscard]] constexpr friend bool operator<(QuantityPoint lhs, QuantityPoint rhs) noexcept
    {
        return lhs._value < rhs._value;
    }

    [[nodiscard]] constexpr friend bool operator>(QuantityPoint lhs, QuantityPoint rhs) noexcept
    {
        return lhs._value > rhs._value;
    }

    [[nodiscard]] constexpr friend bool operator<=(QuantityPoint lhs, QuantityPoint rhs) noexcept
    {
        return lhs._value <= rhs._value;
    }

    [[nodiscard]] constexpr friend bool operator>=(QuantityPoint lhs, QuantityPoint rhs) noexcept
    {
        return lhs._value >= rhs._value;
    }
};

} // namespace uom
