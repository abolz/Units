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

namespace impl {

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
// Dimension
//==================================================================================================

template <int64_t Num = 1, int64_t Den = 1>
using Dimension = typename std::ratio<Num, Den>::type;

template <typename D1, typename D2>
using MulDimensions = typename std::ratio_multiply<D1, D2>::type;

template <typename D1, typename D2>
using DivDimensions = typename std::ratio_divide<D1, D2>::type;

#if 0

// Returns the largest E, such that BASE^E divides PRODUCT.
inline constexpr int64_t exponent_of(int64_t base, int64_t product)
{
    UNITS_ASSERT(base >= 1);
    UNITS_ASSERT(product >= 0);

    int64_t e = 0;
    for (;;) {
        if (product % base != 0)
            break;
        product /= base;
        ++e;
    }

    return e;
}

template <int64_t Num, int64_t Den>
inline constexpr int64_t exponent_of(int64_t base, Dimension<Num, Den>)
{
    const int64_t e_num = exponent_of(base, Dimension<Num, Den>::num);
    const int64_t e_den = exponent_of(base, Dimension<Num, Den>::den);
    return e_num - e_den;
}

#endif

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
            return FNV1a(__FUNCSIG__);
#else
            return FNV1a(__PRETTY_FUNCTION__);
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
                static constexpr int64_t e1 = H1::exponent;
                static constexpr int64_t e2 = MT == MergeType::div ? -H2::exponent : H2::exponent;

                static constexpr uint64_t h1 = impl::TypeId<T1>();
                static constexpr uint64_t h2 = impl::TypeId<T2>();
                static_assert(std::is_same_v<T1, T2> || h1 != h2,
                    "collision detected");

                if constexpr (h1 < h2)
                {
                    using F = Factor<T1, e1>;
                    return Complex<F>{} + impl::Merge<MT>(impl::Tail(lhs), rhs);
                }
                else if constexpr (h1 > h2)
                {
                    using F = Factor<T2, e2>;
                    return Complex<F>{} + impl::Merge<MT>(lhs, impl::Tail(rhs));
                }
                else
                {
                    //
                    // FIXME:
                    // For 'Simple' types the exponent should always be == 1
                    //

//                  static constexpr int64_t e = e1 + e2;
                    static constexpr int64_t e = e1 + e2 == 0 ? 0 : (std::is_same_v<Simple, T1> ? 1 : (e1 + e2));
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
    struct MulTags
    {
        using type = Unwrap< decltype( impl::Merge<impl::MergeType::mul>(Wrap<T1>{}, Wrap<T2>{}) ) >;
    };

    template <typename T1, typename T2>
    struct DivTags
    {
        using type = Unwrap< decltype( impl::Merge<impl::MergeType::div>(Wrap<T1>{}, Wrap<T2>{}) ) >;
    };

#if 1
    // (reduce compile-time???)
    template <>
    struct MulTags<Simple, Simple> { using type = Simple; };

    // (reduce compile-time???)
    template <>
    struct DivTags<Simple, Simple> { using type = Simple; };
#endif

} // namespace kinds

template <typename T1, typename T2>
using MulTags = typename kinds::MulTags<T1, T2>::type;

template <typename T1, typename T2>
using DivTags = typename kinds::DivTags<T1, T2>::type;

template <typename K1, typename K2>
using MulKinds = Kind< MulTags<typename K1::tag, typename K2::tag>,
                       MulDimensions<typename K1::dimension, typename K2::dimension> >;

template <typename K1, typename K2>
using DivKinds = Kind< DivTags<typename K1::tag, typename K2::tag>,
                       DivDimensions<typename K1::dimension, typename K2::dimension> >;

namespace kinds {

    //----------------------------------------------------------------------
    // Base kinds

    // Some prime numbers:
    // 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, ...
    //                             ^^  ^^  ^^

    namespace dim
    {
        using One               = Dimension< 1>; // 1
        using Length            = Dimension< 2>; // Metre m
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
    using Length            = Kind< Simple, dim::Length            >; // Metre m
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

    // m^2
    using Area
        = MulKinds<Length, Length>;

    // m^2/m
    struct AreaPerLength
        : Kind<AreaPerLength, DivKinds<Area, Length>::dimension> {};

    // m^3
    using Volume
        = MulKinds<Area, Length>;

    // Steradian sr = rad^2
#if 1
    struct SolidAngle
        : Kind<SolidAngle, MulKinds<PlaneAngle, PlaneAngle>::dimension> {};
#else
    using SolidAngle
        = MulKinds<PlaneAngle, PlaneAngle>;
#endif

    // m/s
    using Velocity
        = DivKinds<Length, Time>;

    // (m/s)/s = m/s^2
    using Acceleration
        = DivKinds<Velocity, Time>;

    // kg (m/s)
    using Momentum
        = MulKinds<Mass, Velocity>;

    // Newton N = kg m/s^2
    using Force
        = MulKinds<Mass, Acceleration>;

    // Joule J = N m = kg m^2/s^2
    using Energy
        = MulKinds<Force, Length>;

    // Torque = N m/rad = J/rad = kg m^2/(s^2 rad)
    using Torque
        = DivKinds<Energy, PlaneAngle>;

    // J s = kg m^2/s
    using Action
        = MulKinds<Energy, Time>;

    // Watt W = J/s = kg m^2/s^3
    using Power
        = DivKinds<Energy, Time>;

    // Watt W = J/s = kg m^2/s^3 = reactive power
    struct ElectricPower
        : Kind<ElectricPower, Power::dimension> {};

    // Pascal Pa = N/m^2 = kg/(m s^2)
    using Pressure
        = DivKinds<Force, Area>;

    // kg/m^2
    using MassPerArea
        = DivKinds<Mass, Area>;

    // kg/m^3
    using MassPerVolume
        = DivKinds<Mass, Volume>;

    // m^3/m^2
    struct VolumePerArea
        : Kind<VolumePerArea, DivKinds<Volume, Area>::dimension> {};

    // Hertz Hz = 1/s
    using Frequency
        = DivKinds<One, Time>;

    // rad/s
    using AngularVelocity
        = DivKinds<PlaneAngle, Time>;

    // rad/s^2
    using AngularAcceleration
        = DivKinds<AngularVelocity, Time>;

    // Coulomb C = A s
    using ElectricCharge
        = MulKinds<ElectricCurrent, Time>;

    // Volt V = W/A = kg m^2/(s^3 A)
#if 1
    using ElectricPotentialDifference
        = DivKinds<ElectricPower, ElectricCurrent>;
#else
    using ElectricPotentialDifference
        = DivKinds<Power, ElectricCurrent>;
#endif

    // Farad F = C/V = s^4 A/(kg m^2)
    using Capacitance
        = DivKinds<ElectricCharge, ElectricPotentialDifference>;

    // Ohm = V/A = kg m^2/(s^3 A^2)
    using ElectricResistance
        = DivKinds<ElectricPotentialDifference, ElectricCurrent>;

    // Siemens = A/V = s^3 A^2/(kg m^2)
    using ElectricConductance
        = DivKinds<ElectricCurrent, ElectricPotentialDifference>;

    // Dose J/kg = m^2/s^2
    using Dose
        = DivKinds<Power, Mass>;

    // Gray Gy = J/kg = m^2/s^2
    struct AbsorbedDose
        : Kind<AbsorbedDose, Dose::dimension> {};

    // Sievert Sv = J/kg = m^2/s^2
    struct DoseEquivalent
        : Kind<DoseEquivalent, Dose::dimension> {};

    // Lumen lm = cd sr
    using LuminousFlux
        = MulKinds<LuminousIntensity, SolidAngle>;

    // Talbot lm s = cd sr s
    using LuminousEnergy
        = MulKinds<LuminousFlux, Time>;

    // Nit = cd/m^2 = lm/(m^2 sr)
    using Luminance
        = DivKinds<LuminousIntensity, Area>;

    // Lux lx = lm/m^2
    using Illuminance
        = DivKinds<LuminousFlux, Area>;
}

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
    template <typename C, int64_t E>
    struct PowConversion;

    template <typename C>
    struct PowConversion<C, 0>
    {
        using type = Conversion<Ratio<1>>;
    };

    template <typename C>
    struct PowConversion<C, 1>
    {
        using type = C;
    };

    template <typename C, int64_t E>
    struct PowConversion
    {
        static_assert(E >= 2);
        using type = MulConversions<C, typename PowConversion<C, E - 1>::type>;
    };

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

template <typename C, int64_t E>
using PowConversion = typename impl::PowConversion<C, E>::type;

//==================================================================================================
// Unit
//==================================================================================================

template <typename C, typename K/*, ctstring Symbol*/>
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
    template <typename C, typename U, typename K = typename U::kind>
    using ScaledUnit = typename Unit<MulConversions<C, typename U::conversion>, K>::type;

    //--------------------------------------------------------------------------
    // One

    // 1
    using One               = Unit<Conversion<Ratio<1>>, kinds::One>;
    using Dimensionless     = Unit<Conversion<Ratio<1>>, kinds::One>;
    using Percent           = ScaledUnit<Conversion<Ratio<1, 100>>, One>;
    using Permill           = ScaledUnit<Conversion<Ratio<1, 1000>>, One>;

    //--------------------------------------------------------------------------
    // Length

    using Metre             = Unit<Conversion<Ratio<1>>, kinds::Length>;
    using Millimetre        = ScaledUnit<Conversion<Ratio<1, 1000>>, Metre>;
    using Centimetre        = ScaledUnit<Conversion<Ratio<1, 100>>, Metre>;
    using Decimetre         = ScaledUnit<Conversion<Ratio<1, 10>>, Metre>;
    using Kilometre         = ScaledUnit<Conversion<Ratio<1000>>, Metre>;

    using Inch              = ScaledUnit<Conversion<Ratio<254, 100>>, Centimetre>;    // (international)
    using Foot              = ScaledUnit<Conversion<Ratio<12>>, Inch>;                // (international)
    using Yard              = ScaledUnit<Conversion<Ratio<3>>, Foot>;                 // (international)
    using Mile              = ScaledUnit<Conversion<Ratio<1760>>, Yard>;              // (international)

    //--------------------------------------------------------------------------
    // Area

    using SquareMillimetre  = Unit<PowConversion<Millimetre::conversion, 2>, kinds::Area>;
    using SquareCentimetre  = Unit<PowConversion<Centimetre::conversion, 2>, kinds::Area>;
    using SquareMetre       = Unit<PowConversion<Metre::conversion, 2>, kinds::Area>;
    using SquareKilometre   = Unit<PowConversion<Kilometre::conversion, 2>, kinds::Area>;

    //--------------------------------------------------------------------------
    // Volume

    using CubicMillimetre   = Unit<PowConversion<Millimetre::conversion, 3>, kinds::Volume>;
    using CubicCentimetre   = Unit<PowConversion<Centimetre::conversion, 3>, kinds::Volume>;
    using CubicDecimetre    = Unit<PowConversion<Decimetre::conversion, 3>, kinds::Volume>;
    using CubicMetre        = Unit<PowConversion<Metre::conversion, 3>, kinds::Volume>;

    //--------------------------------------------------------------------------
    // Time

    using Second            = Unit<Conversion<Ratio<1>>, kinds::Time>;
    using Millisecond       = ScaledUnit<Conversion<Ratio<1, 1000>>, Second>;
    using Minute            = ScaledUnit<Conversion<Ratio<60>>, Second>;
    using Hour              = ScaledUnit<Conversion<Ratio<60>>, Minute>;
    using Day               = ScaledUnit<Conversion<Ratio<24>>, Hour>;
    using Week              = ScaledUnit<Conversion<Ratio<7>>, Day>;
    using Year              = ScaledUnit<Conversion<Ratio<146097, 400>>, Day>;
    using Month             = ScaledUnit<Conversion<Ratio<1, 12>>, Year>;

    //--------------------------------------------------------------------------
    // Frequency

    using Hertz             = Unit<DivConversions<One::conversion, Second::conversion>, kinds::Frequency>;
    using Kilohertz         = ScaledUnit<Conversion<Ratio<1000>>, Hertz>;
    using Megahertz         = ScaledUnit<Conversion<Ratio<1000>>, Kilohertz>;
    using Gigahertz         = ScaledUnit<Conversion<Ratio<1000>>, Megahertz>;

    //--------------------------------------------------------------------------
    // Mass

    using Kilogram          = Unit<Conversion<Ratio<1>>, kinds::Mass>;
    using Gram              = ScaledUnit<Conversion<Ratio<1, 1000>>, Kilogram>;
    using Milligram         = ScaledUnit<Conversion<Ratio<1, 1000>>, Gram>;
    using Tonne             = ScaledUnit<Conversion<Ratio<1000>>, Kilogram>;

    //--------------------------------------------------------------------------
    // Velocity

    using MetrePerSecond    = Unit<DivConversions<Metre::conversion, Second::conversion>, kinds::Velocity>;

    using KilometrePerHour  = Unit<DivConversions<Kilometre::conversion, Hour::conversion>, kinds::Velocity>;

    //--------------------------------------------------------------------------
    // Acceleration

    using MetrePerSecondSquared
        = Unit<DivConversions<MetrePerSecond::conversion, Second::conversion>, kinds::Acceleration>;

    //--------------------------------------------------------------------------
    // Temperature

    using Kelvin            = Unit<Conversion<Ratio<1>>, kinds::Temperature>;
    using Celsius           = ScaledUnit<Conversion<Ratio<1>>, Kelvin>;
    using Fahrenheit        = ScaledUnit<Conversion<Ratio<5, 9>>, Kelvin>;
    using Rankine           = ScaledUnit<Conversion<Ratio<5, 9>>, Kelvin>;

    //--------------------------------------------------------------------------
    // Plane angle

    using Radian            = Unit<Conversion<Ratio<1>>, kinds::PlaneAngle>;
    using Degree            = ScaledUnit<Conversion<Ratio<1, 180>, /* pi^ */ 1>, Radian>;
    using Gon               = ScaledUnit<Conversion<Ratio<1, 200>, /* pi^ */ 1>, Radian>;
    using Revolution        = ScaledUnit<Conversion<Ratio<2,   1>, /* pi^ */ 1>, Radian>;

    using ArcMinute         = ScaledUnit<Conversion<Ratio<1, 60>>, Degree>;
    using ArcSecond         = ScaledUnit<Conversion<Ratio<1, 60>>, ArcMinute>;

    using ReciprocalRadian  = DivUnits<One, Radian>;

    //--------------------------------------------------------------------------
    // Solid angle

    using Steradian         = Unit<Conversion<Ratio<1>>, kinds::SolidAngle>;
    using SquareDegree      = Unit<PowConversion<Degree::conversion, 2>, kinds::SolidAngle>; // sq.deg = deg^2

    //--------------------------------------------------------------------------
    // Force

    using Newton            = Unit<MulConversions<Kilogram::conversion, MetrePerSecondSquared::conversion>, kinds::Force>;
    using Kilonewton        = ScaledUnit<Conversion<Ratio<1000>>, Newton>;

    //--------------------------------------------------------------------------
    // Energy

    using Joule             = Unit<MulConversions<Newton::conversion, Metre::conversion>, kinds::Energy>;

    //--------------------------------------------------------------------------
    // Torque

    using NewtonMetre       = Unit<DivConversions<Joule::conversion, Radian::conversion>, kinds::Torque>;

    //--------------------------------------------------------------------------
    // Power

    using Watt              = Unit<DivConversions<Joule::conversion, Second::conversion>, kinds::Power>;
    using Kilowatt          = ScaledUnit<Conversion<Ratio<1000>>, Watt>;

    using VoltAmpere        = Unit<Watt::conversion, kinds::ElectricPower>;

    //--------------------------------------------------------------------------
    // Data

    using Bit               = Unit<Conversion<Ratio<1>>, kinds::Bit>;
    using Nibble            = ScaledUnit<Conversion<Ratio<4>>, Bit>;
    using Byte              = ScaledUnit<Conversion<Ratio<8>>, Bit>;
    using Kilobyte          = ScaledUnit<Conversion<Ratio<1000>>, Byte>;
    using Megabyte          = ScaledUnit<Conversion<Ratio<1000>>, Kilobyte>;
    using Gigabyte          = ScaledUnit<Conversion<Ratio<1000>>, Megabyte>;

} // namespace units

//==================================================================================================
// Quantity (value + compile-time unit)
//==================================================================================================

template <typename U>
class Quantity;

namespace impl
{
    template <typename T>
    inline constexpr bool IsQuantity = false;

    template <typename U>
    inline constexpr bool IsQuantity<Quantity<U>> = true;

} // namespace impl

#if 0
template <typename Q, typename Tag>
using QuantityT // a.k.a. Change-Kind
    = Quantity<Unit<typename Q::conversion, Kind<Tag, typename Q::dimension>>>;
#endif

template <typename Q, typename Tag>
using Tagged // a.k.a. Change-Kind
    = Quantity<Unit<typename Q::conversion, Kind<Tag, typename Q::dimension>>>;

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
        : _count(DivConversions<C2, conversion>{}(q.count_unsafe()))
    {
    }

    template <typename C2, typename K2, EnableExplicitConversion<kind, K2, int> = 0>
    constexpr explicit Quantity(Quantity<Unit<C2, K2>> q) noexcept
        : _count(DivConversions<C2, conversion>{}(q.count_unsafe()))
    {
    }

    [[nodiscard]] constexpr simplified_type simplify() const noexcept
    {
        return simplified_type(_count);
    }

    template <typename Q, EnableExplicitConversion<typename Q::kind, kind, int> = 0>
    [[nodiscard]] constexpr Q as() const noexcept
    {
        return Q(*this);
    }

    // count_unsafe   ??
    // count_whatever ??
    // count_any      ??
    // count_this     ??
    // count_current  ??
    [[nodiscard]] constexpr double count_unsafe() const noexcept
    {
        return _count;
    }

    template <typename Q, EnableExplicitConversion<typename Q::kind, kind, int> = 0>
    [[nodiscard]] constexpr double count() const noexcept
    {
        if constexpr (std::is_same_v<Q, Quantity>)
            return count_unsafe();
        else
            return Q(*this).count_unsafe();
    }

    //--------------------------------------------------------------------------
    // Arithmetic

    [[nodiscard]] constexpr friend auto operator+(Quantity q) noexcept
    {
        return q;
    }

    [[nodiscard]] constexpr friend auto operator-(Quantity q) noexcept
    {
        return Quantity(-q.count_unsafe());
    }

    [[nodiscard]] constexpr friend auto operator+(Quantity lhs, Quantity rhs) noexcept
    {
        return Quantity(lhs.count_unsafe() + rhs.count_unsafe());
    }

    [[nodiscard]] constexpr friend auto operator-(Quantity lhs, Quantity rhs) noexcept
    {
        return Quantity(lhs.count_unsafe() - rhs.count_unsafe());
    }

    template <typename U2>
    [[nodiscard]] constexpr friend auto operator*(Quantity lhs, Quantity<U2> rhs) noexcept
    {
        return Quantity<MulUnits<unit, U2>>(lhs.count_unsafe() * rhs.count_unsafe());
    }

    template <typename U2>
    [[nodiscard]] constexpr friend auto operator/(Quantity lhs, Quantity<U2> rhs) noexcept
    {
        return Quantity<DivUnits<unit, U2>>(lhs.count_unsafe() / rhs.count_unsafe());
    }

    [[nodiscard]] constexpr friend auto operator*(Quantity lhs, scalar_type rhs) noexcept
    {
        return Quantity(lhs.count_unsafe() * rhs);
    }

    [[nodiscard]] constexpr friend auto operator/(Quantity lhs, scalar_type rhs) noexcept
    {
        return Quantity(lhs.count_unsafe() / rhs);
    }

    [[nodiscard]] constexpr friend auto operator*(scalar_type lhs, Quantity rhs) noexcept
    {
        return Quantity(lhs * rhs.count_unsafe());
    }

    [[nodiscard]] constexpr friend auto operator/(scalar_type lhs, Quantity rhs) noexcept
    {
        return Quantity<DivUnits<units::One, unit>>(lhs / rhs.count_unsafe());
    }

    constexpr friend Quantity& operator+=(Quantity& lhs, Quantity rhs) noexcept
    {
        lhs._count += rhs.count_unsafe();
        return lhs;
    }

    constexpr friend Quantity& operator-=(Quantity& lhs, Quantity rhs) noexcept
    {
        lhs._count -= rhs.count_unsafe();
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
        const auto x = Q(lhs).count_unsafe();
        const auto y = Q(rhs).count_unsafe();
        if (x < y)
            return -1;
        if (x > y)
            return +1;
        return 0;
    }

    template <typename C2, typename Q = CommonQuantity<C2>>
    [[nodiscard]] constexpr friend bool operator==(Quantity lhs, Quantity<Unit<C2, kind>> rhs) noexcept
    {
        return Q(lhs).count_unsafe() == Q(rhs).count_unsafe();
    }

    template <typename C2, typename Q = CommonQuantity<C2>>
    [[nodiscard]] constexpr friend bool operator!=(Quantity lhs, Quantity<Unit<C2, kind>> rhs) noexcept
    {
        return Q(lhs).count_unsafe() != Q(rhs).count_unsafe();
    }

    template <typename C2, typename Q = CommonQuantity<C2>>
    [[nodiscard]] constexpr friend bool operator<(Quantity lhs, Quantity<Unit<C2, kind>> rhs) noexcept
    {
        return Q(lhs).count_unsafe() < Q(rhs).count_unsafe();
    }

    template <typename C2, typename Q = CommonQuantity<C2>>
    [[nodiscard]] constexpr friend bool operator>(Quantity lhs, Quantity<Unit<C2, kind>> rhs) noexcept
    {
        return Q(lhs).count_unsafe() > Q(rhs).count_unsafe();
    }

    template <typename C2, typename Q = CommonQuantity<C2>>
    [[nodiscard]] constexpr friend bool operator<=(Quantity lhs, Quantity<Unit<C2, kind>> rhs) noexcept
    {
        return Q(lhs).count_unsafe() <= Q(rhs).count_unsafe();
    }

    template <typename C2, typename Q = CommonQuantity<C2>>
    [[nodiscard]] constexpr friend bool operator>=(Quantity lhs, Quantity<Unit<C2, kind>> rhs) noexcept
    {
        return Q(lhs).count_unsafe() >= Q(rhs).count_unsafe();
    }
};

#if 0
template <typename K, typename U>
    // requires K::dimension == U::dimension
constexpr QuantityT<Quantity<U>, K> kind_cast(Quantity<U> q) noexcept
{
    return QuantityT<Quantity<U>, K>(q.count_unsafe());
}
#endif

//--------------------------------------------------------------------------------------------------
// Typedefs
//--------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// One

using One               = Quantity<units::One>;
using Dimensionless     = Quantity<units::Dimensionless>;
using Percent           = Quantity<units::Percent>;
using Permill           = Quantity<units::Permill>;

//------------------------------------------------------------------------------
// Length

using Millimetres       = Quantity<units::Millimetre>;
using Centimetres       = Quantity<units::Centimetre>;
using Decimetres        = Quantity<units::Decimetre>;
using Metres            = Quantity<units::Metre>;
using Kilometres        = Quantity<units::Kilometre>;

using Inches            = Quantity<units::Inch>;
using Feet              = Quantity<units::Foot>;
using Yards             = Quantity<units::Yard>;
using Miles             = Quantity<units::Mile>;

//------------------------------------------------------------------------------
// Area

using SquareCentimetres = Quantity<units::SquareCentimetre>;
using SquareMetres      = Quantity<units::SquareMetre>;
using SquareKilometres  = Quantity<units::SquareKilometre>;

//------------------------------------------------------------------------------
// Volume

using CubicCentimetres  = Quantity<units::CubicCentimetre>;
using CubicDecimetres   = Quantity<units::CubicDecimetre>;
using CubicMetres       = Quantity<units::CubicMetre>;

//------------------------------------------------------------------------------
// Time

using Milliseconds      = Quantity<units::Millisecond>;
using Seconds           = Quantity<units::Second>;
using Minutes           = Quantity<units::Minute>;
using Hours             = Quantity<units::Hour>;
using Days              = Quantity<units::Day>;
using Weeks             = Quantity<units::Week>;
using Months            = Quantity<units::Month>;
using Years             = Quantity<units::Year>;

//------------------------------------------------------------------------------
// Frequency

using Hertz             = Quantity<units::Hertz>;
using Kilohertz         = Quantity<units::Kilohertz>;
using Megahertz         = Quantity<units::Megahertz>;
using Gigahertz         = Quantity<units::Gigahertz>;

//------------------------------------------------------------------------------
// Mass

using Grams             = Quantity<units::Gram>;
using Kilograms         = Quantity<units::Kilogram>;
using Milligrams        = Quantity<units::Milligram>;
using Tonnes            = Quantity<units::Tonne>;

//------------------------------------------------------------------------------
// Velocity

using MetresPerSecond   = Quantity<units::MetrePerSecond>; // decltype(1_m/1_s);

using KilometresPerHour = Quantity<units::KilometrePerHour>; // decltype(1_km/1_h);

//------------------------------------------------------------------------------
// Temperature

using Kelvin            = Quantity<units::Kelvin>;
using Celsius           = Quantity<units::Celsius>;
using Fahrenheit        = Quantity<units::Fahrenheit>;
using Rankine           = Quantity<units::Rankine>;

//------------------------------------------------------------------------------
// Plane angle

using Radians           = Quantity<units::Radian>;
using Degrees           = Quantity<units::Degree>;
using Gons              = Quantity<units::Gon>;
using Revolutions       = Quantity<units::Revolution>;

using ArcMinutes        = Quantity<units::ArcMinute>;
using ArcSeconds        = Quantity<units::ArcSecond>;

using ReciprocalRadians = Quantity<units::ReciprocalRadian>;

//------------------------------------------------------------------------------
// Solid angle

using Steradians        = Quantity<units::Steradian>;
using SquareDegrees     = Quantity<units::SquareDegree>;

//--------------------------------------------------------------------------
// Force

using Newtons           = Quantity<units::Newton>;
using Kilonewtons       = Quantity<units::Kilonewton>;

//--------------------------------------------------------------------------
// Energy

using Joules            = Quantity<units::Joule>;

//--------------------------------------------------------------------------
// Torque

using NewtonMetres      = Quantity<units::NewtonMetre>;

//--------------------------------------------------------------------------
// Power

using Watts             = Quantity<units::Watt>;
using Kilowatts         = Quantity<units::Kilowatt>;

using Vars              = Quantity<units::VoltAmpere>;
using Kilovars          = Quantity<units::ScaledUnit<Conversion<Ratio<1000>>, units::VoltAmpere>>;

//--------------------------------------------------------------------------
// Data

using Bits              = Quantity<units::Bit>;
using Nibbles           = Quantity<units::Nibble>;
using Bytes             = Quantity<units::Byte>;
using Kilobytes         = Quantity<units::Kilobyte>;
using Megabytes         = Quantity<units::Megabyte>;
using Gigabytes         = Quantity<units::Gigabyte>;

//--------------------------------------------------------------------------
// Area per Length

#if 1

using SquareCentimetresPerMetre
    = Tagged<decltype(SquareCentimetres{} / Metres{}), kinds::AreaPerLength>;

using SquareMetresPerMetre
    = Tagged<decltype(SquareMetres{}      / Metres{}), kinds::AreaPerLength>;

#endif

//==================================================================================================
//
//==================================================================================================

#if 0
inline constexpr Dimensionless operator+(Dimensionless lhs, double rhs) noexcept
{
    return Dimensionless(lhs.count_unsafe() + rhs);
}

inline constexpr Dimensionless operator-(Dimensionless lhs, double rhs) noexcept
{
    return Dimensionless(lhs.count_unsafe() - rhs);
}

inline constexpr Dimensionless operator+(double lhs, Dimensionless rhs) noexcept
{
    return Dimensionless(lhs + rhs.count_unsafe());
}

inline constexpr Dimensionless operator-(double lhs, Dimensionless rhs) noexcept
{
    return Dimensionless(lhs - rhs.count_unsafe());
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
