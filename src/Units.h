// Copyright Alexander Bolz 2019
// Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at https://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <type_traits>

#ifndef UNITS_ASSERT
#define UNITS_ASSERT(X) assert(X)
#endif

namespace units {

// Number type for exponents in Exponents<>
// TODO(???): rationals?
using Exponent = int;

// Integer type for (non-negative) compile-time fractions
using Natural = int64_t;

//==================================================================================================
// Compile-time units
//==================================================================================================

//namespace impl
//{
    template <Exponent Length = 0, Exponent Time = 0, Exponent Mass = 0, Exponent LuminousIntensity = 0>
    struct Exponents;

    struct DimensionTag {};

    template <typename Dim, typename Exp>
    struct Dimension : DimensionTag {
        using dimension = Dim;
        using exponents = Exp;
    };

    struct UnknownDimensionTag : DimensionTag {};

    template <typename Dim, typename Exp>
    struct UnknownDimension : UnknownDimensionTag {
        using dimension = Dim;
        using exponents = Exp;
    };

    //template <typename T>
    //struct IsExponents : std::false_type {};

    //template <Exponent L, Exponent T, Exponent M>
    //struct IsExponents<Exponents<L, T, M>>
    //    : std::true_type
    //{
    //};

    template <typename T>
    using IsDimension = std::is_base_of<DimensionTag, T>;

    template <typename T>
    using IsUnknownDimension = std::is_base_of<UnknownDimensionTag, T>;
//}

template <Natural Num, Natural Den = 1>
struct Rational;

template <typename Mul, typename Dim>
struct Unit;

template <typename U>
class Quantity;

//namespace impl
//{
    template <typename T>
    struct IsRational : std::false_type {};

    template <Natural Num, Natural Den>
    struct IsRational<Rational<Num, Den>>
        : std::true_type
    {
    };

    template <typename T>
    struct IsUnit : std::false_type {};

    template <typename Mul, typename Dim>
    struct IsUnit<Unit<Mul, Dim>>
        : std::bool_constant< IsRational<Mul>::value && IsDimension<Dim>::value >
    {
    };

    template <typename T>
    struct IsQuantity : std::false_type {};

    template <typename U>
    struct IsQuantity<Quantity<U>>
        : IsUnit<U>
    {
    };
//}

//--------------------------------------------------------------------------------------------------
// Exponents
//--------------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------
// SI base units
//--------------------------------------------------------------------------------------------------
//    Base quantity:                  Name:
//      time                            second (s)
//      length                          metre (m)
//      mass                            kilogram (kg)
//      electric current                ampere (A)
//      thermodynamic temperature       kelvin (K)
//      amount of substance             mole (mol)
//      luminous intensity              candela (cd)

//namespace impl
//{
    template <Exponent L, Exponent M, Exponent T, Exponent J>
    struct Exponents {
        static constexpr Exponent length = L;
        static constexpr Exponent mass = M;
        static constexpr Exponent time = T;
        static constexpr Exponent luminous_intensity = J;

        //--------------------------------------------------------------------------
        // Type constructors

        constexpr friend auto operator-(Exponents) noexcept {
            return Exponents<-L, -M, -T, -J>{};
        }

        template <Exponent L2, Exponent M2, Exponent T2, Exponent J2>
        constexpr friend auto operator+(Exponents, Exponents<L2, M2, T2, J2>) noexcept {
            return Exponents<L + L2, M + M2, T + T2, J + J2>{};
        }

        template <Exponent L2, Exponent M2, Exponent T2, Exponent J2>
        constexpr friend auto operator-(Exponents, Exponents<L2, M2, T2, J2>) noexcept {
            return Exponents<L - L2, M - M2, T - T2, J - J2>{};
        }
    };
//}

namespace dimensions
{
    struct Length :              Dimension< Length,              Exponents<1, 0, 0, 0>> {};
    struct Mass :                Dimension< Mass,                Exponents<0, 1, 0, 0>> {};
    struct Time :                Dimension< Time,                Exponents<0, 0, 1, 0>> {};
    // struct LuminousIntensity :   Dimension< LuminousIntensity,   Exponents<0, 0, 0, 1>> {};
    struct Area :                Dimension< Area,                decltype(Length::exponents{} + Length::exponents{})> {};
    struct Volume :              Dimension< Volume,              decltype(Length::exponents{} + Area::exponents{})> {};
    struct PlaneAngle :          Dimension< PlaneAngle,          decltype(Length::exponents{} - Length::exponents{})> {};
    struct SolidAngle :          Dimension< SolidAngle,          decltype(Area::exponents{} - Area::exponents{})> {};
    struct Velocity :            Dimension< Velocity,            decltype(Length::exponents{} - Time::exponents{})> {};
    struct Acceleration :        Dimension< Acceleration,        decltype(Velocity::exponents{} - Time::exponents{})> {};
    struct Impulse :             Dimension< Impulse,             decltype(Mass::exponents{} + Velocity::exponents{})> {}; // XXX: (Force * Time) ?
    struct Frequency :           Dimension< Frequency,           decltype(-Time::exponents{})> {};
    struct Density :             Dimension< Density,             decltype(Mass::exponents{} - Volume::exponents{})> {};
    // struct SurfaceDensity :      Dimension< SurfaceDensity,      decltype(Mass::exponents{} - Area::exponents{})> {};
    // struct MassConcentration :   Dimension< MassConcentration,   decltype(Mass::exponents{} - Volume::exponents{})> {};
    // struct SpecificVolume :      Dimension< SpecificVolume,      decltype(Volume::exponents{} - Mass::exponents{})> {};
    struct Force :               Dimension< Force,               decltype(Mass::exponents{} + Acceleration::exponents{})> {};
    struct Energy :              Dimension< Energy,              decltype(Force::exponents{} + Length::exponents{})> {};
    struct Torque :              Dimension< Torque,              decltype(Force::exponents{} + Length::exponents{})> {};
    // struct SpecificEnergy :      Dimension< SpecificEnergy,      decltype(Energy::exponents{} - Mass::exponents{})> {};
    struct Power :               Dimension< Power,               decltype(Energy::exponents{} - Time::exponents{})> {};
    struct Pressure :            Dimension< Pressure,            decltype(Force::exponents{} - Area::exponents{})> {};
    struct AngularVelocity :     Dimension< AngularVelocity,     decltype(PlaneAngle::exponents{} - Time::exponents{})> {};
    struct AngularAcceleration : Dimension< AngularAcceleration, decltype(AngularVelocity::exponents{} - Time::exponents{})> {};
    // struct Wavenumber :          Dimension< Wavenumber,          decltype(-Length::exponents{}) > {};
    // struct SurfaceTension :      Dimension< SurfaceTension,      decltype(Force::exponents{} - Length::exponents{})> {};
    // struct RadiantFlux :         Dimension< RadiantFlux,         decltype(Energy::exponents{} - Time::exponents{})> {};
    // struct LuminousFlux :        Dimension< LuminousFlux,        decltype(LuminousIntensity::exponents{} + SolidAngle::exponents{})> {};
    // struct Illuminance :         Dimension< Illuminance,         decltype(LuminousFlux::exponents{} - Area::exponents{})> {};
    // struct Irradiance :          Dimension< Irradiance,          decltype(RadiantFlux::exponents{} - Area::exponents{})> {};
    // struct RadiantEmiitance :    Dimension< RadiantEmiitance,    decltype(RadiantFlux::exponents{} - Area::exponents{})> {};
    // struct Radiance :            Dimension< Radiance,            decltype(RadiantFlux::exponents{} - SolidAngle::exponents{})> {};

    template <typename Dim>
    struct Inv
        : Dimension<Inv<Dim>, decltype(-(typename Dim::exponents{}))> {};

    template <typename Dim1, typename Dim2>
    struct Mul
        : Dimension<Mul<Dim1, Dim2>, decltype(typename Dim1::exponents{} + typename Dim2::exponents{})> {};

    template <typename Dim1, typename Dim2>
    struct Div
        : Dimension<Div<Dim1, Dim2>, decltype(typename Dim1::exponents{} - typename Dim2::exponents{})> {};

    template <typename Dim, typename Exp>
    constexpr auto Inverse(Dimension<Dim, Exp>) {
        return Inv<Dimension<Dim, Exp>>{};
    }

    // Known (and unabiguous) inverses
    constexpr auto Inverse(Time) { return Frequency{}; }

    template <typename Dim1, typename Exp1, typename Dim2, typename Exp2>
    constexpr auto operator*(Dimension<Dim1, Exp1>, Dimension<Dim2, Exp2>) {
        return Mul<Dimension<Dim1, Exp1>, Dimension<Dim2, Exp2>>{};
    }

    // Known (and unabiguous) products
    constexpr auto operator*(Length, Length) { return Area{}; }
    constexpr auto operator*(Length, Area) { return Volume{}; }
    constexpr auto operator*(Area, Length) { return Volume{}; }
    constexpr auto operator*(Mass, Velocity) { return Impulse{}; }
    constexpr auto operator*(Mass, Acceleration) { return Force{}; }
    constexpr auto operator*(Acceleration, Mass) { return Force{}; }
    constexpr auto operator*(Force, Time) { return Impulse{}; } // XXX

    template <typename Dim1, typename Exp1, typename Dim2, typename Exp2>
    constexpr auto operator/(Dimension<Dim1, Exp1>, Dimension<Dim2, Exp2>) {
        return Div<Dimension<Dim1, Exp1>, Dimension<Dim2, Exp2>>{};
    }

#if 0
    // dimension / dimension => dimension-less
    template <typename Dim, typename Exp>
    constexpr auto operator/(Dimension<Dim, Exp>, Dimension<Dim, Exp>) {
        return Dimension<Dim, Exponents<>>{};
    }
#endif

    // Known (and unabiguous) quaotients
    constexpr auto operator/(Length, Time) { return Velocity{}; }
    constexpr auto operator/(Velocity, Time) { return Acceleration{}; }
    constexpr auto operator/(Mass, Volume) { return Density{}; }
    constexpr auto operator/(Impulse, Time) { return Force{}; } // XXX
//  constexpr auto operator/(Mass, Area) { return SurfaceDensity{}; }
    constexpr auto operator/(Energy, Time) { return Power{}; }
    constexpr auto operator/(PlaneAngle, Time) { return AngularVelocity{}; }
    constexpr auto operator/(AngularVelocity, Time) { return AngularAcceleration{}; }

    // No:
//  constexpr auto operator/(Length, Length) { return PlaneAngle{}; }
//  constexpr auto operator/(Area, Area) { return SolidAngle{}; }

    template <typename Dim, typename Exp>
    constexpr auto Pow2(Dimension<Dim, Exp> d) noexcept { return d * d; }

    template <typename Dim, typename Exp>
    constexpr auto Pow3(Dimension<Dim, Exp> d) noexcept { return d * d * d; }
}

//--------------------------------------------------------------------------------------------------
// Rational
//--------------------------------------------------------------------------------------------------

//namespace impl
//{
    constexpr auto Gcd(Natural a, Natural b) noexcept {
        UNITS_ASSERT(a > 0); // static_assert
        UNITS_ASSERT(b > 0); // static_assert
        while (b > 0) {
            const auto r = a % b;
            a = b;
            b = r;
        }
        UNITS_ASSERT(a > 0); // static_assert
        return a;
    }

    constexpr auto Lcm(Natural a, Natural b) noexcept {
        const auto g = Gcd(a, b);
        UNITS_ASSERT(g > 0); // static_assert
        return (a / g) * b;
    }
//}

// template <Natural Num, Natural Den = 1, Natural Pow10 = 0, Natural PowPi = 0>
// struct Conversion;

template <Natural Num, Natural Den>
struct Rational {
    static_assert(Num > 0, "invalid argument");
    static_assert(Den > 0, "invalid argument");
    static_assert(Gcd(Num, Den) == 1, "use Ratio<> to construct (simplified) Rational's");

    static constexpr Natural num = Num;
    static constexpr Natural den = Den;
//  static constexpr double value = static_cast<double>(Num) / Den;

    constexpr double operator()(double x) const noexcept {
        return x * *this;
    }

    //--------------------------------------------------------------------------
    // Arithmetic

    constexpr friend double operator*(double lhs, Rational) noexcept {
        if constexpr (Num == 1 && Den == 1)
            return lhs;
        else if constexpr (Num == 1)
            return lhs / Den;
        else if constexpr (Den == 1)
            return lhs * Num;
        else
//          return (lhs * Num) / Den;
            return lhs * (static_cast<double>(Num) / Den);
    }

    constexpr friend double operator*(Rational, double rhs) noexcept {
        return rhs * Rational{};
    }

    constexpr friend double operator/(double lhs, Rational) noexcept {
        return lhs * Inverse(Rational{});
    }

    constexpr friend double operator/(Rational, double rhs) noexcept {
#if 1
        return (static_cast<double>(Num) / Den) / rhs;
#else
        if constexpr (Num == 1 && Den == 1)
            return 1 / rhs;
        else if constexpr (Num == 1)
            return 1 / (Den * rhs);
        else if constexpr (Den == 1)
            return Num / rhs;
        else
//          return Num / (Den * rhs);
            return (static_cast<double>(Num) / Den) / rhs;
#endif
    }

    //--------------------------------------------------------------------------
    // Type constructors

    constexpr friend auto Inverse(Rational) noexcept {
        return Rational<Den, Num>{};
    }

    template <Natural Num2, Natural Den2>
    constexpr friend auto operator*(Rational, Rational<Num2, Den2>) noexcept {
        constexpr auto S = Gcd(Num, Den2);
        constexpr auto T = Gcd(Den, Num2);
        return Rational<(Num / S) * (Num2 / T), (Den / T) * (Den2 / S)>{};
    }

    template <Natural Num2, Natural Den2>
    constexpr friend auto operator/(Rational, Rational<Num2, Den2>) noexcept {
        return Rational{} * Rational<Den2, Num2>{};
    }

    constexpr friend auto Pow2(Rational r) noexcept { return r * r; }
    constexpr friend auto Pow3(Rational r) noexcept { return r * r * r; }
};

template <Natural Num, Natural Den>
using Ratio
    = Rational<Num / Gcd(Num, Den), Den / Gcd(Num, Den)>;

template <typename Rat1, typename Rat2> // TODO: CommonConversion
using CommonRatio
    = Rational<Gcd(Rat1::num, Rat2::num), Lcm(Rat1::den, Rat2::den)>;

//--------------------------------------------------------------------------------------------------
// Unit
//--------------------------------------------------------------------------------------------------

template <typename Conv, typename Dim>
struct Unit {
    static_assert( IsRational<Conv>::value, "Conv must be a Rational type" ); // TODO... IsConversion!
    static_assert( IsDimension<Dim>::value, "Dim must be a Dimension type" );

    using conversion_type = Conv;
    using dimension_type = Dim;

    static constexpr conversion_type conversion{};
    static constexpr dimension_type dimension{};

    //--------------------------------------------------------------------------
    // Type constructors

    constexpr friend auto Inverse(Unit) noexcept {
        return Unit<decltype(Inverse(Conv{})), decltype(Inverse(Dim{}))>{};
    }

    template <typename Conv2, typename Dim2>
    constexpr friend auto operator*(Unit, Unit<Conv2, Dim2>) noexcept {
        return Unit<decltype(Conv{} * Conv2{}), decltype(Dim{} * Dim2{})>{};
    }

    template <Natural Num, Natural Den>
    constexpr friend auto operator*(Rational<Num, Den>, Unit) noexcept {
        return Unit<decltype(Rational<Num, Den>{} * Conv{}), Dim>{};
    }

    template <Natural Num, Natural Den>
    constexpr friend auto operator*(Unit, Rational<Num, Den>) noexcept {
        return Unit<decltype(Conv{} * Rational<Num, Den>{}), Dim>{};
    }

    template <typename Conv2, typename Dim2>
    constexpr friend auto operator/(Unit, Unit<Conv2, Dim2>) noexcept {
        return Unit<decltype(Conv{} / Conv2{}), decltype(Dim{} / Dim2{})>{};
    }

    template <Natural Num, Natural Den>
    constexpr friend auto operator/(Rational<Num, Den>, Unit) noexcept {
        return Unit<decltype(Rational<Num, Den>{} / Conv{}), decltype(Inverse(Dim{}))>{};
    }

    template <Natural Num, Natural Den>
    constexpr friend auto operator/(Unit, Rational<Num, Den>) noexcept {
        return Unit<decltype(Conv{} / Rational<Num, Den>{}), Dim>{};
    }

    constexpr friend auto Pow2(Unit) noexcept { return Unit<decltype(Pow2(Conv{})), decltype(Pow2(Dim{}))>{}; }
    constexpr friend auto Pow3(Unit) noexcept { return Unit<decltype(Pow3(Conv{})), decltype(Pow3(Dim{}))>{}; }
};

//--------------------------------------------------------------------------------------------------
// Quantity (value + compile-time unit)
//--------------------------------------------------------------------------------------------------

template <typename UnitType>
class Quantity {
    static_assert( IsUnit<UnitType>::value, "UnitType must be a Unit<> type");

public:
    using quantity_type = Quantity;
    using unit_type = UnitType;
    using conversion_type = typename unit_type::conversion_type;
    using dimension_type = typename unit_type::dimension_type;

private:
    double count_ = 0;

private:
    template <typename Dim>
    using IsCompatible = std::is_same<typename dimension_type::exponents, typename Dim::exponents>;

public:
    static constexpr conversion_type conversion{};
    static constexpr dimension_type dimension{};
    static constexpr unit_type unit{};

    constexpr Quantity() noexcept = default;
    constexpr Quantity(const Quantity&) noexcept = default;
    constexpr Quantity& operator=(const Quantity&) noexcept = default;

    constexpr explicit Quantity(double c) noexcept
        : count_(c)
    {
    }

    template <typename Conv>
    constexpr /*implicit*/ Quantity(Quantity<Unit<Conv, dimension_type>> q) noexcept
        : Quantity(q.convert_to(Quantity{}).count())
    {
    }

    template <typename Conv, typename Dim, std::enable_if_t< IsCompatible<Dim>::value, int > = 0>
    constexpr explicit Quantity(Quantity<Unit<Conv, Dim>> q) noexcept
        : Quantity(Quantity<Unit<Conv, dimension_type>>(q.count()))
    {
    }

    template <typename Conv, typename Dim, std::enable_if_t< !IsCompatible<Dim>::value, int > = 1>
    constexpr explicit Quantity(Quantity<Unit<Conv, Dim>> q) noexcept
        = delete;

    constexpr double count() const noexcept {
        return count_;
    }

    constexpr double value() const noexcept {
        return count_ * conversion;
    }

    template <typename Dim = dimension_type, std::enable_if_t< !IsUnknownDimension<Dim>::value, int > = 0>
    constexpr auto convert_to(unit_type) const noexcept {
        return *this;
    }

    template <typename Conv, typename Dim = dimension_type, std::enable_if_t< !IsUnknownDimension<Dim>::value, int > = 0>
    constexpr auto convert_to(Unit<Conv, dimension_type>) const noexcept {
        constexpr auto convert = conversion / Conv{};
        return Quantity<Unit<Conv, dimension_type>>(convert(count()));
    }

    template <typename Dim = dimension_type, std::enable_if_t< !IsUnknownDimension<Dim>::value, int > = 0>
    constexpr auto convert_to(quantity_type) const noexcept {
        return *this;
    }

    template <typename Conv, typename Dim = dimension_type, std::enable_if_t< !IsUnknownDimension<Dim>::value, int > = 0>
    constexpr auto convert_to(Quantity<Unit<Conv, dimension_type>>) const noexcept {
        constexpr auto convert = conversion / Conv{};
        return Quantity<Unit<Conv, dimension_type>>(convert(count()));
    }

    //------------------------------------------------------------------------------
    // Arithmetic

    constexpr friend auto operator+(Quantity q) noexcept {
        return q;
    }

    constexpr friend auto operator-(Quantity q) noexcept {
        return Quantity(-q.count());
    }

    //
    // TBD:
    //
    // Should it really be possible to add ratio's of the same kind?
    // E.g.: 1_m/1_m + 1_cm/1_cm
    // ???
    //

    template <typename Conv>
    constexpr friend auto operator+(Quantity lhs, Quantity<Unit<Conv, dimension_type>> rhs) noexcept {
        using U = Unit<CommonRatio<conversion_type, Conv>, dimension_type>;
        return Quantity<U>(lhs.convert_to(U{}).count() + rhs.convert_to(U{}).count());
    }

    template <typename Conv, typename Dim2>
    constexpr friend auto operator+(Quantity, Quantity<Unit<Conv, Dim2>>) noexcept
        = delete;

    template <typename Conv>
    constexpr friend auto operator-(Quantity lhs, Quantity<Unit<Conv, dimension_type>> rhs) noexcept {
        using U = Unit<CommonRatio<conversion_type, Conv>, dimension_type>;
        return Quantity<U>(lhs.convert_to(U{}).count() - rhs.convert_to(U{}).count());
    }

    template <typename Conv, typename Dim2>
    constexpr friend auto operator-(Quantity, Quantity<Unit<Conv, Dim2>>) noexcept
        = delete;

    template <typename Conv, typename Dim2>
    constexpr friend auto operator*(Quantity lhs, Quantity<Unit<Conv, Dim2>> rhs) noexcept {
        using U = Unit<decltype(conversion * Conv{}), decltype(dimension * Dim2{})>;
        return Quantity<U>(lhs.count() * rhs.count());
    }

    constexpr friend auto operator*(Quantity lhs, double rhs) noexcept {
        return Quantity(lhs.count() * rhs);
    }

    constexpr friend auto operator*(double lhs, Quantity rhs) noexcept {
        return Quantity(lhs * rhs.count());
    }

    template <typename Conv, typename Dim2>
    constexpr friend auto operator/(Quantity lhs, Quantity<Unit<Conv, Dim2>> rhs) noexcept {
        using U = Unit<decltype(conversion / Conv{}), decltype(dimension / Dim2{})>;
        return Quantity<U>(lhs.count() / rhs.count());
    }

    constexpr friend auto operator/(Quantity lhs, double rhs) noexcept {
        UNITS_ASSERT(rhs != 0);
        return Quantity(lhs.count() / rhs);
    }

    constexpr friend auto operator/(double lhs, Quantity rhs) noexcept {
        using U = Unit<decltype(Inverse(conversion)), decltype(Inverse(dimension))>;
        return Quantity<U>(lhs / rhs.count());
    }

    //------------------------------------------------------------------------------
    // Assignment operators

    template <typename Conv>
    constexpr friend Quantity& operator+=(Quantity& lhs, Quantity<Unit<Conv, dimension_type>> rhs) noexcept {
#if 0
        return lhs = lhs + rhs.convert_to(unit);
#else
        return lhs = lhs + rhs; // (up to two conversions... but probably the more consistent form)
#endif
    }

    template <typename Conv>
    constexpr friend Quantity& operator-=(Quantity& lhs, Quantity<Unit<Conv, dimension_type>> rhs) noexcept {
#if 0
        return lhs = lhs - rhs.convert_to(unit);
#else
        return lhs = lhs - rhs; // (up to two conversions... but probably the more consistent form)
#endif
    }

    template <typename Conv, typename Dim2>
    constexpr friend Quantity& operator*=(Quantity&, Quantity<Unit<Conv, Dim2>>) noexcept
        = delete;

    constexpr friend Quantity& operator*=(Quantity& lhs, double rhs) noexcept {
        return lhs = lhs * rhs;
    }

    template <typename Conv, typename Dim2>
    constexpr friend Quantity& operator/=(Quantity&, Quantity<Unit<Conv, Dim2>>) noexcept
        = delete;

    constexpr friend Quantity& operator/=(Quantity& lhs, double rhs) noexcept {
        return lhs = lhs / rhs;
    }

    //------------------------------------------------------------------------------
    // Comparisons

    template <typename Conv>
    constexpr friend bool operator==(Quantity lhs, Quantity<Unit<Conv, dimension_type>> rhs) noexcept {
        using U = Unit<CommonRatio<conversion_type, Conv>, dimension_type>;
        return lhs.convert_to(U{}).count() == rhs.convert_to(U{}).count();
    }

    template <typename Conv, typename Dim2>
    constexpr friend bool operator==(Quantity, Quantity<Unit<Conv, Dim2>>) noexcept
        = delete;

    template <typename Conv>
    constexpr friend bool operator!=(Quantity lhs, Quantity<Unit<Conv, dimension_type>> rhs) noexcept {
        return !(lhs == rhs);
    }

    template <typename Conv, typename Dim2>
    constexpr friend bool operator!=(Quantity, Quantity<Unit<Conv, Dim2>>) noexcept
        = delete;

    template <typename Conv>
    constexpr friend bool operator<(Quantity lhs, Quantity<Unit<Conv, dimension_type>> rhs) noexcept {
        using U = Unit<CommonRatio<conversion_type, Conv>, dimension_type>;
        return lhs.convert_to(U{}).count() < rhs.convert_to(U{}).count();
    }

    template <typename Conv, typename Dim2>
    constexpr friend bool operator<(Quantity, Quantity<Unit<Conv, Dim2>>) noexcept
        = delete;

    template <typename Conv>
    constexpr friend bool operator>(Quantity lhs, Quantity<Unit<Conv, dimension_type>> rhs) noexcept {
        return rhs < lhs;
    }

    template <typename Conv, typename Dim2>
    constexpr friend bool operator>(Quantity, Quantity<Unit<Conv, Dim2>>) noexcept
        = delete;

    template <typename Conv>
    constexpr friend bool operator<=(Quantity lhs, Quantity<Unit<Conv, dimension_type>> rhs) noexcept {
        return !(rhs < lhs);
    }

    template <typename Conv, typename Dim2>
    constexpr friend bool operator<=(Quantity, Quantity<Unit<Conv, Dim2>>) noexcept
        = delete;

    template <typename Conv>
    constexpr friend bool operator>=(Quantity lhs, Quantity<Unit<Conv, dimension_type>> rhs) noexcept {
        return !(lhs < rhs);
    }

    template <typename Conv, typename Dim2>
    constexpr friend bool operator>=(Quantity, Quantity<Unit<Conv, Dim2>>) noexcept
        = delete;
};

template <typename Q, typename Conv, typename Dim>
constexpr auto quantity_cast(Quantity<Unit<Conv, Dim>> q) noexcept -> decltype( q.convert_to(Q::unit) ) {
    return q.convert_to(Q::unit);
}

//==================================================================================================
// Typedefs
//==================================================================================================

namespace prefixes
{
//  using Exa   = Rational<1000000000000000000, 1>;
//  using Peta  = Rational<1000000000000000, 1>;
//  using Tera  = Rational<1000000000000, 1>;
    using Giga  = Rational<1000000000, 1>;
    using Mega  = Rational<1000000, 1>;
    using Kilo  = Rational<1000, 1>;
    using Hecto = Rational<100, 1>;
    using Deca  = Rational<10, 1>;
    using One   = Rational<1, 1>;
    using Deci  = Rational<1, 10>;
    using Centi = Rational<1, 100>;
    using Milli = Rational<1, 1000>;
    using Micro = Rational<1, 1000000>;
    using Nano  = Rational<1, 1000000000>;
//  using Pico  = Rational<1, 1000000000000>;
//  using Femto = Rational<1, 1000000000000000>;
//  using Atto  = Rational<1, 1000000000000000000>;
}

//namespace length
//{
    using Meter       = Unit< prefixes::One, dimensions::Length >;
    using Nanometer   = decltype(prefixes::Nano{} * Meter{});
    using Micrometer  = decltype(prefixes::Micro{} * Meter{});
    using Millimeter  = decltype(prefixes::Milli{} * Meter{});
    using Centimeter  = decltype(prefixes::Centi{} * Meter{});
    using Decimeter   = decltype(prefixes::Deci{} * Meter{});
    using Kilometer   = decltype(prefixes::Kilo{} * Meter{});

    using Nanometers  = Quantity< Nanometer >;
    using Micrometers = Quantity< Micrometer >;
    using Millimeters = Quantity< Millimeter >;
    using Centimeters = Quantity< Centimeter >;
    using Decimeters  = Quantity< Decimeter >;
    using Meters      = Quantity< Meter >;
    using Kilometers  = Quantity< Kilometer >;

    namespace literals
    {
        constexpr auto operator""_nm(long double x) noexcept {
            return Nanometers{static_cast<double>(x)};
        }
        constexpr auto operator""_nm(unsigned long long x) noexcept {
            return Nanometers{static_cast<double>(x)};
        }
        constexpr auto operator""_mm(long double x) noexcept {
            return Millimeters{static_cast<double>(x)};
        }
        constexpr auto operator""_mm(unsigned long long x) noexcept {
            return Millimeters{static_cast<double>(x)};
        }
        constexpr auto operator""_cm(long double x) noexcept {
            return Centimeters{static_cast<double>(x)};
        }
        constexpr auto operator""_cm(unsigned long long x) noexcept {
            return Centimeters{static_cast<double>(x)};
        }
        constexpr auto operator""_dm(long double x) noexcept {
            return Decimeters{static_cast<double>(x)};
        }
        constexpr auto operator""_dm(unsigned long long x) noexcept {
            return Decimeters{static_cast<double>(x)};
        }
        constexpr auto operator""_m(long double x) noexcept {
            return Meters{static_cast<double>(x)};
        }
        constexpr auto operator""_m(unsigned long long x) noexcept {
            return Meters{static_cast<double>(x)};
        }
        constexpr auto operator""_km(long double x) noexcept {
            return Kilometers{static_cast<double>(x)};
        }
        constexpr auto operator""_km(unsigned long long x) noexcept {
            return Kilometers{static_cast<double>(x)};
        }
    }

    using Inch   = decltype(Ratio<254, 100>{} * Centimeter{}); // (international)
    using Foot   = decltype(Ratio<12, 1>{} * Inch{});             // (international)
    using Yard   = decltype(Ratio<3, 1>{} * Foot{});              // (international)
    using Mile   = decltype(Ratio<1760, 1>{} * Yard{});           // (international)

    using Inches = Quantity< Inch >;
    using Feet   = Quantity< Foot >;
    using Yards  = Quantity< Yard >;
    using Miles  = Quantity< Mile >;

    namespace literals
    {
        constexpr auto operator""_in(long double x) noexcept {
            return Inches{static_cast<double>(x)};
        }
        constexpr auto operator""_in(unsigned long long x) noexcept {
            return Inches{static_cast<double>(x)};
        }
        constexpr auto operator""_ft(long double x) noexcept {
            return Feet{static_cast<double>(x)};
        }
        constexpr auto operator""_ft(unsigned long long x) noexcept {
            return Feet{static_cast<double>(x)};
        }
        constexpr auto operator""_yd(long double x) noexcept {
            return Yards{static_cast<double>(x)};
        }
        constexpr auto operator""_yd(unsigned long long x) noexcept {
            return Yards{static_cast<double>(x)};
        }
        constexpr auto operator""_mi(long double x) noexcept {
            return Miles{static_cast<double>(x)};
        }
        constexpr auto operator""_mi(unsigned long long x) noexcept {
            return Miles{static_cast<double>(x)};
        }
    }
//}

//namespace time
//{
    using Second       = Unit< prefixes::One, dimensions::Time >;
    using Nanosecond   = decltype(prefixes::Nano{} * Second{});
    using Microsecond  = decltype(prefixes::Micro{} * Second{});
    using Millisecond  = decltype(prefixes::Milli{} * Second{});
    using Minute       = decltype(Ratio<60, 1>{} * Second{});
    using Hour         = decltype(Ratio<60, 1>{} * Minute{});
    //using Day          = decltype(Ratio<24, 1>{} * Hour{});
    //using Week         = decltype(Ratio<7, 1>{} * Day{});
    //using Year         = decltype(Ratio<146097, 400>{} * Day{});
    //using Month        = decltype(Ratio<1, 12>{} * Year{});

    using Nanoseconds  = Quantity< Nanosecond >;
    using Microseconds = Quantity< Microsecond >;
    using Milliseconds = Quantity< Millisecond >;
    using Seconds      = Quantity< Second >;
    using Minutes      = Quantity< Minute >;
    using Hours        = Quantity< Hour >;
    //using Days         = Quantity< Day >;
    //using Weeks        = Quantity< Week >;
    //using Years        = Quantity< Year >;
    //using Months       = Quantity< Month >;

    namespace literals
    {
        constexpr auto operator""_ns(long double x) noexcept {
            return Nanoseconds{static_cast<double>(x)};
        }
        constexpr auto operator""_ns(unsigned long long x) noexcept {
            return Nanoseconds{static_cast<double>(x)};
        }
        constexpr auto operator""_ms(long double x) noexcept {
            return Milliseconds{static_cast<double>(x)};
        }
        constexpr auto operator""_ms(unsigned long long x) noexcept {
            return Milliseconds{static_cast<double>(x)};
        }
        constexpr auto operator""_s(long double x) noexcept {
            return Seconds{static_cast<double>(x)};
        }
        constexpr auto operator""_s(unsigned long long x) noexcept {
            return Seconds{static_cast<double>(x)};
        }
        constexpr auto operator""_min(long double x) noexcept {
            return Minutes{static_cast<double>(x)};
        }
        constexpr auto operator""_min(unsigned long long x) noexcept {
            return Minutes{static_cast<double>(x)};
        }
        constexpr auto operator""_h(long double x) noexcept {
            return Hours{static_cast<double>(x)};
        }
        constexpr auto operator""_h(unsigned long long x) noexcept {
            return Hours{static_cast<double>(x)};
        }
        //constexpr auto operator""_d(long double x) noexcept {
        //    return Days{static_cast<double>(x)};
        //}
        //constexpr auto operator""_d(unsigned long long x) noexcept {
        //    return Days{static_cast<double>(x)};
        //}
        //constexpr auto operator""_y(long double x) noexcept {
        //    return Years{static_cast<double>(x)};
        //}
        //constexpr auto operator""_y(unsigned long long x) noexcept {
        //    return Years{static_cast<double>(x)};
        //}
    }
//}

//namespace mass
//{
    //using Milligram  = Unit< prefixes::Micro, dimensions::Mass >;
    using Gram       = Unit< prefixes::Milli, dimensions::Mass >;
    using Kilogram   = Unit< prefixes::One, dimensions::Mass >;
    //using Tonne      = Unit< prefixes::Kilo, dimensions::Mass >;

    //using Milligrams = Quantity< Milligram >;
    using Grams      = Quantity< Gram >;
    using Kilograms  = Quantity< Kilogram >;
    //using Tonnes     = Quantity< Tonne >;

    namespace literals
    {
        //constexpr auto operator""_mg(long double x) noexcept {
        //    return Milligrams{static_cast<double>(x)};
        //}
        //constexpr auto operator""_mg(unsigned long long x) noexcept {
        //    return Milligrams{static_cast<double>(x)};
        //}
        constexpr auto operator""_g(long double x) noexcept {
            return Grams{static_cast<double>(x)};
        }
        constexpr auto operator""_g(unsigned long long x) noexcept {
            return Grams{static_cast<double>(x)};
        }
        constexpr auto operator""_kg(long double x) noexcept {
            return Kilograms{static_cast<double>(x)};
        }
        constexpr auto operator""_kg(unsigned long long x) noexcept {
            return Kilograms{static_cast<double>(x)};
        }
        //constexpr auto operator""_t(long double x) noexcept {
        //    return Tonnes{static_cast<double>(x)};
        //}
        //constexpr auto operator""_t(unsigned long long x) noexcept {
        //    return Tonnes{static_cast<double>(x)};
        //}
    }
//}

//namespace area
//{
    //using SquareMillimeter  = Unit< decltype(Pow2(Millimeter::conversion)), dimensions::Area >;
    using SquareCentimeter  = Unit< decltype(Pow2(Centimeter::conversion)), dimensions::Area >;
    //using SquareDecimeter   = Unit< decltype(Pow2(Decimeter::conversion)), dimensions::Area >;
    using SquareMeter       = Unit< decltype(Pow2(Meter::conversion)), dimensions::Area >;
    //using Hectare           = Unit< decltype(Ratio<10000>{} * SquareMeter::conversion), dimensions::Area >;

    //using SquareMillimeters = Quantity< SquareMillimeter >;
    using SquareCentimeters = Quantity< SquareCentimeter >;
    //using SquareDecimeters  = Quantity< SquareDecimeter >;
    using SquareMeters      = Quantity< SquareMeter >;
    //using Hectares          = Quantity< Hectare >;

    namespace literals
    {
        //constexpr auto operator""_sqmm(long double x) noexcept {
        //    return SquareMillimeters{static_cast<double>(x)};
        //}
        //constexpr auto operator""_sqmm(unsigned long long x) noexcept {
        //    return SquareMillimeters{static_cast<double>(x)};
        //}
        constexpr auto operator""_sqcm(long double x) noexcept {
            return SquareCentimeters{static_cast<double>(x)};
        }
        constexpr auto operator""_sqcm(unsigned long long x) noexcept {
            return SquareCentimeters{static_cast<double>(x)};
        }
        //constexpr auto operator""_sqdm(long double x) noexcept {
        //    return SquareDecimeters{static_cast<double>(x)};
        //}
        //constexpr auto operator""_sqdm(unsigned long long x) noexcept {
        //    return SquareDecimeters{static_cast<double>(x)};
        //}
        constexpr auto operator""_sqm(long double x) noexcept {
            return SquareMeters{static_cast<double>(x)};
        }
        constexpr auto operator""_sqm(unsigned long long x) noexcept {
            return SquareMeters{static_cast<double>(x)};
        }
        //constexpr auto operator""_ha(long double x) noexcept {
        //    return Hectares{static_cast<double>(x)};
        //}
        //constexpr auto operator""_ha(unsigned long long x) noexcept {
        //    return Hectares{static_cast<double>(x)};
        //}
    }

    //using SquareInch   = decltype(Pow2(Inch{}));                 // (international)
    //using SquareFoot   = decltype(Pow2(Foot{}));                 // (international)
    //using SquareYard   = decltype(Pow2(Yard{}));                 // (international)
    //using SquareMile   = decltype(Pow2(Mile{}));                 // (international)
    //using Acre         = decltype(Ratio<4840>{} * SquareYard{}); // (international)

    //using SquareInches = Quantity< SquareInch >;
    //using SquareFeet   = Quantity< SquareFoot >;
    //using SquareYards  = Quantity< SquareYard >;
    //using SquareMiles  = Quantity< SquareMile >;
    //using Acres        = Quantity< Acre >;

    namespace literals
    {
        //constexpr auto operator""_sqin(long double x) noexcept {
        //    return SquareInches{static_cast<double>(x)};
        //}
        //constexpr auto operator""_sqin(unsigned long long x) noexcept {
        //    return SquareInches{static_cast<double>(x)};
        //}
        //constexpr auto operator""_sqft(long double x) noexcept {
        //    return SquareFeet{static_cast<double>(x)};
        //}
        //constexpr auto operator""_sqft(unsigned long long x) noexcept {
        //    return SquareFeet{static_cast<double>(x)};
        //}
        //constexpr auto operator""_sqyd(long double x) noexcept {
        //    return SquareYards{static_cast<double>(x)};
        //}
        //constexpr auto operator""_sqyd(unsigned long long x) noexcept {
        //    return SquareYards{static_cast<double>(x)};
        //}
        //constexpr auto operator""_sqmi(long double x) noexcept {
        //    return SquareMiles{static_cast<double>(x)};
        //}
        //constexpr auto operator""_sqmi(unsigned long long x) noexcept {
        //    return SquareMiles{static_cast<double>(x)};
        //}
        //constexpr auto operator""_ac(long double x) noexcept {
        //    return Acres{static_cast<double>(x)};
        //}
        //constexpr auto operator""_ac(unsigned long long x) noexcept {
        //    return Acres{static_cast<double>(x)};
        //}
    }
//}

//namespace volume
//{
    using CubicMeter  = Unit< decltype(Pow3(Meter::conversion)), dimensions::Volume >;

    using CubicMeters = Quantity< CubicMeter >;

    namespace literals
    {
        constexpr auto operator""_cbm(long double x) noexcept {
            return CubicMeters{static_cast<double>(x)};
        }
        constexpr auto operator""_cbm(unsigned long long x) noexcept {
            return CubicMeters{static_cast<double>(x)};
        }
    }

    using Liter       = decltype(prefixes::Milli{} * CubicMeter{}); // = CubicDecimeter
    //using Milliliter  = decltype(prefixes::Milli{} * Liter{});

    using Liters      = Quantity< Liter >;
    //using Milliliters = Quantity< Milliliter >;

    namespace literals
    {
        constexpr auto operator""_l(long double x) noexcept {
            return Liters{static_cast<double>(x)};
        }
        constexpr auto operator""_l(unsigned long long x) noexcept {
            return Liters{static_cast<double>(x)};
        }
        //constexpr auto operator""_ml(long double x) noexcept {
        //    return Milliliters{static_cast<double>(x)};
        //}
        //constexpr auto operator""_ml(unsigned long long x) noexcept {
        //    return Milliliters{static_cast<double>(x)};
        //}
    }
//}

////namespace frequency
////{
//    using Hertz  = Unit< decltype(Inverse(Second::conversion)), dimensions::Frequency >;
//
//    using Hertzs = Quantity< Hertz >;
//
//    namespace literals
//    {
//        constexpr auto operator""_Hz(long double x) noexcept {
//            return Hertzs{static_cast<double>(x)};
//        }
//        constexpr auto operator""_Hz(unsigned long long x) noexcept {
//            return Hertzs{static_cast<double>(x)};
//        }
//    }
////}

//namespace velocity
//{
    using MeterPerSecond     = Unit< decltype(Meter::conversion / Second::conversion), dimensions::Velocity >;
    using KilometerPerHour   = Unit< decltype(Kilometer::conversion / Hour::conversion), dimensions::Velocity >;

    using MetersPerSeconds   = Quantity< MeterPerSecond >;
    using KilometersPerHours = Quantity< KilometerPerHour >;

    namespace literals
    {
        constexpr auto operator""_mps(long double x) noexcept {
            return MetersPerSeconds{static_cast<double>(x)};
        }
        constexpr auto operator""_mps(unsigned long long x) noexcept {
            return MetersPerSeconds{static_cast<double>(x)};
        }
        constexpr auto operator""_kmph(long double x) noexcept {
            return KilometersPerHours{static_cast<double>(x)};
        }
        constexpr auto operator""_kmph(unsigned long long x) noexcept {
            return KilometersPerHours{static_cast<double>(x)};
        }
    }

    //using MilePerSecond   = decltype(Mile{} / Second{});
    //using MilePerHour     = decltype(Mile{} / Hour{});

    //using MilesPerSeconds = Quantity< MilePerSecond >;
    //using MilesPerHours   = Quantity< MilePerHour >;

    //namespace literals
    //{
    //    constexpr auto operator""_mips(long double x) noexcept {
    //        return MilesPerSeconds{static_cast<double>(x)};
    //    }
    //    constexpr auto operator""_mips(unsigned long long x) noexcept {
    //        return MilesPerSeconds{static_cast<double>(x)};
    //    }
    //    constexpr auto operator""_miph(long double x) noexcept {
    //        return MilesPerHours{static_cast<double>(x)};
    //    }
    //    constexpr auto operator""_miph(unsigned long long x) noexcept {
    //        return MilesPerHours{static_cast<double>(x)};
    //    }
    //}
//}

////namespace acceleration
////{
//    using MeterPerSquareSecond   = Unit< decltype(MeterPerSecond::conversion / Second::conversion), dimensions::Acceleration >;
//
//    using MetersPerSquareSeconds = Quantity< MeterPerSquareSecond >;
//
//    namespace literals
//    {
//        constexpr auto operator""_mpsqs(long double x) noexcept {
//            return MetersPerSquareSeconds{static_cast<double>(x)};
//        }
//        constexpr auto operator""_mpsqs(unsigned long long x) noexcept {
//            return MetersPerSquareSeconds{static_cast<double>(x)};
//        }
//    }
////}

//namespace density
//{
//}

////namespace force
////{
//    using Newton  = Unit< decltype(Kilogram::conversion * MeterPerSquareSecond::conversion), dimensions::Force >;
//
//    using Newtons = Quantity< Newton >;
//
//    namespace literals
//    {
//        constexpr auto operator""_N(long double x) noexcept {
//            return Newtons{static_cast<double>(x)};
//        }
//        constexpr auto operator""_N(unsigned long long x) noexcept {
//            return Newtons{static_cast<double>(x)};
//        }
//    }
////}

////namespace impulse
////{
//    using NewtonSecond   = Unit< decltype(Newton::conversion / Second::conversion), dimensions::Impulse >;
//
//    using NewtonsSeconds = Quantity< NewtonSecond >;
//
//    namespace literals
//    {
//        constexpr auto operator""_Ns(long double x) noexcept {
//            return NewtonsSeconds{static_cast<double>(x)};
//        }
//        constexpr auto operator""_Ns(unsigned long long x) noexcept {
//            return NewtonsSeconds{static_cast<double>(x)};
//        }
//    }
////}

////namespace pressure
////{
//    using Pascal  = Unit< decltype(Kilogram::conversion / SquareMeter::conversion), dimensions::Pressure >;
//
//    using Pascals = Quantity< Pascal >;
//
//    namespace literals
//    {
//        constexpr auto operator""_Pa(long double x) noexcept {
//            return Pascals{static_cast<double>(x)};
//        }
//        constexpr auto operator""_Pa(unsigned long long x) noexcept {
//            return Pascals{static_cast<double>(x)};
//        }
//    }
////}

////namespace energy
////{
//    using Joule  = Unit< decltype(Newton::conversion * Meter::conversion), dimensions::Energy >;
//
//    using Joules = Quantity< Joule >;
//
//    namespace literals
//    {
//        constexpr auto operator""_J(long double x) noexcept {
//            return Joules{static_cast<double>(x)};
//        }
//        constexpr auto operator""_J(unsigned long long x) noexcept {
//            return Joules{static_cast<double>(x)};
//        }
//    }
////}

////namespace torque
////{
//    using NewtonMeter  = Unit< decltype(Newton::conversion * Meter::conversion), dimensions::Torque >;
//
//    using NewtonMeters = Quantity< NewtonMeter >;
//
//    namespace literals
//    {
//        constexpr auto operator""_Nm(long double x) noexcept {
//            return NewtonMeters{static_cast<double>(x)};
//        }
//        constexpr auto operator""_Nm(unsigned long long x) noexcept {
//            return NewtonMeters{static_cast<double>(x)};
//        }
//    }
////}

//namespace plane_angle
//{
    using Radian  = Unit< Ratio<1, 1>, dimensions::PlaneAngle >;

    using Radians = Quantity< Radian >;

    namespace literals
    {
        constexpr auto operator""_rad(long double x) noexcept {
            return Radians{static_cast<double>(x)};
        }
        constexpr auto operator""_rad(unsigned long long x) noexcept {
            return Radians{static_cast<double>(x)};
        }
    }
//}

////namespace solid_angle
////{
//    using Steradian  = Unit< Ratio<1, 1>, dimensions::SolidAngle >;
//
//    using Steradians = Quantity< Steradian >;
//
//    namespace literals
//    {
//        constexpr auto operator""_sr(long double x) noexcept {
//            return Steradians{static_cast<double>(x)};
//        }
//        constexpr auto operator""_sr(unsigned long long x) noexcept {
//            return Steradians{static_cast<double>(x)};
//        }
//    }
////}

//==================================================================================================
// Math
//==================================================================================================

// ... as friends in Quantity?!?!

//==================================================================================================
// Input/Output
//==================================================================================================

// ... as friends in Quantity/Unit?!?!

} // namespace units
