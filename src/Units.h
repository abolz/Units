// Copyright Alexander Bolz 2019
// Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at https://www.boost.org/LICENSE_1_0.txt)

#pragma once

#define UNITS_HAS_MATH() 0
#define UNITS_PRIME_DIMENSION() 0

#include "Ratio.h"
#include "Torsor.h"

#if UNITS_HAS_MATH()
#include <cmath>
#endif
#include <chrono>

#ifndef UNITS_ASSERT
#define UNITS_ASSERT(X) assert(X)
#endif

namespace sc {

//==================================================================================================
// Compile-time units
//==================================================================================================

template <typename K, typename D>
struct Kind
{
    using type = K; // **NO**: Kind
    using kind = K;
    using dimension = D;
};

template <typename Tag, typename K>
struct TaggedKind : Kind< Tag, typename K::dimension > {};

template <typename R, int64_t Exp>
struct Conversion;

template <typename C, typename K>
struct Unit;

template <typename U>
class Quantity;

namespace units
{
    struct DefaultZeroPoint
    {
        static constexpr double value = 0.0;
    };

    template <typename U>
    struct ZeroPoint : DefaultZeroPoint
    {
    };
}

template <typename U>
class QuantityPoint;

//--------------------------------------------------------------------------------------------------
// Dimension
//--------------------------------------------------------------------------------------------------

#if UNITS_PRIME_DIMENSION()
// Doesn't support rational exponents...

template <int64_t Num = 1, int64_t Den = 1>
using Dimension = Ratio<Num, Den>;

namespace dim // Base quantities
{
    // Some prime numbers:
    // 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97

    // WARNING: _1 MUST BE == 1!
    template <int _1> inline constexpr int64_t Length            = 2;
    template <int _1> inline constexpr int64_t Mass              = 3;
    template <int _1> inline constexpr int64_t Time              = 5;
    template <int _1> inline constexpr int64_t ElectricCurrent   = 7;
    template <int _1> inline constexpr int64_t Temperature       = 11;
    template <int _1> inline constexpr int64_t AmountOfSubstance = 13;
    template <int _1> inline constexpr int64_t LuminousIntensity = 17;
    template <int _1> inline constexpr int64_t Bit               = 19;
    template <int _1> inline constexpr int64_t Currency          = 23;
    template <int _1> inline constexpr int64_t Pixel             = 29;
    template <int _1> inline constexpr int64_t Dot               = 31;
}

template <typename D1, typename D2>
using MulDimensions = MulRatios<D1, D2>;

template <typename D1, typename D2>
using DivDimensions = DivRatios<D1, D2>;

#else

template <typename...>
struct Dimension
{
};

namespace dim // Base quantities
{
    // Some prime numbers:
    // 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97

    // Length^(n/d), Mass^(n/d), etc...

    template <int Num, int Den = 1> struct Length            { static constexpr int64_t id =  2; }; // Meter m
    template <int Num, int Den = 1> struct Mass              { static constexpr int64_t id =  3; }; // Kilogram kg
    template <int Num, int Den = 1> struct Time              { static constexpr int64_t id =  5; }; // Second s
    template <int Num, int Den = 1> struct ElectricCurrent   { static constexpr int64_t id =  7; }; // Ampere A
    template <int Num, int Den = 1> struct Temperature       { static constexpr int64_t id = 11; }; // Kelvin K
    template <int Num, int Den = 1> struct AmountOfSubstance { static constexpr int64_t id = 13; }; // Mole mol
    template <int Num, int Den = 1> struct LuminousIntensity { static constexpr int64_t id = 17; }; // Candela cd
    template <int Num, int Den = 1> struct Bit               { static constexpr int64_t id = 19; };
    template <int Num, int Den = 1> struct Currency          { static constexpr int64_t id = 23; }; // TODO: Euro, Dollar, etc...
    template <int Num, int Den = 1> struct Pixel             { static constexpr int64_t id = 29; };
    template <int Num, int Den = 1> struct Dot               { static constexpr int64_t id = 31; };
}

namespace impl
{
    template <typename T, typename U>
    struct IsSameBaseDimension
        : std::false_type
    {
    };

    template <template <int, int> class D, int Num1, int Den1, int Num2, int Den2>
    struct IsSameBaseDimension<D<Num1, Den1>, D<Num2, Den2>>
        : std::true_type
    {
    };

    template <template <int, int> class D, int Num, int Den>
    constexpr auto NegateDimension(D<Num, Den>) noexcept {
        return D<-Num, Den>{};
    }

    template <typename... Un>
    constexpr auto InvertDimension(Dimension<Un...>) {
        return Dimension<decltype(NegateDimension(Un{}))...>{};
    }

    template <typename... Ln, typename... Rn>
    constexpr auto Concat(Dimension<Ln...>, Dimension<Rn...>) {
        return Dimension<Ln..., Rn...>{};
    }

    constexpr auto Merge(Dimension<>, Dimension<>) {
        return Dimension<>{};
    }

    template <typename L1, typename... Ln>
    constexpr auto Merge(Dimension<L1, Ln...> lhs, Dimension<>) {
        return lhs;
    }

    template <typename R1, typename... Rn>
    constexpr auto Merge(Dimension<>, Dimension<R1, Rn...> rhs) {
        return rhs;
    }

    template <
        template <int, int> class L1, int L1Num, int L1Den, typename... Ln,
        template <int, int> class R1, int R1Num, int R1Den, typename... Rn
        >
    constexpr auto Merge(Dimension<L1<L1Num, L1Den>, Ln...> lhs, Dimension<R1<R1Num, R1Den>, Rn...> rhs)
    {
        static_assert(L1Den > 0, "invalid denominator");
        static_assert(R1Den > 0, "invalid denominator");

        constexpr int id1 = L1<L1Num, L1Den>::id;
        constexpr int id2 = R1<R1Num, R1Den>::id;
        if constexpr (id1 < id2)
        {
            return Concat(Dimension<L1<L1Num, L1Den>>{}, Merge(Dimension<Ln...>{}, rhs));
        }
        else if constexpr (id2 < id1)
        {
            return Concat(Dimension<R1<R1Num, R1Den>>{}, Merge(lhs, Dimension<Rn...>{}));
        }
        else
        {
            static_assert(IsSameBaseDimension<L1<L1Num, L1Den>, R1<R1Num, R1Den>>::value,
                "the 'id' of a base dimensions must be globally unique");

            using Sum = AddRatios< Ratio<L1Num, L1Den>, Ratio<R1Num, R1Den> >;
            if constexpr (Sum::num != 0)
                return Concat(Dimension<L1<Sum::num, Sum::den>>{}, Merge(Dimension<Ln...>{}, Dimension<Rn...>{}));
            else
                return Merge(Dimension<Ln...>{}, Dimension<Rn...>{});
        }
    }
}

//
// TODO:
// MulDimensions<D1, D2, Dn...>
//
template <typename D1, typename D2>
using MulDimensions = decltype(impl::Merge(D1{}, D2{}));

template <typename D1, typename D2>
using DivDimensions = decltype(impl::Merge(D1{}, impl::InvertDimension(D2{})));

#endif

namespace kinds
{
    template <typename K1, typename K2>
    struct Product
        : Kind< Product<K1, K2>, MulDimensions<typename K1::dimension, typename K2::dimension> >
    {
    };

    template <typename K1, typename K2>
    struct Quotient
        : Kind< Quotient<K1, K2>, DivDimensions<typename K1::dimension, typename K2::dimension> >
    {
    };

    //--------------------------------------------------------------------------
    // Base kinds:

    // 1
    struct One
        : Kind< One,
                Dimension<> > {};

    // Meter m
    struct Length
        : Kind< Length,
                Dimension<dim::Length<1>> > {};

    // Kilogram kg
    struct Mass
        : Kind< Mass,
                Dimension<dim::Mass<1>> > {};

    // Second s
    struct Time
        : Kind< Time,
                Dimension<dim::Time<1>> > {};

    // Ampere A
    struct ElectricCurrent
        : Kind< ElectricCurrent,
                Dimension<dim::ElectricCurrent<1>> > {};

    // Kelvin K
    struct Temperature
        : Kind< Temperature,
                Dimension<dim::Temperature<1>> > {};

    // Celsius °C
    struct TemperatureCelsius
        : Kind< TemperatureCelsius,
                Dimension<dim::Temperature<1>> > {};

    // Mole mol
    struct AmountOfSubstance
        : Kind< AmountOfSubstance,
                Dimension<dim::AmountOfSubstance<1>> > {};

    // Candela cd
    struct LuminousIntensity
        : Kind< LuminousIntensity,
                Dimension<dim::LuminousIntensity<1>> > {};

    // Radian rad = 1
    struct PlaneAngle
        : Kind< PlaneAngle,
                Dimension<> > {};

    // Steradian sr = 1 (= rad^2)
    struct SolidAngle
        : Kind< SolidAngle,
                Dimension<> > {};

    // Bit b
    struct Bit
        : Kind< Bit,
                Dimension<dim::Bit<1>> > {};

    // Pixel px
    struct Pixel
        : Kind< Pixel,
                Dimension<dim::Pixel<1>> > {};

    //--------------------------------------------------------------------------
    // Derived kinds:

    // m^2
    struct Area
        : Kind< Area,
                MulDimensions<Length::dimension, Length::dimension> > {};

    // m^3
    struct Volume
        : Kind< Volume,
                MulDimensions<Length::dimension, Area::dimension> > {};

    // m/s
    struct Velocity
        : Kind< Velocity,
                DivDimensions<Length::dimension, Time::dimension> > {};

    // m/s^2
    struct Acceleration
        : Kind< Acceleration,
                DivDimensions<Velocity::dimension, Time::dimension> > {};

    // Hertz Hz = 1/s
    struct Frequency
        : Kind< Frequency,
                DivDimensions<One::dimension, Time::dimension> > {};

    // kg/m^3
    struct Density
        : Kind< Density,
                DivDimensions<Mass::dimension, Volume::dimension> > {};

    // kg/m^2
    struct SurfaceDensity
        : Kind< SurfaceDensity,
                DivDimensions<Mass::dimension, Area::dimension> > {};

    // kg m/s
    struct Impulse
        : Kind< Impulse,
                MulDimensions<Mass::dimension, Velocity::dimension> > {};

    // Newton N = kg m/s^2
    struct Force
        : Kind< Force,
                MulDimensions<Mass::dimension, Acceleration::dimension> > {};

    // Joule J = N m = kg m^2/s^2
    struct Energy
        : Kind< Energy,
                MulDimensions<Force::dimension, Length::dimension> > {};

    // Torque (Moment of force) N m = kg m^2/s^2
    struct Torque
        : Kind< Torque,
                MulDimensions<Force::dimension, Length::dimension> > {};

    // Watt W = J/s = kg m^2/s^3
    struct Power
        : Kind< Power,
                DivDimensions<Energy::dimension, Time::dimension> > {};

    // Watt W = J/s = kg m^2/s^3
    struct RadiantFlux
        : Kind< RadiantFlux,
                DivDimensions<Energy::dimension, Time::dimension> > {};

    // Pascal Pa = N/m^2 = kg/(m s^2)
    struct Pressure
        : Kind< Pressure,
                DivDimensions<Force::dimension, Area::dimension> > {};

    // rad/s = 1/s
    struct AngularVelocity
        : Kind< AngularVelocity,
                DivDimensions<PlaneAngle::dimension, Time::dimension> > {};

    // rad/s^2 = 1/s^2
    struct AngularAcceleration
        : Kind< AngularAcceleration,
                DivDimensions<AngularVelocity::dimension, Time::dimension> > {};

    // J s = kg m^2 / s
    struct Action
        : Kind< Action,
                MulDimensions<Energy::dimension, Time::dimension> > {};

    // Coulomb C = s A
    struct ElectricCharge
        : Kind< ElectricCharge,
                MulDimensions<Time::dimension, ElectricCurrent::dimension> > {};

    // Volt V = W/A = kg m^2/(s^3 A)
    struct ElectricPotentialDifference
        : Kind< ElectricPotentialDifference,
                DivDimensions<Power::dimension, ElectricCurrent::dimension> > {};

    // Farad F = C/V = s^4 A/(kg m^2)
    struct Capacitance
        : Kind< Capacitance,
                DivDimensions<ElectricCharge::dimension, ElectricPotentialDifference::dimension> > {};

    // Ohm = V/A = kg m^2/(s^3 A^2)
    struct ElectricResistance
        : Kind< ElectricResistance,
                DivDimensions<ElectricPotentialDifference::dimension, ElectricCurrent::dimension> > {};

    // Siemens = A/V = s^3 A^2/(kg m^2)
    struct ElectricConductance
        : Kind< ElectricConductance,
                DivDimensions<ElectricCurrent::dimension, ElectricPotentialDifference::dimension> > {};

    // Gray Gy = J/kg = m^2/s^2
    struct AbsorbedDose
        : Kind< AbsorbedDose,
                DivDimensions<Energy::dimension, Mass::dimension> > {};

    // Sievert Sv = J/kg = m^2/s^2
    struct DoseEquivalent
        : Kind< DoseEquivalent,
                DivDimensions<Energy::dimension, Mass::dimension> > {};

    // lm s = cd sr s = talbot
    // Luminous energy is the perceived energy of light.
    struct LuminousEnergy
        : Kind< LuminousEnergy,
                MulDimensions<LuminousIntensity::dimension, MulDimensions<SolidAngle::dimension, Time::dimension>> > {};

    // Lumen lm = cd sr
    // Luminous flux/power is the change of luminous energy with time.
    struct LuminousFlux
        : Kind< LuminousFlux,
                MulDimensions<LuminousIntensity::dimension, SolidAngle::dimension> > {};

    // nit = cd/m^2 = lm/(m^2 sr)
    // Luminance is the density of luminous intensity with respect to projected area in a specified
    // direction at a specified point on a real or imaginary surface.
    struct Luminance
        : Kind< Luminance,
                DivDimensions<LuminousIntensity::dimension, Area::dimension> > {};

    // Lux lx = lm/m^2
    // Illuminance is the density of incident luminous flux with respect to area at a point on a real
    // or imaginary surface.
    struct Illuminance
        : Kind< Illuminance,
                DivDimensions<LuminousFlux::dimension, Area::dimension> > {};
}

template <typename K1, typename K2>
using MulKinds = typename kinds::Product<K1, K2>::type;

template <typename K1, typename K2>
using DivKinds = typename kinds::Quotient<K1, K2>::type;

namespace kinds
{
#if 0
    // 1 * 1 = 1
    template <> struct Product<One, One> { using type = One; };
    // A * 1 = A
    template <typename A> struct Product<A, One> { using type = A; };
    // 1 * A = A
    template <typename A> struct Product<One, A> { using type = A; };
    // A / 1 = A
    template <typename A> struct Quotient<A, One> { using type = A; };
#endif

#if 0
    // A / A = 1
    // i.e. all ratio's are created equal...
    template <typename A> struct Quotient<A, A> { using type = One; };
#endif

#if 0
    // A * (B / A) = B
    template <typename A, typename B> struct Product< A, Quotient<B, A> > { using type = B; };
    // (B / A) * A = B
    template <typename A, typename B> struct Product< Quotient<B, A>, A > { using type = B; };
    // (A * B) / A = B
    template <typename A, typename B> struct Quotient< Product<A, B>, A > { using type = B; };
    // (B * A) / A = B
    template <typename A, typename B> struct Quotient< Product<B, A>, A > { using type = B; };
#endif

    //==========================================================================
    // DANGER ZONE!
    //
    // Whether or not the simplifications here are valid is (very) domain
    // specific!
    //==========================================================================

    //--------------------------------------------------------------------------
    // Simplify **unambiguous** products

    template <> struct Product  < Length,               Length                  > { using type = Area; };
    template <> struct Product  < Area,                 Length                  > { using type = Volume; };
    template <> struct Product  < Length,               Area                    > { using type = Volume; };
    template <> struct Product  < Mass,                 Velocity                > { using type = Impulse; };
    template <> struct Product  < Velocity,             Mass                    > { using type = Impulse; };
    template <> struct Product  < Mass,                 Acceleration            > { using type = Force; };
    template <> struct Product  < Acceleration,         Mass                    > { using type = Force; };

#if 0
    // sr = rad^2
    template <> struct Product<PlaneAngle, PlaneAngle> { using type = SolidAngle; };
#endif

    // lm = cd sr
    template <> struct Product< LuminousIntensity, SolidAngle > { using type = LuminousFlux; };
    // lm = sr cd
    template <> struct Product< SolidAngle, LuminousIntensity > { using type = LuminousFlux; };
    // lm s
    template <> struct Product< LuminousFlux, Time > { using type = LuminousEnergy; };
    // s lm
    template <> struct Product< Time, LuminousFlux > { using type = LuminousEnergy; };

    //--------------------------------------------------------------------------
    // Simplify **unambiguous** quotients

    template <> struct Quotient < Length,               Time                    > { using type = Velocity; };
    template <> struct Quotient < Velocity,             Time                    > { using type = Acceleration; };
    template <> struct Quotient < Mass,                 Area                    > { using type = SurfaceDensity; };
    template <> struct Quotient < Mass,                 Volume                  > { using type = Density; };
    template <> struct Quotient < Force,                Area                    > { using type = Pressure; };
    template <> struct Quotient < PlaneAngle,           Time                    > { using type = AngularVelocity; };
    template <> struct Quotient < AngularVelocity,      Time                    > { using type = AngularAcceleration; };
    template <> struct Quotient < Impulse,              Time                    > { using type = Force; };
    template <> struct Quotient < Energy,               Time                    > { using type = Power; };

    template <> struct Quotient < Length,               MulKinds<Time, Time>    > { using type = Acceleration; };
    template <> struct Quotient < PlaneAngle,           MulKinds<Time, Time>    > { using type = AngularAcceleration; };

#if 0
    // sr / rad = rad
    template <> struct Quotient<SolidAngle, PlaneAngle> { using type = PlaneAngle; };
#endif

    // lm = W / s
    template <> struct Quotient < LuminousEnergy, Time > { using type = LuminousFlux; };
    // cd = lm / sr
    template <> struct Quotient < LuminousFlux, SolidAngle > { using type = LuminousIntensity; };
    // nt = cd / m^2
    template <> struct Quotient < LuminousIntensity, Area > { using type = Luminance; };
    // lx = lm / m^2
    template <> struct Quotient < LuminousFlux, Area > { using type = Illuminance; };
    // nt = lx / sr
    template <> struct Quotient < Illuminance, SolidAngle > { using type = Luminance; };
#if 1
    // lx = lm / (m^2 sr)
    template <> struct Quotient < LuminousFlux, MulKinds<Area, SolidAngle> > { using type = Luminance; };
    // lx = lm / (sr m^2)
    template <> struct Quotient < LuminousFlux, MulKinds<SolidAngle, Area> > { using type = Luminance; };
#endif
}

//--------------------------------------------------------------------------------------------------
// Conversion
//  value = R * pi^Exp = (Num / Den) * pi^Exp
//--------------------------------------------------------------------------------------------------

template <int64_t Num, int64_t Den = 1, int64_t Exp = 0>
using Conversion_t
    = Conversion< Ratio<Num, Den>, Exp >;

template <typename R, int64_t Exp = 0>
struct Conversion
{
    static constexpr int64_t num = R::num;
    static constexpr int64_t den = R::den;
    static constexpr int64_t exp = Exp;
    using ratio = typename R::type;

    static constexpr int64_t Two53 = 9007199254740992; // == 2^53

    static_assert(num >= -Two53,
        "invalid argument");
    static_assert(num <= Two53,
        "invalid argument");
    static_assert(den > 0,
        "invalid argument");
    static_assert(den <= Two53,
        "invalid argument");
    static_assert(impl::Abs(Exp) <= 4,
        "argument out of range (sorry, not implemented...)");

    // Returns: (x * num / den) * pi^exp
    [[nodiscard]] constexpr double operator()(double x) const noexcept {
        constexpr double Powers[] = {
            1,                 // pi^0
            3.141592653589793, // pi^1
            9.869604401089358, // pi^2
            31.00627668029982, // pi^3
            97.40909103400244, // pi^4
        };

        const double scaled = ratio{}(x);
        if constexpr (exp == 0)
            return scaled;
        else if constexpr (exp > 0)
            return scaled * Powers[exp];
        else
            return scaled / Powers[-exp];
    }
};

template <typename C1, typename C2 /* = C1 */>
using MulConversions = Conversion< MulRatios<typename C1::ratio, typename C2::ratio>, C1::exp + C2::exp >;

template <typename C1, typename C2>
using DivConversions = Conversion< DivRatios<typename C1::ratio, typename C2::ratio>, C1::exp - C2::exp >;

template <typename C>
using SquareConversion = MulConversions<C, C>;

template <typename C>
using CubicConversion = MulConversions<SquareConversion<C>, C>;

namespace impl
{
    template <typename C1, typename C2>
    struct CommonConversion
    {
        // SFINAE: missing 'type'
    };

    template <typename R1, typename R2, int64_t CommonExp>
    struct CommonConversion< Conversion<R1, CommonExp>, Conversion<R2, CommonExp> >
    {
        using type = Conversion< CommonRatio<R1, R2>, CommonExp >;
    };
}

template <typename C1, typename C2> // SFINAE
using CommonConversion = typename impl::CommonConversion<C1, C2>::type;

//--------------------------------------------------------------------------------------------------
// Unit
//--------------------------------------------------------------------------------------------------

template <typename C, typename K>
struct Unit
{
    using type = Unit;
    using conversion = C;
    using kind = K;
    using dimension = typename kind::dimension;
};

template <typename U1, typename U2>
using MulUnits = Unit<MulConversions<typename U1::conversion, typename U2::conversion>, MulKinds<typename U1::kind, typename U2::kind>>;

template <typename U1, typename U2>
using DivUnits = Unit<DivConversions<typename U1::conversion, typename U2::conversion>, DivKinds<typename U1::kind, typename U2::kind>>;

template <typename S, typename U, typename K = typename U::kind >
using ScaledUnit = Unit< MulConversions<S, typename U::conversion>, K >;

template <typename Tag, typename U>
using TaggedUnit = Unit< typename U::conversion, TaggedKind< Tag, typename U::kind > >;

//--------------------------------------------------------------------------------------------------
// Quantity (value + compile-time unit)
//--------------------------------------------------------------------------------------------------

template <typename C, typename K>
class Quantity<Unit<C, K>>
{
    //template <typename U2> friend class Quantity;

public:
    using type       = Quantity;
//  using point      = QuantityPoint<type>;
    using unit       = Unit<C, K>;
    using conversion = C;
    using kind       = K; // **NO**: typename K::kind
    using dimension  = typename kind::dimension;

private:
    double count_ = 0;

private:
    template <typename C1>
    using IsIntegral = std::bool_constant< C1::den == 1 && C1::exp == 0 >;

    template <typename C1, typename C2>
    using Divides = IsIntegral< DivConversions<C2, C1> >;

    template <typename K1, typename K2>
    using IsSameDimension = std::is_same<typename K1::dimension, typename K2::dimension>;

    template <typename K1>
    using IsDimensionless = std::is_same<typename K1::dimension, Dimension<>>;

public:
    constexpr Quantity() noexcept = default;
    constexpr Quantity(const Quantity&) noexcept = default;
    constexpr Quantity& operator=(const Quantity&) noexcept = default;

    constexpr explicit Quantity(double c) noexcept
        : count_(c)
    {
    }

    template <typename C2,
        std::enable_if_t< Divides<C, C2>::value, int > = 0>
    constexpr Quantity(Quantity<Unit<C2, K>> q) noexcept
        : count_(DivConversions<C2, C>{}(q.count()))
    {
    }

    template <typename C2, typename K2,
        std::enable_if_t< IsSameDimension<K, K2>::value, int > = 0>
    constexpr explicit Quantity(Quantity<Unit<C2, K2>> q) noexcept
        : count_(DivConversions<C2, C>{}(q.count()))
    {
    }

    [[nodiscard]] constexpr double count() const noexcept {
        return count_;
    }

    //[[nodiscard]] constexpr double value() const noexcept {
    //    return conversion{}(count_);
    //}

    template <typename C2, typename K2,
        std::enable_if_t< IsSameDimension<K, K2>::value, int > = 0>
    [[nodiscard]] constexpr auto convert_to(Unit<C2, K2>) const noexcept {
        return Quantity<Unit<C2, K>>(DivConversions<C, C2>{}(count()));
    }

    template <typename C2, typename K2,
        std::enable_if_t< IsSameDimension<K, K2>::value, int > = 0>
    [[nodiscard]] constexpr auto convert_to(Quantity<Unit<C2, K2>>) const noexcept {
        return Quantity<Unit<C2, K>>(DivConversions<C, C2>{}(count()));
    }

    //------------------------------------------------------------------------------
    // Arithmetic

    [[nodiscard]] constexpr friend auto operator+(Quantity q) noexcept {
        return q;
    }

    [[nodiscard]] constexpr friend auto operator-(Quantity q) noexcept {
        return Quantity(-q.count());
    }

    template <typename C2,
        typename Q = Quantity<Unit<CommonConversion<C, C2>, K>>>
    [[nodiscard]] constexpr friend auto operator+(Quantity lhs, Quantity<Unit<C2, K>> rhs) noexcept {
        return Q(Q(lhs).count() + Q(rhs).count());
    }

    template <typename C2,
        typename Q = Quantity<Unit<CommonConversion<C, C2>, K>>>
    [[nodiscard]] constexpr friend auto operator-(Quantity lhs, Quantity<Unit<C2, K>> rhs) noexcept {
        return Q(Q(lhs).count() - Q(rhs).count());
    }

    //
    // NB:
    // This operator is **not** commutative!
    //
    template <typename U2>
    [[nodiscard]] constexpr friend auto operator*(Quantity lhs, Quantity<U2> rhs) noexcept {
        return Quantity<MulUnits<unit, U2>>(lhs.count() * rhs.count());
    }

    template <typename U2>
    [[nodiscard]] constexpr friend auto operator/(Quantity lhs, Quantity<U2> rhs) noexcept {
        return Quantity<DivUnits<unit, U2>>(lhs.count() / rhs.count());
    }

    [[nodiscard]] constexpr friend auto operator*(Quantity lhs, double rhs) noexcept {
        return Quantity(lhs.count() * rhs);
    }

    [[nodiscard]] constexpr friend auto operator/(Quantity lhs, double rhs) noexcept {
        return Quantity(lhs.count() / rhs);
    }

    [[nodiscard]] constexpr friend auto operator*(double lhs, Quantity rhs) noexcept {
        return Quantity(lhs * rhs.count());
    }

    [[nodiscard]] constexpr friend auto operator/(double lhs, Quantity rhs) noexcept {
        return Quantity<DivUnits<Unit<Conversion_t<1>, kinds::One>, unit>>(lhs / rhs.count());
    }

    //------------------------------------------------------------------------------
    // Assignment operators

    template <typename C2,
        std::enable_if_t< Divides<C, C2>::value, int > = 0>
    constexpr friend Quantity& operator+=(Quantity& lhs, Quantity<Unit<C2, K>> rhs) noexcept {
        lhs.count_ += DivConversions<C2, C>{}(rhs.count());
        return lhs;
    }

    template <typename C2,
        std::enable_if_t< Divides<C, C2>::value, int > = 0>
    constexpr friend Quantity& operator-=(Quantity& lhs, Quantity<Unit<C2, K>> rhs) noexcept {
        lhs.count_ -= DivConversions<C2, C>{}(rhs.count());
        return lhs;
    }

    constexpr friend Quantity& operator*=(Quantity& lhs, double rhs) noexcept {
        return lhs = lhs * rhs;
    }

    constexpr friend Quantity& operator/=(Quantity& lhs, double rhs) noexcept {
        return lhs = lhs / rhs;
    }

    //------------------------------------------------------------------------------
    // Comparisons

private:
    static constexpr int CompareRep(double lhs, double rhs) noexcept {
        return (lhs < rhs) ? -1 : (rhs < lhs);
    }

public:
    template <typename C2,
        typename Q = Quantity<Unit<CommonConversion<C, C2>, K>>>
    [[nodiscard]] constexpr friend int Compare(Quantity lhs, Quantity<Unit<C2, K>> rhs) noexcept {
        return CompareRep(Q(lhs).count(), Q(rhs).count());
    }

    // lhs == rhs <==> lhs - rhs == 0
    // lhs >= rhs <==> lhs - rhs >= 0
    // etc...
    // If you can add them, you can compare them...
    //
    // But we use the correct comparison operator here, to support NaNs for floating-point representations...

    template <typename C2,
        typename Q = Quantity<Unit<CommonConversion<C, C2>, K>>>
    [[nodiscard]] constexpr friend bool operator==(Quantity lhs, Quantity<Unit<C2, K>> rhs) noexcept {
        return Q(lhs).count() == Q(rhs).count();
    }

    template <typename C2,
        typename Q = Quantity<Unit<CommonConversion<C, C2>, K>>>
    [[nodiscard]] constexpr friend bool operator!=(Quantity lhs, Quantity<Unit<C2, K>> rhs) noexcept {
        return Q(lhs).count() != Q(rhs).count();
    }

    template <typename C2,
        typename Q = Quantity<Unit<CommonConversion<C, C2>, K>>>
    [[nodiscard]] constexpr friend bool operator<(Quantity lhs, Quantity<Unit<C2, K>> rhs) noexcept {
        return Q(lhs).count() < Q(rhs).count();
    }

    template <typename C2,
        typename Q = Quantity<Unit<CommonConversion<C, C2>, K>>>
    [[nodiscard]] constexpr friend bool operator>(Quantity lhs, Quantity<Unit<C2, K>> rhs) noexcept {
        return Q(lhs).count() > Q(rhs).count();
    }

    template <typename C2,
        typename Q = Quantity<Unit<CommonConversion<C, C2>, K>>>
    [[nodiscard]] constexpr friend bool operator<=(Quantity lhs, Quantity<Unit<C2, K>> rhs) noexcept {
        return Q(lhs).count() <= Q(rhs).count();
    }

    template <typename C2,
        typename Q = Quantity<Unit<CommonConversion<C, C2>, K>>>
    [[nodiscard]] constexpr friend bool operator>=(Quantity lhs, Quantity<Unit<C2, K>> rhs) noexcept {
        return Q(lhs).count() >= Q(rhs).count();
    }

    //------------------------------------------------------------------------------
    // Math functions

#if UNITS_HAS_MATH()

    friend auto Abs(Quantity q) noexcept {
        return Quantity(std::fabs(q.count()));
    }

    friend auto Floor(Quantity q) noexcept {
        return Quantity(std::floor(q.count()));
    }

    friend auto Ceil(Quantity q) noexcept {
        return Quantity(std::ceil(q.count()));
    }

    friend auto Round(Quantity q) noexcept {
        return Quantity(std::round(q.count()));
    }

    friend auto Min(Quantity lhs, Quantity rhs) noexcept {
        return Quantity(std::fmin(lhs.count(), rhs.count()));
    }

    friend auto Max(Quantity lhs, Quantity rhs) noexcept {
        return Quantity(std::fmax(lhs.count(), rhs.count()));
    }

    template <typename U2>
    constexpr friend auto Fma(Quantity a, Quantity<U2> b, decltype(a * b) c) noexcept -> decltype(a * b)
    {
        //
        // NB:
        //
        // This function only allows arguments as the third parametre, which are implicitly
        // convertible to the exact type of the multiplication (a * b). This is more restrictive
        // than the expression (a * b + c). But otherwise, we would need to convert the
        // product (a * b) and the third argument (c) to their common type before the addition,
        // which would defeat the application of std::fma.
        //
        using Q = decltype(a * b);
#if 1 // is_constexpr_evaluated
        return Q(a.count() * b.count() + c.count());
#else
        return Q(std::fma(a.count(), b.count(), c.count()));
#endif
    }

    // exp(a b) = exp(a) exp(b)
    //  check!

    template <
        typename K2 = K, std::enable_if_t< IsDimensionless<K2>::value, int > = 0>
    friend auto Exp(Quantity q) noexcept {
        return Quantity(std::exp(q.count()));
    }

    // log_a(t) = log_b(t) / log_b(a)
    //  ... !?!?

    template <
        typename K2 = K, std::enable_if_t< IsDimensionless<K2>::value, int > = 0>
    friend auto Log(Quantity q) noexcept {
        return Quantity(std::log(q.count()));
    }

#endif
};

template <typename Tag, typename Q>
using TaggedQuantity = Quantity< TaggedUnit< Tag, typename Q::unit > >;

template <typename T, typename U> // T = Quantity or Unit
[[nodiscard]] constexpr auto quantity_cast(Quantity<U> q) noexcept -> decltype( q.convert_to(T{}) ) {
    return q.convert_to(T{});
}

template <typename NewKind, typename C, typename K>
[[nodiscard]] constexpr auto kind_cast(Quantity<Unit<C, K>> q) {
    static_assert( std::is_same<typename K::dimension, typename NewKind::dimension>::value,
        "incompatible dimensions" );
    return Quantity<Unit<C, NewKind>>(q.count());
}

template <typename C, typename K>
[[nodiscard]] constexpr auto remove_conversion(Quantity<Unit<C, K>> q) noexcept {
    return Quantity<Unit<Conversion_t<1>, K>>(C{}(q.count()));
}

//==================================================================================================
// Typedefs
//==================================================================================================

//namespace prefixes
//{
//    using Exa   = Conversion_t<1000000000000000000>;
//    using Peta  = Conversion_t<1000000000000000>;
//    using Tera  = Conversion_t<1000000000000>;
//    using Giga  = Conversion_t<1000000000>;
//    using Mega  = Conversion_t<1000000>;
//    using Kilo  = Conversion_t<1000>;
//    using Hecto = Conversion_t<100>;
//    using Deca  = Conversion_t<10>;
//    using One   = Conversion_t<1>;
//    using Deci  = Conversion_t<1, 10>;
//    using Centi = Conversion_t<1, 100>;
//    using Milli = Conversion_t<1, 1000>;
//    using Micro = Conversion_t<1, 1000000>;
//    using Nano  = Conversion_t<1, 1000000000>;
//    using Pico  = Conversion_t<1, 1000000000000>;
//    using Femto = Conversion_t<1, 1000000000000000>;
//    using Atto  = Conversion_t<1, 1000000000000000000>;
//}

//--------------------------------------------------------------------------------------------------
// Length

namespace units
{
#if 0
    struct Nanometre
        : Unit< Nanometre, Conversion_t<1, 1000000000>, kinds::Length > {};
    struct Micrometre
        : Unit< Micrometre, Conversion_t<1, 1000000>, kinds::Length > {};
    struct Millimetre
        : Unit< Millimetre, Conversion_t<1, 1000>, kinds::Length > {};
    struct Centimetre
        : Unit< Centimetre, Conversion_t<1, 100>, kinds::Length > {};
    struct Decimetre
        : Unit< Decimetre, Conversion_t<1, 10>, kinds::Length > {};
    struct Metre
        : Unit< Metre, Conversion_t<1>, kinds::Length > {};
    struct Hectometre
        : Unit< Hectometre, Conversion_t<100>, kinds::Length > {};
    struct Kilometre
        : Unit< Kilometre, Conversion_t<1000>, kinds::Length > {};
#else
    using Nanometre  = Unit< Conversion_t<1, 1000000000>, kinds::Length >;
    using Micrometre = Unit< Conversion_t<1, 1000000>, kinds::Length >;
    using Millimetre = Unit< Conversion_t<1, 1000>, kinds::Length >;
    using Centimetre = Unit< Conversion_t<1, 100>, kinds::Length >;
    using Decimetre  = Unit< Conversion_t<1, 10>, kinds::Length >;
    using Metre      = Unit< Conversion_t<1>, kinds::Length >;
    using Hectometre = Unit< Conversion_t<100>, kinds::Length >;
    using Kilometre  = Unit< Conversion_t<1000>, kinds::Length >;
#endif
}

using Nanometres  = Quantity< units::Nanometre >;
using Micrometres = Quantity< units::Micrometre >;
using Millimetres = Quantity< units::Millimetre >;
using Centimetres = Quantity< units::Centimetre >;
using Decimetres  = Quantity< units::Decimetre >;
using Metres      = Quantity< units::Metre >;
using Hectometres = Quantity< units::Hectometre >;
using Kilometres  = Quantity< units::Kilometre >;

namespace literals
{
    [[nodiscard]] constexpr auto operator""_nm(long double x) noexcept {
        return Nanometres{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_nm(unsigned long long x) noexcept {
        return Nanometres{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_mm(long double x) noexcept {
        return Millimetres{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_mm(unsigned long long x) noexcept {
        return Millimetres{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_cm(long double x) noexcept {
        return Centimetres{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_cm(unsigned long long x) noexcept {
        return Centimetres{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_dm(long double x) noexcept {
        return Decimetres{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_dm(unsigned long long x) noexcept {
        return Decimetres{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_m(long double x) noexcept {
        return Metres{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_m(unsigned long long x) noexcept {
        return Metres{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_hm(long double x) noexcept {
        return Hectometres{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_hm(unsigned long long x) noexcept {
        return Hectometres{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_km(long double x) noexcept {
        return Kilometres{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_km(unsigned long long x) noexcept {
        return Kilometres{static_cast<double>(x)};
    }
}

namespace units
{
    using Inch = ScaledUnit< Conversion_t<254, 100>, Centimetre >; // (international)
    using Foot = ScaledUnit< Conversion_t<12>, Inch >;             // (international)
    using Yard = ScaledUnit< Conversion_t<3>, Foot >;              // (international)
    using Mile = ScaledUnit< Conversion_t<1760>, Yard >;           // (international)
}

using Inches = Quantity< units::Inch >;
using Feet   = Quantity< units::Foot >;
using Yards  = Quantity< units::Yard >;
using Miles  = Quantity< units::Mile >;

namespace literals
{
    [[nodiscard]] constexpr auto operator""_in(long double x) noexcept {
        return Inches{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_in(unsigned long long x) noexcept {
        return Inches{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_ft(long double x) noexcept {
        return Feet{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_ft(unsigned long long x) noexcept {
        return Feet{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_yd(long double x) noexcept {
        return Yards{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_yd(unsigned long long x) noexcept {
        return Yards{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_mi(long double x) noexcept {
        return Miles{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_mi(unsigned long long x) noexcept {
        return Miles{static_cast<double>(x)};
    }
}

//--------------------------------------------------------------------------------------------------
// Area

namespace units
{
    using SquareCentimetre = Unit< SquareConversion<Centimetre::conversion>, kinds::Area >;
    using SquareDecimetre  = Unit< SquareConversion<Decimetre::conversion>, kinds::Area >;
    using SquareMetre      = Unit< SquareConversion<Metre::conversion>, kinds::Area >;
    using SquareKilometre  = Unit< SquareConversion<Kilometre::conversion>, kinds::Area >;

    using Hectare = Unit< SquareConversion<Hectometre::conversion>, kinds::Area >;
}

using SquareCentimetres = Quantity< units::SquareCentimetre >;
using SquareDecimetres  = Quantity< units::SquareDecimetre >;
using SquareMetres      = Quantity< units::SquareMetre >;
using SquareKilometres  = Quantity< units::SquareKilometre >;

//--------------------------------------------------------------------------------------------------
// Volume

namespace units
{
    using CubicCentimetre = Unit< CubicConversion<Centimetre::conversion>, kinds::Volume >;
    using CubicDecimetre  = Unit< CubicConversion<Decimetre::conversion>, kinds::Volume >;
    using CubicMetre      = Unit< CubicConversion<Metre::conversion>, kinds::Volume >;

    using Litre = CubicDecimetre;
}

using CubicCentimetres = Quantity< units::CubicCentimetre >;
using CubicDecimetres  = Quantity< units::CubicDecimetre >;
using CubicMetres      = Quantity< units::CubicMetre >;

//--------------------------------------------------------------------------------------------------
// Time

namespace units
{
    using Second      = Unit< Conversion_t<1>, kinds::Time >;
    using Millisecond = ScaledUnit< Conversion_t<1, 1000>, Second >;
    using Microsecond = ScaledUnit< Conversion_t<1, 1000>, Millisecond >;
    using Nanosecond  = ScaledUnit< Conversion_t<1, 1000>, Microsecond >;
    using Minute      = ScaledUnit< Conversion_t<60>, Second >;
    using Hour        = ScaledUnit< Conversion_t<60>, Minute >;
    //using Day         = ScaledUnit< Conversion_t<24>, Hour >;
    //using Week        = ScaledUnit< Conversion_t<7>, Day >;
    //using Year        = ScaledUnit< Conversion_t<146097, 400>, Day >;
    //using Month       = ScaledUnit< Conversion_t<1, 12>, Year >;
}

using Nanoseconds  = Quantity< units::Nanosecond >;
using Microseconds = Quantity< units::Microsecond >;
using Milliseconds = Quantity< units::Millisecond >;
using Seconds      = Quantity< units::Second >;
using Minutes      = Quantity< units::Minute >;
using Hours        = Quantity< units::Hour >;

namespace literals
{
    [[nodiscard]] constexpr auto operator""_ns(long double x) noexcept {
        return Nanoseconds{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_ns(unsigned long long x) noexcept {
        return Nanoseconds{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_ms(long double x) noexcept {
        return Milliseconds{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_ms(unsigned long long x) noexcept {
        return Milliseconds{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_s(long double x) noexcept {
        return Seconds{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_s(unsigned long long x) noexcept {
        return Seconds{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_min(long double x) noexcept {
        return Minutes{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_min(unsigned long long x) noexcept {
        return Minutes{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_h(long double x) noexcept {
        return Hours{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_h(unsigned long long x) noexcept {
        return Hours{static_cast<double>(x)};
    }
}

//--------------------------------------------------------------------------------------------------
// Frequency

namespace units
{
    using Hertz = Unit< DivConversions<Conversion_t<1>, Second::conversion>, kinds::Frequency >;
}

using Hertz = Quantity< units::Hertz >;

namespace literals
{
    [[nodiscard]] constexpr auto operator""_Hz(long double x) noexcept {
        return Hertz{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_Hz(unsigned long long x) noexcept {
        return Hertz{static_cast<double>(x)};
    }
}

//--------------------------------------------------------------------------------------------------
// Velocity

namespace units
{
#if 0
    using MetrePerSecond   = DivUnits<Metre, Second>;
    using KilometrePerHour = DivUnits<Kilometre, Hour>;
    using MilePerHour      = DivUnits<Mile, Hour>;
#else
    using MetrePerSecond
        = Unit< DivConversions<Metre::conversion, Second::conversion>, kinds::Velocity >;
    using KilometrePerHour
        = Unit< DivConversions<Kilometre::conversion, Hour::conversion>, kinds::Velocity >;
    using MilePerHour
        = Unit< DivConversions<Mile::conversion, Hour::conversion>, kinds::Velocity >;
#endif
}

using MetresPerSecond   = Quantity< units::MetrePerSecond >;
using KilometresPerHour = Quantity< units::KilometrePerHour >;
using MilesPerHour      = Quantity< units::MilePerHour >;

namespace literals
{
    [[nodiscard]] constexpr auto operator""_mps(long double x) noexcept {
        return MetresPerSecond{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_mps(unsigned long long x) noexcept {
        return MetresPerSecond{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_kmph(long double x) noexcept {
        return KilometresPerHour{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_kmph(unsigned long long x) noexcept {
        return KilometresPerHour{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_miph(long double x) noexcept {
        return MilesPerHour{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_miph(unsigned long long x) noexcept {
        return MilesPerHour{static_cast<double>(x)};
    }
}

//--------------------------------------------------------------------------------------------------
// Acceleration

namespace units
{
#if 0
    using MetrePerSecondSquared = DivUnits< MetrePerSecond, Second >;
#else
    using MetrePerSecondSquared = Unit< DivConversions<Metre::conversion, SquareConversion<Second::conversion>>, kinds::Acceleration >;
#endif
}

using MetresPerSecondSquared = Quantity< units::MetrePerSecondSquared >;

namespace literals
{
    [[nodiscard]] constexpr auto operator""_mpssq(long double x) noexcept {
        return MetresPerSecondSquared{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_mpssq(unsigned long long x) noexcept {
        return MetresPerSecondSquared{static_cast<double>(x)};
    }
}

//--------------------------------------------------------------------------------------------------
// Mass

namespace units
{
    using Gram     = Unit< Conversion_t<1, 1000>, kinds::Mass >;
    using Kilogram = Unit< Conversion_t<1>, kinds::Mass >;
    using Tonne    = Unit< Conversion_t<1000>, kinds::Mass >;
}

using Grams     = Quantity< units::Gram >;
using Kilograms = Quantity< units::Kilogram >;
using Tonnes    = Quantity< units::Tonne >;

namespace literals
{
    [[nodiscard]] constexpr auto operator""_g(long double x) noexcept {
        return Grams{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_g(unsigned long long x) noexcept {
        return Grams{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_kg(long double x) noexcept {
        return Kilograms{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_kg(unsigned long long x) noexcept {
        return Kilograms{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_t(long double x) noexcept {
        return Tonnes{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_t(unsigned long long x) noexcept {
        return Tonnes{static_cast<double>(x)};
    }
}

//--------------------------------------------------------------------------------------------------
// Force

namespace units
{
#if 0
    using Newton = MulUnits< Kilogram, MetrePerSecondSquared >;
#else
    using Newton = Unit< MulConversions<Kilogram::conversion, MetrePerSecondSquared::conversion>, kinds::Force >;
#endif
}

using Newton = Quantity< units::Newton >;

namespace literals
{
    [[nodiscard]] constexpr auto operator""_N(long double x) noexcept {
        return Newton{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_N(unsigned long long x) noexcept {
        return Newton{static_cast<double>(x)};
    }
}

//--------------------------------------------------------------------------------------------------
// Temperature

#if 0

namespace units
{
    struct Kelvin
        : Unit< Conversion_t<1>, kinds::Temperature > {};
    struct Celsius
        : Unit< Conversion_t<1>, kinds::Temperature > {};
    struct Fahrenheit
        : Unit< Conversion_t<5, 9>, kinds::Temperature > {};

    //using Kelvin     = Unit< Conversion_t<1>, kinds::Temperature >;
    //using Celsius    = ScaledUnit< Conversion_t<1>, Kelvin >;       // NB: implicitly convertible to Kelvin
    //using Fahrenheit = ScaledUnit< Conversion_t<5, 9>, Kelvin >;    // NB: **not** implicitly convertible to Kelvin
    //using Rankine    = ScaledUnit< Conversion_t<5, 9>, Kelvin >;    // NB: **not** implicitly convertible to Kelvin

    //template<> struct ZeroPoint<Kelvin>     { static constexpr Quantity<Kelvin> value{0}; };
    //template<> struct ZeroPoint<Celsius>    { static constexpr Quantity<Kelvin> value{273.15}; }; // 237.15 K == 0 °C
    //template<> struct ZeroPoint<Fahrenheit> { static constexpr Quantity<Kelvin> value{459.67}; }; // 459.67 K == 0 °F
}

using Kelvin     = Quantity< units::Kelvin  >;
using Celsius    = QuantityPoint< units::Celsius >;
using Fahrenheit = QuantityPoint< units::Fahrenheit >;
//using Rankine    = Quantity< units::Rankine >;

inline constexpr Celsius CelsiusZero = Celsius{-273.15};
inline constexpr Fahrenheit FahrenheitZero = Fahrenheit{-459.67};
//inline constexpr Rankine RankineZero = Rankine{0.0};

namespace literals
{
    [[nodiscard]] constexpr auto operator""_K(long double x) noexcept {
        return Kelvin{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_K(unsigned long long x) noexcept {
        return Kelvin{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_degC(long double x) noexcept {
        return Celsius{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_degC(unsigned long long x) noexcept {
        return Celsius{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_degF(long double x) noexcept {
        return Fahrenheit{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_degF(unsigned long long x) noexcept {
        return Fahrenheit{static_cast<double>(x)};
    }
    //[[nodiscard]] constexpr auto operator""_degR(long double x) noexcept {
    //    return Rankine{static_cast<double>(x)};
    //}
    //[[nodiscard]] constexpr auto operator""_degR(unsigned long long x) noexcept {
    //    return Rankine{static_cast<double>(x)};
    //}
}

#endif

//--------------------------------------------------------------------------------------------------
// Amount of substance

namespace units
{
    using Mole = Unit< Conversion_t<1>, kinds::AmountOfSubstance >;
}

using Moles = Quantity< units::Mole >;

namespace literals
{
    [[nodiscard]] constexpr auto operator""_mol(long double x) noexcept {
        return Moles{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_mol(unsigned long long x) noexcept {
        return Moles{static_cast<double>(x)};
    }
}

//--------------------------------------------------------------------------------------------------
// Plane angle

namespace units
{
    using Radian     = Unit< Conversion_t<1>, kinds::PlaneAngle >;
    using Degree     = Unit< Conversion_t<1, 180, /* pi^ */ 1>, kinds::PlaneAngle >;
    using Gon        = Unit< Conversion_t<1, 200, /* pi^ */ 1>, kinds::PlaneAngle >;
    using Revolution = Unit< Conversion_t<2,   1, /* pi^ */ 1>, kinds::PlaneAngle >;
}

using Radians     = Quantity< units::Radian >;
using Degrees     = Quantity< units::Degree >;
using Gons        = Quantity< units::Gon >;
using Revolutions = Quantity< units::Revolution >;

namespace literals
{
    [[nodiscard]] constexpr auto operator""_rad(long double x) noexcept {
        return Radians{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_rad(unsigned long long x) noexcept {
        return Radians{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_deg(long double x) noexcept {
        return Degrees{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_deg(unsigned long long x) noexcept {
        return Degrees{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_gon(long double x) noexcept {
        return Gons{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_gon(unsigned long long x) noexcept {
        return Gons{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_rev(long double x) noexcept {
        return Revolutions{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_rev(unsigned long long x) noexcept {
        return Revolutions{static_cast<double>(x)};
    }
}

//--------------------------------------------------------------------------------------------------
// Solid angle

namespace units
{
    using Steradian    = Unit< Conversion_t<1>, kinds::SolidAngle >;
    using SquareDegree = Unit< SquareConversion<Degree::conversion>, kinds::SolidAngle >; // sq.deg = deg^2
}

using Steradians = Quantity< units::Steradian >;

namespace literals
{
    [[nodiscard]] constexpr auto operator""_sr(long double x) noexcept {
        return Steradians{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_sr(unsigned long long x) noexcept {
        return Steradians{static_cast<double>(x)};
    }
}

//--------------------------------------------------------------------------------------------------
// Data

namespace units
{
    using Bit      = Unit< Conversion_t<1>, kinds::Bit >;
    using Nibble   = ScaledUnit< Conversion_t<4>, Bit >;
    using Byte     = ScaledUnit< Conversion_t<8>, Bit >;
    using Kilobyte = ScaledUnit< Conversion_t<1000>, Byte >;
    using Megabyte = ScaledUnit< Conversion_t<1000>, Kilobyte >;
    using Gigabyte = ScaledUnit< Conversion_t<1000>, Megabyte >;
}

using Bits      = Quantity< units::Bit >;
using Nibbles   = Quantity< units::Nibble >;
using Bytes     = Quantity< units::Byte >;
using Kilobytes = Quantity< units::Kilobyte >;
using Megabytes = Quantity< units::Megabyte >;
using Gigabytes = Quantity< units::Gigabyte >;

namespace literals
{
    [[nodiscard]] constexpr auto operator""_b(long double x) noexcept {
        return Bits{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_b(unsigned long long x) noexcept {
        return Bits{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_B(long double x) noexcept {
        return Bytes{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_B(unsigned long long x) noexcept {
        return Bytes{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_kB(long double x) noexcept {
        return Kilobytes{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_kB(unsigned long long x) noexcept {
        return Kilobytes{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_MB(long double x) noexcept {
        return Megabytes{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_MB(unsigned long long x) noexcept {
        return Megabytes{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_GB(long double x) noexcept {
        return Gigabytes{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_GB(unsigned long long x) noexcept {
        return Gigabytes{static_cast<double>(x)};
    }
}

//--------------------------------------------------------------------------------------------------
// Photometric

namespace units
{
    using Candela = Unit< Conversion_t<1>, kinds::LuminousIntensity >;

#if 0
    using Lumen   = MulUnits< Candela, Steradian >;
    using Talbot  = MulUnits< Lumen, Second >;
    using Nit     = DivUnits< Candela, SquareMetre >;
    using Lux     = DivUnits< Lumen, SquareMetre >;
#else
    using Lumen
        = Unit< MulConversions<Candela::conversion, Steradian::conversion>, kinds::LuminousFlux >;

    using Talbot
        = Unit< MulConversions<Lumen::conversion, Second::conversion>, kinds::LuminousEnergy >;

    using Nit
        = Unit< DivConversions<Candela::conversion, SquareMetre::conversion>, kinds::Luminance >;

    using Lux
        = Unit< DivConversions<Lumen::conversion, SquareMetre::conversion>, kinds::Illuminance >;
#endif
}

using Candelas = Quantity< units::Candela >;
using Talbots  = Quantity< units::Talbot >;
using Lumens   = Quantity< units::Lumen >;
using Nits     = Quantity< units::Nit >;
using Luxs     = Quantity< units::Lux >;

namespace literals
{
    [[nodiscard]] constexpr auto operator""_cd(long double x) noexcept {
        return Candelas{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_cd(unsigned long long x) noexcept {
        return Candelas{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_talbot(long double x) noexcept {
        return Talbots{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_talbot(unsigned long long x) noexcept {
        return Talbots{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_lm(long double x) noexcept {
        return Lumens{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_lm(unsigned long long x) noexcept {
        return Lumens{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_nt(long double x) noexcept {
        return Nits{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_nt(unsigned long long x) noexcept {
        return Nits{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_lx(long double x) noexcept {
        return Luxs{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_lx(unsigned long long x) noexcept {
        return Luxs{static_cast<double>(x)};
    }
}

} // namespace sc
