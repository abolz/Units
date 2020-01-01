// Copyright Alexander Bolz 2019
// Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at https://www.boost.org/LICENSE_1_0.txt)

#pragma once

#define UNITS_DELETE_EVERYTHING_ELSE() 0
#define UNITS_HAS_ANY() 0
#define UNITS_HAS_MATH() 0
#define UNITS_IGNORE_KIND() 0

#include <cassert>
#if UNITS_HAS_MATH()
#include <cmath>
#endif
#include <cstdint>
#include <type_traits>

#ifndef UNITS_ASSERT
#define UNITS_ASSERT(X) assert(X)
#endif

namespace sc {

//==================================================================================================
// Compile-time units
//==================================================================================================

using Exponent = int;
using Natural = int64_t;

struct Any {};

template <typename K, typename D>
struct Kind
{
    using type = Kind;
    using kind = K;
    using dimension = D;
};

template <Natural Num, Natural Den, Exponent Exp>
struct Rational;

template <typename C, typename K>
struct Unit;

template <typename U>
class Quantity;

//template <typename U>
//class QuantityPoint;

//template <typename T>
//auto treat_as_floating_point(T&&) // never implemented!
//    -> std::bool_constant<std::is_floating_point_v<std::remove_cv_t<T>>>;

//--------------------------------------------------------------------------------------------------
// Dimension
//--------------------------------------------------------------------------------------------------

template <typename... BaseDimension>
struct Dimension
{
};

namespace impl
{
    template <typename T, typename U>
    struct IsSameBaseDimension
        : std::false_type
    {
    };

    template <template <int> class D, int Exp1, int Exp2>
    struct IsSameBaseDimension<D<Exp1>, D<Exp2>>
        : std::true_type
    {
    };

    template <int X, template <int> class D, int E>
    constexpr auto MulExponent(D<E>) noexcept {
        return D<E * X>{};
    }

    template <int X, template <int> class D, int E>
    constexpr auto DivExponent(D<E>) noexcept {
        return D<E / X>{};
    }

    template <Natural S, typename... Un>
    constexpr auto ScaleDimension(Dimension<Un...>) noexcept {
        return Dimension<decltype(MulExponent<S>(Un{}))...>{};
    }

    template <typename... Un>
    constexpr auto Invert(Dimension<Un...>)
    {
        return Dimension<decltype(MulExponent<-1>(Un{}))...>{};
    }

    template <typename... Ln, typename... Rn>
    constexpr auto Concat(Dimension<Ln...>, Dimension<Rn...>)
    {
        return Dimension<Ln..., Rn...>{};
    }

    constexpr auto Merge(Dimension<>, Dimension<>)
    {
        return Dimension<>{};
    }

    template <typename L1, typename... Ln>
    constexpr auto Merge(Dimension<L1, Ln...> lhs, Dimension<>)
    {
        return lhs;
    }

    template <typename R1, typename... Rn>
    constexpr auto Merge(Dimension<>, Dimension<R1, Rn...> rhs)
    {
        return rhs;
    }

    template <
        template <int> class L1, int ExpL1, typename... Ln,
        template <int> class R1, int ExpR1, typename... Rn
        >
    constexpr auto Merge(Dimension<L1<ExpL1>, Ln...> lhs, Dimension<R1<ExpR1>, Rn...> rhs)
    {
        constexpr int id1 = L1<ExpL1>::id;
        constexpr int id2 = R1<ExpR1>::id;
        if constexpr (id1 < id2)
        {
            return Concat(Dimension<L1<ExpL1>>{}, Merge(Dimension<Ln...>{}, rhs));
        }
        else if constexpr (id2 < id1)
        {
            return Concat(Dimension<R1<ExpR1>>{}, Merge(lhs, Dimension<Rn...>{}));
        }
        else
        {
            static_assert(IsSameBaseDimension<L1<ExpL1>, R1<ExpR1>>::value,
                "the 'id' of a base dimensions must be globally unique");

            if constexpr (ExpL1 + ExpR1 != 0)
                return Concat(Dimension<L1<ExpL1 + ExpR1>>{}, Merge(Dimension<Ln...>{}, Dimension<Rn...>{}));
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
using DivDimensions = decltype(impl::Merge(D1{}, impl::Invert(D2{})));

namespace dim // Base quantities
{
    // TODO:
    // Rational exponents...

    template <int E> struct Length            { static constexpr int64_t id =  2; }; // Meter m
    template <int E> struct Mass              { static constexpr int64_t id =  3; }; // Kilogram kg
    template <int E> struct Time              { static constexpr int64_t id =  5; }; // Second s
    template <int E> struct ElectricCurrent   { static constexpr int64_t id =  7; }; // Ampere A
    template <int E> struct Temperature       { static constexpr int64_t id = 11; }; // Kelvin K
    template <int E> struct AmountOfSubstance { static constexpr int64_t id = 13; }; // Mole mol
    template <int E> struct LuminousIntensity { static constexpr int64_t id = 17; }; // Candela cd
    template <int E> struct Bit               { static constexpr int64_t id = 19; };
    template <int E> struct Currency          { static constexpr int64_t id = 23; }; // TODO: Euro, Dollar, etc...
    template <int E> struct Pixel             { static constexpr int64_t id = 29; };
    template <int E> struct Dot               { static constexpr int64_t id = 31; };
}

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

    // Typedefs:
    //  weak:     using X = Kind< Any, ... >
    //  strong:  struct X : Kind< X, ... >

    //--------------------------------------------------------------------------
    // Base kinds:

    // 1
    struct One
        : Kind< One, Dimension<> > {};

    // Meter m
    struct Length
        : Kind< Length, Dimension<dim::Length<1>> > {};

    // Kilogram kg
    struct Mass
        : Kind< Mass, Dimension<dim::Mass<1>> > {};

    // Second s
    struct Time
        : Kind< Time, Dimension<dim::Time<1>> > {};

    // Ampere A
    struct ElectricCurrent
        : Kind< ElectricCurrent, Dimension<dim::ElectricCurrent<1>> > {};

    // Kelvin K
    struct Temperature
        : Kind< Temperature, Dimension<dim::Temperature<1>> > {};

    // Mole mol
    struct AmountOfSubstance
        : Kind< AmountOfSubstance, Dimension<dim::AmountOfSubstance<1>> > {};

    // Candela cd
    struct LuminousIntensity
        : Kind< LuminousIntensity, Dimension<dim::LuminousIntensity<1>> > {};

    // Radian rad = 1
    struct PlaneAngle
        : Kind< PlaneAngle, Dimension<> > {};

    // Steradian sr = 1 (= rad^2)
    struct SolidAngle
        : Kind< SolidAngle, Dimension<> > {};

    // Bit b
    struct Bit
        : Kind< Bit, Dimension<dim::Bit<1>> > {};

    // Pixel px
    struct Pixel
        : Kind< Pixel, Dimension<dim::Pixel<1>> > {};

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
// Rational
//--------------------------------------------------------------------------------------------------

namespace impl
{
    template <typename L, typename R>
    constexpr int CompareValues(L lhs, R rhs) noexcept {
        if (lhs == rhs)
            return 0;
        return lhs < rhs ? -1 : 1;
    }

    constexpr int Sgn(Natural x) noexcept {
        return CompareValues(x, static_cast<Natural>(0));
    }

    constexpr Natural Abs(Natural x) noexcept {
        return x < 0 ? -x : x;
    }

    constexpr Natural Min(Natural x, Natural y) noexcept {
        return y < x ? y : x;
    }

    constexpr Natural Gcd(Natural a, Natural b) noexcept {
        UNITS_ASSERT(a >= 1); // static_assert
        UNITS_ASSERT(b >= 1); // static_assert
        while (b > 0) {
            const auto r = a % b;
            a = b;
            b = r;
        }
        UNITS_ASSERT(a >= 1); // static_assert
        return a;
    }

    constexpr Natural Lcm(Natural a, Natural b) noexcept {
        return (a / Gcd(a, b)) * b;
    }

#if 1
    constexpr Natural Power(Natural x, Natural n) noexcept {
        UNITS_ASSERT(x >= 1);
        UNITS_ASSERT(n >= 0);

        Natural p = 1;
        if (x > 1) {
            for ( ; n > 0; --n) {
                UNITS_ASSERT(p <= INT64_MAX / x);
                p *= x;
            }
        }

        return p;
    }

    // Returns x^n > lower (without overflow).
    constexpr bool IsPowerGreaterThan(Natural x, Natural n, Natural lower) noexcept {
        UNITS_ASSERT(x >= 1);
        UNITS_ASSERT(n >= 0);

        const auto lim = INT64_MAX / x;
        const auto max = lim < lower ? lim : lower; // = min(lim, lower)

        Natural p = 1;
        for ( ; n > 0; --n) {
            if (p > max) // p*x will overflow, or p > lower
                return true;
            p *= x;
        }

        return p > lower;
    }

    // Computes the n-th root y of x,
    // i.e. returns the largest y, such that y^n <= x
    constexpr Natural Root(Natural x, Natural n) noexcept {
        UNITS_ASSERT(x >= 1);
        UNITS_ASSERT(n >= 1);

        if (x <= 1 || n <= 1)
            return x;

        Natural lo = 1;
        Natural hi = 1 + x / n;
        // Bernoulli  ==>  x^(1/n) <= 1 + (x - 1)/n < 1 + x/n
        // Since n >= 2, hi will not overflow here.

        for (;;)
        {
            const auto y = lo + (hi - lo) / 2;
            if (y == lo)                           // hi - lo <= 1
                return y;
            else if (IsPowerGreaterThan(y, n, x))  // x < y^n
                hi = y;
            else                                   // y^n <= x
                lo = y;
        }
    }
#endif
}

template <Natural Num, Natural Den = 1, Exponent Exp = 0>
using Ratio = Rational< Num / impl::Gcd(Num, Den), Den / impl::Gcd(Num, Den), Exp >;

template <Natural Num, Natural Den = 1, Exponent Exp = 0>
struct Rational // RatPi. Delicious.
{
    // value = (Num / Den) * pi^Exp

    static constexpr Natural Two53 = 9007199254740992; // == 2^53

    static_assert(Num > 0,
        "invalid argument");
    static_assert(Num <= Two53,
        "invalid argument");
    static_assert(Den > 0,
        "invalid argument");
    static_assert(Den <= Two53,
        "invalid argument");
    static_assert(impl::Abs(Exp) <= 4,
        "argument out of range (sorry, not implemented...)");
    static_assert(impl::Gcd(Num, Den) == 1,
        "use Ratio<> to construct (reduced) Rational's");

    static constexpr Natural num = Num;
    static constexpr Natural den = Den;
    static constexpr Exponent exp = Exp;

    // Returns: x * num / den * pi^exp
    [[nodiscard]] constexpr double operator()(double x) const noexcept {
        return ApplyExp(ApplyRat(x));
    }

    template <Natural N2, Natural D2, Exponent E2>
    [[nodiscard]] constexpr friend auto operator*(Rational /*lhs*/, Rational<N2, D2, E2> /*rhs*/) noexcept {
        constexpr Natural S = impl::Gcd(num, D2);
        constexpr Natural T = impl::Gcd(den, N2);
        return Rational< (num / S) * (N2 / T), (den / T) * (D2 / S), exp + E2 >{};
    }

    template <Natural N2, Natural D2, Exponent E2>
    [[nodiscard]] constexpr friend auto operator/(Rational lhs, Rational<N2, D2, E2> /*rhs*/) noexcept {
        return lhs * Rational<D2, N2, -E2>{};
    }

    //[[nodiscard]] constexpr friend bool operator==(Rational, Rational) noexcept {
    //    return true;
    //}

    //[[nodiscard]] constexpr friend bool operator!=(Rational, Rational) noexcept {
    //    return false;
    //}

private:
    [[nodiscard]] static constexpr double ApplyRat(double x) noexcept {
        if constexpr (num == 1 && den == 1)
            return x;
        else if constexpr (num == 1)
            return x / den;
        else if constexpr (den == 1)
            return x * num;
        else
            return x * num / den;
    }

    [[nodiscard]] static constexpr double ApplyExp(double x) noexcept {
        /**static**/ constexpr double Powers[] = {
            1,                 // pi^0
            3.141592653589793, // pi^1
            9.869604401089358, // pi^2
            31.00627668029982, // pi^3
            97.40909103400244, // pi^4
        };

        if constexpr (exp == 0)
            return x;
        else if constexpr (exp > 0)
            return x * Powers[exp];
        else
            return x / Powers[-exp];
    }
};

//namespace impl
//{
//    template <typename R1, typename R2>
//    struct RationalProduct
//    {
//        static constexpr Natural Num1 = R1::num;
//        static constexpr Natural Den1 = R1::den;
//        static constexpr Natural Num2 = R2::num;
//        static constexpr Natural Den2 = R2::den;
//
//        static constexpr Natural S = Gcd(Num1, Den2);
//        static constexpr Natural T = Gcd(Den1, Num2);
//
//        using type = Rational< (Num1 / S) * (Num2 / T), (Den1 / T) * (Den2 / S), R1::exp + R2::exp >;
//    };
//}

template <typename C1, typename C2 /* = C1 */>
using MulConversions = decltype(C1{} * C2{});

template <typename C1, typename C2>
using DivConversions = decltype(C1{} / C2{});

template <typename C>
using Square = decltype(C{} * C{});

template <typename C>
using Cubic = decltype(C{} * C{} * C{});

namespace impl
{
    template <Natural N>
    using HasSquareRoot = std::bool_constant<Power(Root(N, 2), 2) == N>;

    template <Natural N, Natural D>
    using HasRationalSquareRoot = std::bool_constant<HasSquareRoot<N>::value && HasSquareRoot<D>::value>;
}

namespace impl
{
    template <typename C1, typename C2>
    struct CommonRational;

    template <Natural Num1, Natural Den1, Natural Num2, Natural Den2, Exponent CommonExp>
    struct CommonRational< Rational<Num1, Den1, CommonExp>, Rational<Num2, Den2, CommonExp> >
    {
        using type = Rational< Gcd(Num1, Num2), Lcm(Den1, Den2), CommonExp >;
    };
}

template <typename C1, typename C2>
using CommonRatio = typename impl::CommonRational<C1, C2>::type;

//--------------------------------------------------------------------------------------------------
// Unit
//--------------------------------------------------------------------------------------------------

template <typename C, typename K>
struct Unit
{
    using conversion = C;
    using kind = K;
    using dimension = typename kind::dimension;

    template <typename C2, typename K2>
    [[nodiscard]] constexpr friend auto operator*(Unit /*lhs*/, Unit<C2, K2> /*rhs*/) noexcept {
        return Unit<MulConversions<C, C2>, MulKinds<K, K2>>{};
    }

    template <typename C2, typename K2>
    [[nodiscard]] constexpr friend auto operator/(Unit /*lhs*/, Unit<C2, K2> /*rhs*/) noexcept {
        return Unit<DivConversions<C, C2>, DivKinds<K, K2>>{};
    }

    template <Natural N2, Natural D2, Natural E2>
    [[nodiscard]] constexpr friend auto operator*(Rational<N2, D2, E2> /*lhs*/, Unit /*rhs*/) noexcept {
        return Unit<MulConversions<Rational<N2, D2, E2>, C>, K>{};
    }

    template <Natural N2, Natural D2, Natural E2>
    [[nodiscard]] constexpr friend auto operator*(Unit /*lhs*/, Rational<N2, D2, E2> /*rhs*/) noexcept {
        return Unit<MulConversions<C, Rational<N2, D2, E2>>, K>{};
    }
};

//template <typename U1, typename U2>
//using MulUnits
//    = Unit< MulConversions<typename U1::conversion, typename U2::conversion>, MulKinds<typename U1::kind, typename U2::kind> >;
//
//template <typename U1, typename U2>
//using DivUnits
//    = Unit< DivConversions<typename U1::conversion, typename U2::conversion>, DivKinds<typename U1::kind, typename U2::kind> >;

template <typename U1, typename U2>
using MulUnits = decltype(U1{} * U2{});

template <typename U1, typename U2>
using DivUnits = decltype(U1{} / U2{});

template <typename S, typename U, typename K = typename U::kind >
using ScaledUnit = Unit< MulConversions<S, typename U::conversion>, K >;

//--------------------------------------------------------------------------------------------------
// Quantity (value + compile-time unit)
//--------------------------------------------------------------------------------------------------

template </*typename R, */typename C, typename K>
class Quantity</*R, */Unit<C, K>>
{
    //template <typename U2> friend class Quantity;
    //template <typename U2> friend class QuantityPoint;

public:
    using quantity = Quantity;
    using rep = double;
    using unit = Unit<C, K>;
    using conversion = C;
    using kind = K;
    using dimension = typename kind::dimension;

private:
    double count_ = 0;

private:
    template <typename R>
    using IsNatural
        = std::bool_constant< R::den == 1 && R::exp == 0 >;

    // C1 | C2
    template <typename C1, typename C2>
    using Divides
        = IsNatural< DivConversions<C2, C1> >;

    template <typename C1, typename C2, typename T = int>
    using EnableIfDivides
        = std::enable_if_t< Divides<C1, C2>::value, T >;

    template <typename K1, typename K2>
    using IsSameDimensionType
        = std::is_same< typename K1::dimension, typename K2::dimension >;

    template <typename K1, typename K2, typename T = int>
    using EnableIfCompatible
        = std::enable_if_t< IsSameDimensionType<K1, K2>::value, T >;

    // PRE: K1 == K2
    template <typename C1, typename C2, typename K2 = K>
    using CommonTypeSfinae
        = std::enable_if_t< C1::exp == C2::exp, Quantity<Unit<CommonRatio<C1, C2>, K2>> >;

#if UNITS_HAS_ANY()
    template <typename K2>
    using IsAnyKind
        = std::is_same< typename K2::kind, Any >;

    template <typename C1, typename K1, typename C2, typename K2, typename T = int>
    using EnableIfConvertibleFromAny
        = std::enable_if_t< Divides<C1, C2>::value && IsSameDimensionType<K1, K2>::value && IsAnyKind<K2>::value, T >;
#endif

    template <typename K2, typename T = int>
    using EnableIfDimensionless
        = std::enable_if_t< std::is_same< typename K2::dimension, Dimension<> >::value, T >;

public:
    constexpr Quantity() noexcept = default;
    constexpr Quantity(const Quantity&) noexcept = default;
    constexpr Quantity& operator=(const Quantity&) noexcept = default;

    constexpr explicit Quantity(double c) noexcept
        : count_(c)
    {
    }

    template <typename C2, EnableIfDivides<C, C2> = 0>
    constexpr Quantity(Quantity<Unit<C2, K>> q) noexcept
        : count_(DivConversions<C2, C>{}(q.count()))
    {
    }

#if UNITS_HAS_ANY()
    template <typename C2, typename K2, EnableIfConvertibleFromAny<C, K, C2, K2> = 0>
    constexpr Quantity(Quantity<Unit<C2, K2>> q) noexcept
        : count_(DivConversions<C2, C>{}(q.count()))
    {
    }
#endif

    template <typename C2, typename K2, EnableIfCompatible<K, K2> = 0>
    constexpr explicit Quantity(Quantity<Unit<C2, K2>> q) noexcept
        : count_(DivConversions<C2, C>{}(q.count()))
    {
    }

#if UNITS_DELETE_EVERYTHING_ELSE()
    template <typename U2>
    constexpr Quantity(Quantity<U2> q) noexcept
        = delete;
#endif

    [[nodiscard]] constexpr double count() const noexcept {
        return count_;
    }

#if 1 // DANGER!!!
    [[nodiscard]] constexpr double value() const noexcept {
        return conversion{}(count_);
    }
#endif

    template <typename C2, typename K2, EnableIfCompatible<K, K2> = 0>
    [[nodiscard]] constexpr auto convert_to(Unit<C2, K2>) const noexcept {
        return Quantity<Unit<C2, K>>(DivConversions<C, C2>{}(count()));
    }

    template <typename C2, typename K2, EnableIfCompatible<K, K2> = 0>
    [[nodiscard]] constexpr auto convert_to(Quantity<Unit<C2, K2>>) const noexcept {
        return Quantity<Unit<C2, K>>(DivConversions<C, C2>{}(count()));
    }

#if UNITS_DELETE_EVERYTHING_ELSE()
    template <typename U2>
    constexpr auto convert_to(Quantity<U2>) const noexcept
        = delete;
#endif

    //------------------------------------------------------------------------------
    // Arithmetic

    [[nodiscard]] constexpr friend auto operator+(Quantity q) noexcept {
        return q;
    }

    [[nodiscard]] constexpr friend auto operator-(Quantity q) noexcept {
        return Quantity(-q.count());
    }

    template <typename C2, typename Q = CommonTypeSfinae<C, C2, K>>
    [[nodiscard]] constexpr friend auto operator+(Quantity lhs, Quantity<Unit<C2, K>> rhs) noexcept {
        //static_assert(std::is_convertible<decltype(lhs), Q>::value, "");
        //static_assert(std::is_convertible<decltype(rhs), Q>::value, "");
        return Q(Q(lhs).count() + Q(rhs).count());
    }

    template <typename C2, typename Q = CommonTypeSfinae<C, C2, K>>
    [[nodiscard]] constexpr friend auto operator-(Quantity lhs, Quantity<Unit<C2, K>> rhs) noexcept {
        //static_assert(std::is_convertible<decltype(lhs), Q>::value, "");
        //static_assert(std::is_convertible<decltype(rhs), Q>::value, "");
        return Q(Q(lhs).count() - Q(rhs).count());
    }

    template <typename U2>
    [[nodiscard]] constexpr friend auto operator*(Quantity lhs, Quantity<U2> rhs) noexcept {
        using Q = Quantity<MulUnits<unit, U2>>;
        return Q(lhs.count() * rhs.count());
    }

    template <typename U2>
    [[nodiscard]] constexpr friend auto operator/(Quantity lhs, Quantity<U2> rhs) noexcept {
        using Q = Quantity<DivUnits<unit, U2>>;
        return Q(lhs.count() / rhs.count());
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
        using Q = Quantity<DivUnits<Unit<Ratio<1>, kinds::One>, unit>>;
        return Q(lhs / rhs.count());
    }

#if 0 // UNITS_DELETE_EVERYTHING_ELSE()
    template <typename U2>
    constexpr friend void operator+(Quantity lhs, Quantity<U2> rhs) noexcept
        = delete;

    template <typename U2>
    constexpr friend void operator-(Quantity lhs, Quantity<U2> rhs) noexcept
        = delete;
#endif

    //------------------------------------------------------------------------------
    // Assignment operators

    template <typename C2, EnableIfDivides<C, C2> = 0>
    constexpr friend Quantity& operator+=(Quantity& lhs, Quantity<Unit<C2, K>> rhs) noexcept {
        //static_assert(std::is_convertible<decltype(rhs), Quantity>::value, "");
//      lhs.count_ += Quantity(rhs).count();
//      lhs.count_ += DivConversions<C2, C>{}(rhs.count());
        lhs.count_ += DivConversions<C2, C>::num * rhs.count();
        return lhs;
    }

    template <typename C2, EnableIfDivides<C, C2> = 0>
    constexpr friend Quantity& operator-=(Quantity& lhs, Quantity<Unit<C2, K>> rhs) noexcept {
        //static_assert(std::is_convertible<decltype(rhs), Quantity>::value, "");
//      lhs.count_ -= Quantity(rhs).count();
//      lhs.count_ -= DivConversions<C2, C>{}(rhs.count());
        lhs.count_ -= DivConversions<C2, C>::num * rhs.count();
        return lhs;
    }

    constexpr friend Quantity& operator*=(Quantity& lhs, double rhs) noexcept {
        return lhs = lhs * rhs;
    }

    constexpr friend Quantity& operator/=(Quantity& lhs, double rhs) noexcept {
        return lhs = lhs / rhs;
    }

#if UNITS_DELETE_EVERYTHING_ELSE()
    template <typename U2>
    constexpr friend void operator+=(Quantity& lhs, Quantity<U2> rhs) noexcept
        = delete;

    template <typename U2>
    constexpr friend void operator-=(Quantity& lhs, Quantity<U2> rhs) noexcept
        = delete;
#endif

    //------------------------------------------------------------------------------
    // Comparisons

    template <typename C2, typename Q = CommonTypeSfinae<C, C2, K>>
    [[nodiscard]] constexpr friend int Compare(Quantity lhs, Quantity<Unit<C2, K>> rhs) noexcept {
        //static_assert(std::is_convertible<decltype(lhs), Q>::value, "");
        //static_assert(std::is_convertible<decltype(rhs), Q>::value, "");
        return impl::CompareValues(Q(lhs).count(), Q(rhs).count());
    }

#if UNITS_DELETE_EVERYTHING_ELSE()
    template <typename U2>
    constexpr friend void Compare(Quantity lhs, Quantity<U2> rhs) noexcept
        = delete;
#endif

    // lhs == rhs <==> lhs - rhs == 0
    // lhs >= rhs <==> lhs - rhs >= 0
    // etc...
    // If you can add them, you can compare them...

    template <typename C2, typename Q = CommonTypeSfinae<C, C2, K>>
    [[nodiscard]] constexpr friend bool operator==(Quantity lhs, Quantity<Unit<C2, K>> rhs) noexcept {
        //static_assert(std::is_convertible<decltype(lhs), Q>::value, "");
        //static_assert(std::is_convertible<decltype(rhs), Q>::value, "");
        return Q(lhs).count() == Q(rhs).count();
    }

    template <typename C2, typename Q = CommonTypeSfinae<C, C2, K>>
    [[nodiscard]] constexpr friend bool operator!=(Quantity lhs, Quantity<Unit<C2, K>> rhs) noexcept {
        return !(lhs == rhs);
    }

    template <typename C2, typename Q = CommonTypeSfinae<C, C2, K>>
    [[nodiscard]] constexpr friend bool operator<(Quantity lhs, Quantity<Unit<C2, K>> rhs) noexcept {
        //static_assert(std::is_convertible<decltype(lhs), Q>::value, "");
        //static_assert(std::is_convertible<decltype(rhs), Q>::value, "");
        return Q(lhs).count() < Q(rhs).count();
    }

    template <typename C2, typename Q = CommonTypeSfinae<C, C2, K>>
    [[nodiscard]] constexpr friend bool operator>(Quantity lhs, Quantity<Unit<C2, K>> rhs) noexcept {
        return rhs < lhs;
    }

    template <typename C2, typename Q = CommonTypeSfinae<C, C2, K>>
    [[nodiscard]] constexpr friend bool operator<=(Quantity lhs, Quantity<Unit<C2, K>> rhs) noexcept {
        return !(rhs < lhs);
    }

    template <typename C2, typename Q = CommonTypeSfinae<C, C2, K>>
    [[nodiscard]] constexpr friend bool operator>=(Quantity lhs, Quantity<Unit<C2, K>> rhs) noexcept {
        return !(lhs < rhs);
    }

#if UNITS_DELETE_EVERYTHING_ELSE()
    template <typename U2>
    constexpr friend void operator==(Quantity lhs, Quantity<U2> rhs) noexcept
        = delete;

    template <typename U2>
    constexpr friend void operator!=(Quantity lhs, Quantity<U2> rhs) noexcept
        = delete;

    template <typename U2>
    constexpr friend void operator<(Quantity lhs, Quantity<U2> rhs) noexcept
        = delete;

    template <typename U2>
    constexpr friend void operator>(Quantity lhs, Quantity<U2> rhs) noexcept
        = delete;

    template <typename U2>
    constexpr friend void operator<=(Quantity lhs, Quantity<U2> rhs) noexcept
        = delete;

    template <typename U2>
    constexpr friend void operator>=(Quantity lhs, Quantity<U2> rhs) noexcept
        = delete;
#endif

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

    template <typename K2 = K, EnableIfDimensionless<K2> = 0>
    friend auto Exp(Quantity q) noexcept {
        return Quantity(std::exp(q.count()));
    }

    // log_a(t) = log_b(t) / log_b(a)
    //  ... !?!?

    template <typename K2 = K, EnableIfDimensionless<K2> = 0>
    friend auto Log(Quantity q) noexcept {
        return Quantity(std::log(q.count()));
    }

    template <typename C1>
    using SqrtTypeSfinae = void;

    template <typename Q = SqrtTypeSfinae<C>>
    friend auto Sqrt(Quantity q) noexcept {
    }

#endif
};

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

#if UNITS_HAS_ANY()
template <typename C, typename K>
[[nodiscard]] constexpr auto as_any_quantity(Quantity<Unit<C, K>> q) noexcept {
    return Quantity<Unit<C, Kind<kinds::Any, typename K::dimension>>>(q.count());
}

template <typename C, typename K>
[[nodiscard]] constexpr auto flatten(Quantity<Unit<C, K>> q) noexcept {
    return as_any_quantity(q);
}
#endif

template <typename C, typename K>
[[nodiscard]] constexpr auto remove_conversion(Quantity<Unit<C, K>> q) noexcept {
    return Quantity<Unit<Ratio<1>, K>>(C{}(q.count()));
}

//==================================================================================================
// Typedefs
//==================================================================================================

//namespace prefixes
//{
//    using Exa   = Ratio<1000000000000000000>;
//    using Peta  = Ratio<1000000000000000>;
//    using Tera  = Ratio<1000000000000>;
//    using Giga  = Ratio<1000000000>;
//    using Mega  = Ratio<1000000>;
//    using Kilo  = Ratio<1000>;
//    using Hecto = Ratio<100>;
//    using Deca  = Ratio<10>;
//    using One   = Ratio<1>;
//    using Deci  = Ratio<1, 10>;
//    using Centi = Ratio<1, 100>;
//    using Milli = Ratio<1, 1000>;
//    using Micro = Ratio<1, 1000000>;
//    using Nano  = Ratio<1, 1000000000>;
//    using Pico  = Ratio<1, 1000000000000>;
//    using Femto = Ratio<1, 1000000000000000>;
//    using Atto  = Ratio<1, 1000000000000000000>;
//}

//--------------------------------------------------------------------------------------------------
// Length

namespace units
{
    using Nanometre   = Unit< Ratio<1, 1000000000>, kinds::Length >;
    using Micrometre  = Unit< Ratio<1, 1000000>, kinds::Length >;
    using Millimetre  = Unit< Ratio<1, 1000>, kinds::Length >;
    using Centimetre  = Unit< Ratio<1, 100>, kinds::Length >;
    using Decimetre   = Unit< Ratio<1, 10>, kinds::Length >;
    using Metre       = Unit< Ratio<1>, kinds::Length >;
    using Hectometre  = Unit< Ratio<100>, kinds::Length >;
    using Kilometre   = Unit< Ratio<1000>, kinds::Length >;
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
#if 1
    using Inch = ScaledUnit< Ratio<254, 100>, Centimetre >;
    using Foot = ScaledUnit< Ratio<12>, Inch >;
    using Yard = ScaledUnit< Ratio<3>, Foot >;
    using Mile = ScaledUnit< Ratio<1760>, Yard >;
#else
    using Inch = Unit< MulConversions<Ratio<254, 100>, Centimetre::conversion>, kinds::Length >; // (international)
    using Foot = Unit< MulConversions<Ratio<12>, Inch::conversion>, kinds::Length >;             // (international)
    using Yard = Unit< MulConversions<Ratio<3>, Foot::conversion>, kinds::Length >;              // (international)
    using Mile = Unit< MulConversions<Ratio<1760>, Yard::conversion>, kinds::Length >;           // (international)
#endif
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
    using SquareCentimetre = Unit< Square<Centimetre::conversion>, kinds::Area >;
    using SquareDecimetre  = Unit< Square<Decimetre::conversion>, kinds::Area >;
    using SquareMetre      = Unit< Square<Metre::conversion>, kinds::Area >;
    using SquareKilometre  = Unit< Square<Kilometre::conversion>, kinds::Area >;

    using Hectare = Unit< Square<Hectometre::conversion>, kinds::Area >;
}

using SquareCentimetres = Quantity< units::SquareCentimetre >;
using SquareDecimetres  = Quantity< units::SquareDecimetre >;
using SquareMetres      = Quantity< units::SquareMetre >;
using SquareKilometres  = Quantity< units::SquareKilometre >;

//--------------------------------------------------------------------------------------------------
// Volume

namespace units
{
    using CubicCentimetre = Unit< Cubic<Centimetre::conversion>, kinds::Volume >;
    using CubicDecimetre  = Unit< Cubic<Decimetre::conversion>, kinds::Volume >;
    using CubicMetre      = Unit< Cubic<Metre::conversion>, kinds::Volume >;

    using Litre = CubicDecimetre;
}

//--------------------------------------------------------------------------------------------------
// Time

namespace units
{
    using Second      = Unit< Ratio<1>, kinds::Time >;
    using Millisecond = ScaledUnit< Ratio<1, 1000>, Second >;
    using Microsecond = ScaledUnit< Ratio<1, 1000>, Millisecond >;
    using Nanosecond  = ScaledUnit< Ratio<1, 1000>, Microsecond >;
#if 1
    using Minute      = ScaledUnit< Ratio<60>, Second >;
    using Hour        = ScaledUnit< Ratio<60>, Minute >;
    //using Day         = ScaledUnit< Ratio<24>, Hour >;
    //using Week        = ScaledUnit< Ratio<7>, Day >;
    //using Year        = ScaledUnit< Ratio<146097, 400>, Day >;
    //using Month       = ScaledUnit< Ratio<1, 12>, Year >;
#else
    using Minute      = Unit< MulConversions<Ratio<60>, Second::conversion>, kinds::Time >;
    using Hour        = Unit< MulConversions<Ratio<60>, Minute::conversion>, kinds::Time >;
    //using Day         = Unit< MulConversions<Ratio<24>, Hour::conversion>, kinds::Time >;
    //using Week        = Unit< MulConversions<Ratio<7>, Day::conversion>, kinds::Time >;
    //using Year        = Unit< MulConversions<Ratio<146097, 400>, Day::conversion>, kinds::Time >;
    //using Month       = Unit< MulConversions<Ratio<1, 12>, Year::conversion>, kinds::Time >;
#endif
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
    using Hertz = Unit< DivConversions<Ratio<1>, Second::conversion>, kinds::Frequency >;
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
    using MetrePerSecond
        = Unit< DivConversions<Metre::conversion, Second::conversion>,
                DivKinds<kinds::Length, kinds::Time> >;
    using KilometrePerHour
        = Unit< DivConversions<Kilometre::conversion, Hour::conversion>,
                DivKinds<kinds::Length, kinds::Time> >;
    using MilePerHour
        = Unit< DivConversions<Mile::conversion, Hour::conversion>,
                DivKinds<kinds::Length, kinds::Time> >;
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
//  using MetrePerSecondSquared = Unit< MulConversions<MetrePerSecond::conversion, Second::conversion>, kinds::Acceleration >;
    using MetrePerSecondSquared = Unit< DivConversions<Metre::conversion, Square<Second::conversion>>, kinds::Acceleration >;
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
    using Gram     = Unit< Ratio<1, 1000>, kinds::Mass >;
    using Kilogram = Unit< Ratio<1>, kinds::Mass >;
    using Tonne    = Unit< Ratio<1000>, kinds::Mass >;
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
    using Newton = Unit< MulConversions<Kilogram::conversion, MetrePerSecondSquared::conversion>, kinds::Force >;
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

namespace units
{
    using Kelvin  = Unit< Ratio<1>, kinds::Temperature >;
//  using Celsius = Unit< Ratio<1>, kinds::Temperature >;
}

using Kelvin  = Quantity< units::Kelvin >;
//using Celsius = Quantity< units::Celsius >;

namespace literals
{
    [[nodiscard]] constexpr auto operator""_K(long double x) noexcept {
        return Kelvin{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_K(unsigned long long x) noexcept {
        return Kelvin{static_cast<double>(x)};
    }
//  [[nodiscard]] constexpr auto operator""_degC(long double x) noexcept {
//      return Celsius{static_cast<double>(x)};
//  }
//  [[nodiscard]] constexpr auto operator""_degC(unsigned long long x) noexcept {
//      return Celsius{static_cast<double>(x)};
//  }
}

//--------------------------------------------------------------------------------------------------
// Amount of substance

namespace units
{
    using Mole = Unit< Ratio<1>, kinds::AmountOfSubstance >;
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
    using Radian     = Unit< Ratio<1>, kinds::PlaneAngle >;
    using Degree     = Unit< Ratio<1, 180, /* pi^ */ 1>, kinds::PlaneAngle >;
    using Gon        = Unit< Ratio<1, 200, /* pi^ */ 1>, kinds::PlaneAngle >;
    using Revolution = Unit< Ratio<2,   1, /* pi^ */ 1>, kinds::PlaneAngle >;
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
    using Steradian    = Unit< Ratio<1>, kinds::SolidAngle >;
    using SquareDegree = Unit< Square<Degree::conversion>, kinds::SolidAngle >; // sq.deg = deg^2
}

using Steradians    = Quantity< units::Steradian >;
//using SquareDegrees = Quantity< units::SquareDegree >;

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
    using Bit      = Unit< Ratio<1>, kinds::Bit >;
    using Nibble   = ScaledUnit< Ratio<4>, Bit >;
    using Byte     = ScaledUnit< Ratio<8>, Bit >;
    using Kilobyte = ScaledUnit< Ratio<1000>, Byte >;
    using Megabyte = ScaledUnit< Ratio<1000>, Kilobyte >;
    using Gigabyte = ScaledUnit< Ratio<1000>, Megabyte >;
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
    using Candela
        = Unit< Ratio<1>, kinds::LuminousIntensity >;

    using Lumen
        = Unit< MulConversions<Candela::conversion, Steradian::conversion>, kinds::LuminousFlux >;

//  using LumenSecond
//      = Unit< MulConversions<Lumen::conversion, Second::conversion>, kinds::LuminousEnergy >;
    using Talbot
        = Unit< MulConversions<Lumen::conversion, Second::conversion>, kinds::LuminousEnergy >;

//  using CandelaPerSquareMeter
//      = Unit< DivConversions<Candela::conversion, SquareMeter::conversion>, kinds::Luminance >;
    using Nit
        = Unit< DivConversions<Candela::conversion, SquareMetre::conversion>, kinds::Luminance >;

    using Lux
        = Unit< DivConversions<Lumen::conversion, SquareMetre::conversion>, kinds::Illuminance >;
}

using Candelas = Quantity< units::Candela >;
//using Talbots  = Quantity< units::Talbot >;
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
