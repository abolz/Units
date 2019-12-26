// Copyright Alexander Bolz 2019
// Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at https://www.boost.org/LICENSE_1_0.txt)

#pragma once

#define UNITS_DIMENSIONLESS_ARITHMETIC() 0
#define UNITS_DELETE_EVERYTHING_ELSE() 0
#define UNITS_HAS_ANY() 0
#define UNITS_HAS_MATH() 0
#define UNITS_SI_DIMENSION() 0

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

template <typename K, typename D>
struct Kind
{
    using type = Kind;
    using kind = K;
    using dimension = D;
};

template <Natural Num, Natural Den, int Exp>
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

// SI base units
//    Base quantity:                  Name:
//      time                            second (s)
//      length                          metre (m)
//      mass                            kilogram (kg)
//      electric current                ampere (A)
//      thermodynamic temperature       kelvin (K)
//      amount of substance             mole (mol)
//      luminous intensity              candela (cd)

template <typename... BaseQuantities>
struct Dimension
{
};

namespace impl
{
    template <template <int> class D, int E>
    constexpr auto NegateExp(D<E>) noexcept {
        return D<-E>{};
    }

    template <typename... Un>
    constexpr auto Invert(Dimension<Un...>)
    {
        return Dimension<decltype(NegateExp(Un{}))...>{};
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
            if constexpr (ExpL1 + ExpR1 != 0)
                return Concat(Dimension<L1<ExpL1 + ExpR1>>{}, Merge(Dimension<Ln...>{}, Dimension<Rn...>{}));
            else
                return Merge(Dimension<Ln...>{}, Dimension<Rn...>{});
        }
    }
}

template <typename D1, typename D2>
using MulDimensions = decltype(impl::Merge(D1{}, D2{}));

template <typename D1, typename D2>
using DivDimensions = decltype(impl::Merge(D1{}, impl::Invert(D2{})));

namespace dim // Base quantities
{
    template <int E> struct Length            { static constexpr int64_t id = 1; };
    template <int E> struct Mass              { static constexpr int64_t id = 2; };
    template <int E> struct Time              { static constexpr int64_t id = 3; };
    template <int E> struct ElectricCurrent   { static constexpr int64_t id = 4; };
    template <int E> struct Temperature       { static constexpr int64_t id = 5; };
    template <int E> struct AmountOfSubstance { static constexpr int64_t id = 6; };
    template <int E> struct LuminousIntensity { static constexpr int64_t id = 7; };
    template <int E> struct Bit               { static constexpr int64_t id = 8; };
    template <int E> struct Currency          { static constexpr int64_t id = 9; };
}

namespace kinds
{
    struct Any {};

    struct One :                    Kind< One,                  Dimension<> > {};

    struct Length :                 Kind< Length,               Dimension<dim::Length<1>> > {};
    struct Mass :                   Kind< Mass,                 Dimension<dim::Mass<1>> > {};
    struct Time :                   Kind< Time,                 Dimension<dim::Time<1>> > {};
    struct ElectricCurrent :        Kind< ElectricCurrent,      Dimension<dim::ElectricCurrent<1>> > {};
    struct Temperature :            Kind< Temperature,          Dimension<dim::Temperature<1>> > {};
    struct AmountOfSubstance :      Kind< AmountOfSubstance,    Dimension<dim::AmountOfSubstance<1>> > {};
    struct LuminousIntensity :      Kind< LuminousIntensity,    Dimension<dim::LuminousIntensity<1>> > {};

    struct Area :                   Kind< Area,                 Dimension<dim::Length< 2>                              > > {}; // m^2
    struct Volume :                 Kind< Volume,               Dimension<dim::Length< 3>                              > > {}; // m^3
    struct PlaneAngle :             Kind< PlaneAngle,           Dimension<                                             > > {}; // 1
    struct SolidAngle :             Kind< SolidAngle,           Dimension<                                             > > {}; // 1
    struct Velocity :               Kind< Velocity,             Dimension<dim::Length< 1>,                dim::Time<-1>> > {}; // m s^-1
    struct Acceleration :           Kind< Acceleration,         Dimension<dim::Length< 1>,                dim::Time<-2>> > {}; // m s^-2
    struct Frequency :              Kind< Frequency,            Dimension<                                dim::Time<-1>> > {}; // s^-1
    struct Density :                Kind< Density,              Dimension<dim::Length<-3>, dim::Mass< 1>               > > {}; // kg m^-3
    struct SurfaceDensity :         Kind< SurfaceDensity,       Dimension<dim::Length<-2>, dim::Mass< 1>               > > {}; // kg m^-2
    struct Impulse :                Kind< Impulse,              Dimension<dim::Length< 1>, dim::Mass< 1>, dim::Time<-1>> > {}; // kg m s^-1
    struct Force :                  Kind< Force,                Dimension<dim::Length< 1>, dim::Mass< 1>, dim::Time<-2>> > {}; // kg m s^-2
    struct Energy :                 Kind< Energy,               Dimension<dim::Length< 2>, dim::Mass< 1>, dim::Time<-2>> > {}; // kg m^2 s^-2
    struct Torque :                 Kind< Torque,               Dimension<dim::Length< 2>, dim::Mass< 1>, dim::Time<-2>> > {}; // kg m^2 s^-2
    struct Power :                  Kind< Power,                Dimension<dim::Length< 2>, dim::Mass< 1>, dim::Time<-3>> > {}; // kg m^2 s^-3
    struct Pressure :               Kind< Pressure,             Dimension<dim::Length<-1>, dim::Mass< 1>, dim::Time<-2>> > {}; // kg m^-1 s^-2
    struct AngularVelocity :        Kind< AngularVelocity,      Dimension<                                dim::Time<-1>> > {}; // s^-1
    struct AngularAcceleration :    Kind< AngularAcceleration,  Dimension<                                dim::Time<-2>> > {}; // s^-2

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
}

template <typename K1, typename K2>
using MulKinds = typename kinds::Product<K1, K2>::type;

template <typename K1, typename K2>
using DivKinds = typename kinds::Quotient<K1, K2>::type;

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

    constexpr Natural Abs(Natural x) noexcept {
        return x < 0 ? -x : x;
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

#if 0
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

template <Natural Num, Natural Den, Exponent Exp>
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

    [[nodiscard]] constexpr double operator()(double x) const noexcept {
        return ApplyExp(ApplyRat(x));
    }

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
        constexpr double Powers[] = {
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

namespace impl
{
    template <typename R1, typename R2>
    struct RationalProduct
    {
        static constexpr Natural Num1 = R1::num;
        static constexpr Natural Den1 = R1::den;
        static constexpr Natural Num2 = R2::num;
        static constexpr Natural Den2 = R2::den;

        static constexpr Natural S = Gcd(Num1, Den2);
        static constexpr Natural T = Gcd(Den1, Num2);

        using type = Rational< (Num1 / S) * (Num2 / T), (Den1 / T) * (Den2 / S), R1::exp + R2::exp >;
    };
}

template <typename C1, typename C2>
using MulRatios = typename impl::RationalProduct<C1, C2>::type;

template <typename C1, typename C2>
using DivRatios = typename impl::RationalProduct<C1, Ratio<C2::den, C2::num, -C2::exp>>::type;

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
};

template <typename U1, typename U2>
using MulUnits
    = Unit< MulRatios<typename U1::conversion, typename U2::conversion>, MulKinds<typename U1::kind, typename U2::kind> >;

template <typename U1, typename U2>
using DivUnits
    = Unit< DivRatios<typename U1::conversion, typename U2::conversion>, DivKinds<typename U1::kind, typename U2::kind> >;

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
    template <typename C2>
    using IsNatural
        = std::bool_constant< C2::den == 1 && C2::exp == 0 >;

    // C1 | C2
    template <typename C1, typename C2>
    using Divides
        = IsNatural< DivRatios<C2, C1> >;

    template <typename C1, typename C2, typename T = int>
    using EnableIfDivides
        = std::enable_if_t< Divides<C1, C2>::value, T >;

    template <typename K1, typename K2>
    using IsSameDimensionType
        = std::is_same<typename K1::dimension, typename K2::dimension>;

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
        = std::is_same< typename K2::kind, kinds::Any >;

    template <typename C1, typename K1, typename C2, typename K2, typename T = int>
    using EnableIfConvertibleFromAny
        = std::enable_if_t< Divides<C1, C2>::value && IsSameDimensionType<K1, K2>::value && IsAnyKind<K2>::value, T >;
#endif

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
        : count_(DivRatios<C2, C>{}(q.count()))
    {
    }

#if UNITS_HAS_ANY()
    template <typename C2, typename K2, EnableIfConvertibleFromAny<C, K, C2, K2> = 0>
    constexpr Quantity(Quantity<Unit<C2, K2>> q) noexcept
        : count_(DivRatios<C2, C>{}(q.count()))
    {
    }
#endif

    template <typename C2, typename K2, EnableIfCompatible<K, K2> = 0>
    constexpr explicit Quantity(Quantity<Unit<C2, K2>> q) noexcept
        : count_(DivRatios<C2, C>{}(q.count()))
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
    [[nodiscard]] constexpr auto convert_to(Quantity<Unit<C2, K2>>) const noexcept {
        return Quantity<Unit<C2, K>>(DivRatios<C, C2>{}(count()));
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
        return Q(Q(lhs).count() + Q(rhs).count());
    }

    template <typename C2, typename Q = CommonTypeSfinae<C, C2, K>>
    [[nodiscard]] constexpr friend auto operator-(Quantity lhs, Quantity<Unit<C2, K>> rhs) noexcept {
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
        using U = Unit<Ratio<1>, kinds::One>;
        using Q = Quantity<DivUnits<U, unit>>;
        return Q(lhs / rhs.count());
    }

#if 0//UNITS_DELETE_EVERYTHING_ELSE()
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
        //lhs.count_ += Quantity(rhs).count();
        //lhs.count_ += DivRatios<C2, C>{}(rhs.count());
        lhs.count_ += DivRatios<C2, C>::num * rhs.count();
        return lhs;
    }

    template <typename C2, EnableIfDivides<C, C2> = 0>
    constexpr friend Quantity& operator-=(Quantity& lhs, Quantity<Unit<C2, K>> rhs) noexcept {
        //lhs.count_ -= Quantity(rhs).count();
        //lhs.count_ -= DivRatios<C2, C>{}(rhs.count());
        lhs.count_ -= DivRatios<C2, C>::num * rhs.count();
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
        return Q(lhs).count() == Q(rhs).count();
    }

    template <typename C2, typename Q = CommonTypeSfinae<C, C2, K>>
    [[nodiscard]] constexpr friend bool operator!=(Quantity lhs, Quantity<Unit<C2, K>> rhs) noexcept {
        return !(lhs == rhs);
    }

    template <typename C2, typename Q = CommonTypeSfinae<C, C2, K>>
    [[nodiscard]] constexpr friend bool operator<(Quantity lhs, Quantity<Unit<C2, K>> rhs) noexcept {
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
    constexpr friend auto Fma(Quantity a, Quantity<U2> b, decltype(a * b) c) noexcept
        -> decltype(a * b)
    {
        //
        // NB:
        //
        // This function only allows arguments as the third parameter, which are implicitly
        // convertible to the exact type of the multiplication (a * b). This is more restrictive
        // than the expression (a * b + c). But otherwise, we would need to convert the
        // product (a * b) and the third argument (c) to their common type before the addition,
        // which would defeat the application of std::fma.
        //

        using Q = decltype(a * b);
        return Q(a.count() * b.count() + c.count());
//      return Q(std::fma(a.count(), b.count(), c.count()));
    }

#endif
};

template <typename T, typename U> // T = Quantity (or Unit)
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

//==================================================================================================
// Typedefs
//==================================================================================================

//namespace prefixes
//{
////  using Exa   = Ratio<1000000000000000000>;
////  using Peta  = Ratio<1000000000000000>;
////  using Tera  = Ratio<1000000000000>;
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
////  using Pico  = Ratio<1, 1000000000000>;
////  using Femto = Ratio<1, 1000000000000000>;
////  using Atto  = Ratio<1, 1000000000000000000>;
//}

//--------------------------------------------------------------------------------------------------
// Length

namespace units
{
#if 0
    using Nanometer   = Unit< Ratio<1, 1000000>, kinds::Length >;
    using Micrometer  = Unit< Ratio<1, 1000>, kinds::Length >;
    using Millimeter  = Unit< Ratio<1>, kinds::Length >;
    using Centimeter  = Unit< Ratio<10>, kinds::Length >;
    using Decimeter   = Unit< Ratio<100>, kinds::Length >;
    using Meter       = Unit< Ratio<1000>, kinds::Length >;
    using Kilometer   = Unit< Ratio<1000000>, kinds::Length >;
#else
    using Nanometer   = Unit< Ratio<1, 1000000000>, kinds::Length >;
    using Micrometer  = Unit< Ratio<1, 1000000>, kinds::Length >;
    using Millimeter  = Unit< Ratio<1, 1000>, kinds::Length >;
    using Centimeter  = Unit< Ratio<1, 100>, kinds::Length >;
    using Decimeter   = Unit< Ratio<1, 10>, kinds::Length >;
    using Meter       = Unit< Ratio<1>, kinds::Length >;
    using Kilometer   = Unit< Ratio<1000>, kinds::Length >;
#endif
}

using Nanometers  = Quantity< units::Nanometer >;
using Micrometers = Quantity< units::Micrometer >;
using Millimeters = Quantity< units::Millimeter >;
using Centimeters = Quantity< units::Centimeter >;
using Decimeters  = Quantity< units::Decimeter >;
using Meters      = Quantity< units::Meter >;
using Kilometers  = Quantity< units::Kilometer >;

namespace literals
{
    [[nodiscard]] constexpr auto operator""_nm(long double x) noexcept {
        return Nanometers{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_nm(unsigned long long x) noexcept {
        return Nanometers{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_mm(long double x) noexcept {
        return Millimeters{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_mm(unsigned long long x) noexcept {
        return Millimeters{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_cm(long double x) noexcept {
        return Centimeters{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_cm(unsigned long long x) noexcept {
        return Centimeters{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_dm(long double x) noexcept {
        return Decimeters{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_dm(unsigned long long x) noexcept {
        return Decimeters{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_m(long double x) noexcept {
        return Meters{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_m(unsigned long long x) noexcept {
        return Meters{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_km(long double x) noexcept {
        return Kilometers{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_km(unsigned long long x) noexcept {
        return Kilometers{static_cast<double>(x)};
    }
}

namespace units
{
    using Inch = Unit< MulRatios<Ratio<254, 100>, Centimeter::conversion>, kinds::Length >; // (international)
    using Foot = Unit< MulRatios<Ratio<12>, Inch::conversion>, kinds::Length >;             // (international)
    using Yard = Unit< MulRatios<Ratio<3>, Foot::conversion>, kinds::Length >;              // (international)
    using Mile = Unit< MulRatios<Ratio<1760>, Yard::conversion>, kinds::Length >;           // (international)
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
// Time

namespace units
{
    using Nanosecond  = Unit< Ratio<1, 1000000000>, kinds::Time >;
    using Microsecond = Unit< Ratio<1, 1000000>, kinds::Time >;
    using Millisecond = Unit< Ratio<1, 1000>, kinds::Time >;
    using Second      = Unit< Ratio<1>, kinds::Time >;
    using Minute      = Unit< MulRatios<Ratio<60>, Second::conversion>, kinds::Time >;
    using Hour        = Unit< MulRatios<Ratio<60>, Minute::conversion>, kinds::Time >;
    // using Day         = Unit< MulRatios<Ratio<24>, Hour::conversion>, kinds::Time >;
    // using Week        = Unit< MulRatios<Ratio<7>, Day::conversion>, kinds::Time >;
    // using Year        = Unit< MulRatios<Ratio<146097, 400>, Day::conversion>, kinds::Time >;
    // using Month       = Unit< MulRatios<Ratio<1, 12>, Year::conversion>, kinds::Time >;
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
// Velocity

namespace units
{
    using MeterPerSecond
        = Unit< DivRatios<Meter::conversion, Second::conversion>,
                //DivKinds<kinds::Length, kinds::Time> >;
                kinds::Velocity >;
    using KilometerPerHour
        = Unit< DivRatios<Kilometer::conversion, Hour::conversion>,
                //DivKinds<kinds::Length, kinds::Time> >;
                kinds::Velocity >;
    using MilePerHour
        = Unit< DivRatios<Mile::conversion, Hour::conversion>,
                //DivKinds<kinds::Length, kinds::Time> >;
                kinds::Velocity >;
}

using MetersPerSecond   = Quantity< units::MeterPerSecond >;
using KilometersPerHour = Quantity< units::KilometerPerHour >;
using MilesPerHour      = Quantity< units::MilePerHour >;

namespace literals
{
    [[nodiscard]] constexpr auto operator""_mps(long double x) noexcept {
        return MetersPerSecond{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_mps(unsigned long long x) noexcept {
        return MetersPerSecond{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_kmph(long double x) noexcept {
        return KilometersPerHour{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_kmph(unsigned long long x) noexcept {
        return KilometersPerHour{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_miph(long double x) noexcept {
        return MilesPerHour{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_miph(unsigned long long x) noexcept {
        return MilesPerHour{static_cast<double>(x)};
    }
}

//--------------------------------------------------------------------------------------------------
// Mass

namespace units
{
    using Gram     = Unit< Ratio<1, 1000>, kinds::Mass >;
    using Kilogram = Unit< Ratio<1>, kinds::Mass >;
    using Ton      = Unit< Ratio<1000>, kinds::Mass >;
}

using Grams     = Quantity< units::Gram >;
using Kilograms = Quantity< units::Kilogram >;
using Tons      = Quantity< units::Ton >;

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
        return Tons{static_cast<double>(x)};
    }
    [[nodiscard]] constexpr auto operator""_t(unsigned long long x) noexcept {
        return Tons{static_cast<double>(x)};
    }
}

//--------------------------------------------------------------------------------------------------
// Plane angle

namespace units
{
    using Radian     = Unit< Ratio<1>, kinds::PlaneAngle >;
    using Degree     = Unit< Ratio<1, 180, /* pi^ */ 1>, kinds::PlaneAngle >;
    using Gon        = Unit< Ratio<1, 200, /* pi^ */ 1>, kinds::PlaneAngle >;
//  using Revolution = Unit< Ratio<2,   1, /* pi^ */ 1>, kinds::PlaneAngle >;
    using Revolution = Unit< MulRatios<Ratio<360>, Degree::conversion>, kinds::PlaneAngle >;
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
    using Steradian = Unit< Ratio<1>, kinds::SolidAngle >;
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

namespace kinds
{
    struct Bit : Kind< Bit, Dimension<dim::Bit<1>> > {};
}

namespace units
{
    using Bit      = Unit< Ratio<1>, kinds::Bit >;
    using Nibble   = Unit< MulRatios<Ratio<4>, Bit::conversion>, kinds::Bit >;
    using Byte     = Unit< MulRatios<Ratio<8>, Bit::conversion>, kinds::Bit >;
#if 1
    using Kilobyte = Unit< MulRatios<Ratio<1024>, Byte::conversion>, kinds::Bit >;
    using Megabyte = Unit< MulRatios<Ratio<1024>, Kilobyte::conversion>, kinds::Bit >;
    using Gigabyte = Unit< MulRatios<Ratio<1024>, Megabyte::conversion>, kinds::Bit >;
#else
    using Kilobyte = Unit< MulRatios<Ratio<1000>, Byte::conversion>, kinds::Bit >;
    using Megabyte = Unit< MulRatios<Ratio<1000>, Kilobyte::conversion>, kinds::Bit >;
    using Gigabyte = Unit< MulRatios<Ratio<1000>, Megabyte::conversion>, kinds::Bit >;
#endif
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

} // namespace sc
