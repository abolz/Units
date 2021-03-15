#include "Unit.h"
#include "UnitLiterals.h"
#include "UnitMath.h"

#include <type_traits>

using namespace uom;
using namespace uom::literals;

template <typename T> struct Incomplet;

template <typename L, typename R> using Assign = decltype(std::declval<L&>() = std::declval<R>());

template <typename L, typename R> using Add = decltype(std::declval<L>() + std::declval<R>());
template <typename L, typename R> using Sub = decltype(std::declval<L>() - std::declval<R>());
template <typename L, typename R> using Mul = decltype(std::declval<L>() * std::declval<R>());
template <typename L, typename R> using Div = decltype(std::declval<L>() / std::declval<R>());

template <typename L, typename R> using AssignAdd = decltype(std::declval<L&>() += std::declval<R>());
template <typename L, typename R> using AssignSub = decltype(std::declval<L&>() -= std::declval<R>());
template <typename L, typename R> using AssignMul = decltype(std::declval<L&>() *= std::declval<R>());
template <typename L, typename R> using AssignDiv = decltype(std::declval<L&>() /= std::declval<R>());

namespace details
{
    template <typename T>
    struct Void {
        using type = void;
    };

    template <typename T>
    using Void_t = typename Void<T>::type;

    template <typename AlwaysVoid, template <typename...> class Op, typename... Args>
    struct Detector {
        static constexpr bool value = false;
    };

    template <template <typename...> class Op, typename... Args>
    struct Detector<Void_t<Op<Args...>>, Op, Args...> {
        static constexpr bool value = true;
    };
}

template <template <typename...> class Op, typename... Args>
inline constexpr bool Compiles = details::Detector<void, Op, Args...>::value;

template <typename T1, typename T2>
inline constexpr bool IsSame = std::is_same_v<std::remove_const_t<std::remove_reference_t<T1>>,
                                              std::remove_const_t<std::remove_reference_t<T2>>>;

static constexpr void fun_mm(Millimetres mm)
{
}
static constexpr void fun_mm(Centimetres cm)
{
}
//static constexpr void fun_mm(Metres m)
//{
//}

static constexpr void fun_cm(Centimetres cm)
{
}
static constexpr void fun_cm(Metres m)
{
}

static constexpr void test()
{
    {
        constexpr auto t0 = 1_m + 1_m;
        static_assert(t0.count<Metres>() == 2.0);
        constexpr auto t1 = 1_cm + 1_cm;
        static_assert(t1.count<Millimetres>() == 20.0);
        constexpr auto t2 = 1_cm + 1_mm;
        static_assert(IsSame<Millimetres, decltype(t2)>);
        static_assert(t2.count_unsafe() == 11.0);
        constexpr auto t3 = 1_mm + 1_cm;
        static_assert(IsSame<Millimetres, decltype(t3)>);
        static_assert(t3.count_unsafe() == 11.0);
        static_assert(t2 == t3);
        constexpr auto t4 = 1_mm - 1_cm;
        static_assert(IsSame<Millimetres, decltype(t4)>);
        static_assert(t4.count_unsafe() == -9.0);
        static_assert(t4 < t3);
        constexpr auto t5 = 1_cm - 1_mm;
        static_assert(IsSame<Millimetres, decltype(t5)>);
        static_assert(t5.count_unsafe() == 9.0);
        static_assert(t5 < 1_cm);
    }
    {
        constexpr auto t0 = 1_m * 1_s;
        constexpr auto t1 = 1_s * 1_m;
        static_assert( std::is_same_v<decltype(t0), decltype(t1)>);
        constexpr auto ts = t0 + t1;
        constexpr auto tp = t0 * t1;
    }
    {
        fun_mm(1_mm);
        fun_mm(1_cm);
//      fun_mm(1_m);

//      fun_cm(1_mm);
    }
#if 0
    {
        constexpr auto t0 = 1_m / 1_m;
        constexpr auto t1 = 1_s / 1_s;
        static_assert(!IsSame<decltype(t0), decltype(t1)>);
        static_assert(!Compiles<Add, decltype(t0), decltype(t1)>);
        static_assert(!Compiles<Sub, decltype(t0), decltype(t1)>);
        static_assert( std::is_constructible_v<decltype(t0), decltype(t1)>);
        static_assert(!std::is_convertible_v<decltype(t1), decltype(t0)>);
        static_assert( std::is_constructible_v<decltype(t1), decltype(t0)>);
        static_assert(!std::is_convertible_v<decltype(t0), decltype(t1)>);
    }
#endif

    {
        //constexpr auto t0 = 1_m2;
        //constexpr auto t1 = 1_m * 1_m;
        //static_assert( std::is_same_v<decltype(t0), decltype(t1)>);
        //constexpr auto ts = t0 + t1;
        //constexpr auto tp = t0 * t1;
    }

#if 0
    {
        constexpr auto t0 = 1_m / 1_s; // 1_mps;
        constexpr auto t1 = 1_mi / 1_h; // 1_miph;
        static_assert( std::is_constructible_v<MetresPerSecond, decltype(t1)>);
        static_assert(!std::is_convertible_v<decltype(t1), MetresPerSecond>);
        constexpr auto t2 = MetresPerSecond(t1);
        constexpr auto t3 = 1_km / 1_h;
        //MetresPerSecond t4;
        //t4 = t3;
        constexpr MetresPerSecond t4(t3);
        constexpr KilometresPerHour t5(1_mps);
        constexpr auto t6 = 1_m / 1_s + MetresPerSecond(1_km / 1_h);
    }
#endif

    {
//      constexpr Revolutions r0 = 1_gon;
        constexpr Gons r1 = 1_rev;
        static_assert(r1.count_unsafe() == 400.0);
        constexpr Degrees r2 = 1_rev;
        static_assert(r2.count_unsafe() == 360.0);
        static_assert(!std::is_convertible_v<Gons, Revolutions>);
        static_assert( std::is_constructible_v<Revolutions, Gons>);

        constexpr auto t0 = 1_rad;
        constexpr auto t1 = 1_deg;
        static_assert( std::is_constructible_v<Radians, Degrees>);
        static_assert(!std::is_convertible_v<Degrees, Radians>);
        constexpr auto t2 = Radians(1_deg);
        constexpr auto t3 = 1_gon;
        constexpr auto t4 = 1_rev;
        static_assert( std::is_constructible_v<Gons, Degrees>);
        static_assert(!std::is_convertible_v<Degrees, Gons>);
        static_assert( std::is_constructible_v<Revolutions, Degrees>);
        static_assert(!std::is_convertible_v<Degrees, Revolutions>);

        static_assert( Compiles<Add, Radians, Radians>);
        static_assert(!Compiles<Add, Radians, Degrees>);
        static_assert(!Compiles<Add, Radians, Gons>);
        static_assert(!Compiles<Add, Radians, Revolutions>);

        static_assert( Compiles<Add, Degrees, Degrees>);
        static_assert(!Compiles<Add, Degrees, Radians>);
        static_assert(!Compiles<Add, Degrees, Gons>);
        static_assert( Compiles<Add, Degrees, Revolutions>);
        constexpr auto s0 = 1_deg; // + Degrees(1_gon);
        constexpr auto s1 = 1_deg + 1_rev;
        constexpr auto s2 = s0 + s1;
        constexpr auto s3 = s0 - s1;

        static_assert( Compiles<Add, Gons, Gons>);
        static_assert(!Compiles<Add, Gons, Radians>);
        static_assert(!Compiles<Add, Gons, Degrees>);
        static_assert( Compiles<Add, Gons, Revolutions>);

        //constexpr auto s4 = 1_deg + 1_gon;
        //constexpr auto s5 = 1_gon + 1_deg;
        //static_assert(std::is_same_v<decltype(s4), decltype(s5)>);

        constexpr auto p0 = 1_deg * 1_gon;
        constexpr auto p1 = 1_gon * 1_deg;
        static_assert(IsSame<decltype(p0), decltype(p1)>);

        constexpr auto q0 = 1_deg / 1_gon;
        constexpr auto q1 = 1_gon / 1_deg;
        static_assert(!IsSame<decltype(q0), decltype(q1)>);
        //constexpr auto sq0 = q0 + q1;
        //constexpr auto sq1 = q0 - q1;
        //Incomplet<decltype(q0)>{};
        //Incomplet<decltype(q1)>{};

        static_assert( Compiles<Add, Revolutions, Revolutions>);
        static_assert(!Compiles<Add, Revolutions, Radians>);
        static_assert( Compiles<Add, Revolutions, Degrees>);
        static_assert( Compiles<Add, Revolutions, Gons>);
    }

    {
        constexpr Metres t1(1_m);
        constexpr auto t2 = 1_cm;
        constexpr auto t3 = 1_mm;
        constexpr auto ts = t1 + t2 + t3;
        static_assert( std::is_convertible_v<Metres, Centimetres>);
        static_assert( std::is_convertible_v<Metres, Millimetres>);
        static_assert( std::is_convertible_v<Centimetres, Millimetres>);
        static_assert(!std::is_convertible_v<Centimetres, Inches>);
        static_assert(!std::is_convertible_v<Centimetres, Feet>);
        static_assert(!std::is_convertible_v<Metres, Miles>);
        constexpr auto tx = Centimetres(1_in);
        static_assert(tx.count_unsafe() == 2.54);

        static_assert(compare(1_m, 1_m) == 0);
        static_assert(compare(1_in, 1_cm) > 0);
        static_assert(compare(1_mm, 1_ft) < 0);

        constexpr auto s0 = 1_yd;
        static_assert(s0.count_unsafe() == 1.0);
        constexpr Feet x0 = s0;
        constexpr auto s1 = 1_yd + 1_ft;
        static_assert(s1.count_unsafe() == 4.0);
        constexpr Inches x1 = 1_ft;
        static_assert(x1.count_unsafe() == 12.0);
        constexpr auto s2 = 1_yd + 1_ft + 1_in;
        static_assert(s2.count_unsafe() == 49.0);
        static_assert(!std::is_convertible_v<Inches, Yards>);
        static_assert( std::is_constructible_v<Yards, Inches>);
        constexpr Yards x2 = Yards(1_in);
        constexpr Inches x3 = 1_yd;
    }

    {
        //struct ReinforcementContentKind
        //    : Kind<ReinforcementContentKind, Dimension<dim::Length> /* cm^2/m */> {};

        //using SquareCentimetresPerMetre = QuantityT<decltype(1_cm2 / 1_m), class AreaPerLength>;
        //using SquareMetresPerMetre      = QuantityT<decltype(1_m2 / 1_m), class AreaPerLength>;

        using SimplifedSquareMetresPerMetre = decltype((1_m*1_m)/1_m);

        static_assert(IsSame<Metres, SimplifedSquareMetresPerMetre>);
        static_assert( std::is_convertible_v<SimplifedSquareMetresPerMetre, Metres>);
        static_assert( std::is_constructible_v<Metres, SimplifedSquareMetresPerMetre>);

        static_assert(!std::is_convertible_v<SimplifedSquareMetresPerMetre, SquareMetresPerMetre>);
        static_assert( std::is_constructible_v<SquareMetresPerMetre, SimplifedSquareMetresPerMetre>);

        constexpr auto t0 = 1_m * 1_m;
        constexpr auto t1 = t0 / 1_m;
        static_assert( std::is_convertible_v<decltype(t1), Metres>);
        static_assert( std::is_constructible_v<Metres, decltype(t1)>);
        constexpr auto t2 = (1_m * 1_m) / 1_m;
        static_assert(IsSame<decltype(t2), Metres>);
        //constexpr decltype(t2) t22 = t1;
        constexpr decltype(t1) t22 = t2;
        //using ReinfCont = decltype(t1);
        //using ReinfCont = decltype(quantity_kind_cast<kinds::ReinforcementContent>(1_cm2 / 1_m));
        //constexpr ReinfCont t3(t2);
        //constexpr ReinfCont t3(quantity_kind_cast<ReinfCont>(t2));
        constexpr SquareCentimetresPerMetre t4(1_cm2 / 1_m);
        static_assert(!std::is_convertible_v<decltype(1_cm2/1_m), SquareCentimetresPerMetre>);
        static_assert( std::is_constructible_v<SquareCentimetresPerMetre, decltype(1_cm2/1_m)>);
        constexpr decltype(1_cm2 / 1_m) t44(t4);
        static_assert(t4.count_unsafe() == 1.0);
        static_assert(!std::is_assignable_v<SquareCentimetresPerMetre, decltype((1_cm / 1_m) * 1_cm)>);
        //static_assert( std::is_assignable_v<ReinfCont, decltype((1_cm / 1_m) * 1_cm)>);
        static_assert(!std::is_assignable_v<decltype(t1), Metres>);

        constexpr SquareCentimetresPerMetre t5(1_m2 / 1_cm);
        static_assert(t5.count_unsafe() == 1'000'000.0);
    }
    {
        static_assert( Compiles<AssignAdd, Millimetres, Centimetres>);
        static_assert(!Compiles<AssignAdd, Centimetres, Millimetres>);
        static_assert((1_cm + 1_mm).count_unsafe() == 11.0);
    }
    {
        constexpr auto energy = 1_J;          // J
        constexpr auto torque = 1_J/1_rad; // J/rad
        constexpr auto radian = energy / torque;
        //static_assert(!std::is_convertible_v<decltype(radian), Radians>);
        //static_assert( std::is_constructible_v<Radians, decltype(radian)>);
    }
    {
        using TonsPerCubicMetre = decltype(1_t/1_m3);
        using GramsPerCubicDecimetre = decltype(1_g/1_dm3);
        constexpr TonsPerCubicMetre w(1.0);
        constexpr Grams     m0 = w * 1_m3;
        constexpr Kilograms m1 = w * 1_m3;
        constexpr Tonnes    m2 = w * 1_m3;
        constexpr Tonnes    m3 = Tonnes(w * 1_dm3);
        constexpr Kilograms m4 = w * 1_dm3;
        constexpr Kilograms m5 = GramsPerCubicDecimetre(1.0) * 1_m3;
        constexpr Tonnes    m6 = Tonnes(GramsPerCubicDecimetre(1.0) * 1_m3);
    }
}

static void test2()
{
    {
        Millimetres mm;
        mm += 1_mm;
        mm += 1_cm;
        mm += 1_m;
        assert(mm == 1000_mm + 10_mm + 1_mm);

        Centimetres cm;
        static_assert(!Compiles<AssignAdd, Centimetres, Millimetres>);
        cm += 1_cm;
        cm += 1_m;
        assert(cm == 100_cm + 1_cm);

        Metres m;
        static_assert(!Compiles<AssignAdd, Metres, Millimetres>);
        static_assert(!Compiles<AssignAdd, Metres, Centimetres>);
        m += 1_m;
        assert(m == 1_m);
    }
    {
        Millimetres mm{1.0};
        mm += 1_cm;
        assert(mm.count_unsafe() == 11.0);
    }
    {
        using Width  = Tagged<Millimetres, class _width>;
        using Height = Tagged<Millimetres, class _height>;

        static_assert(!Compiles<Add, Width, Millimetres>);

        constexpr Width w1(1_cm);
        constexpr Height h1(1_cm + 1_mm);
        Width w;
        w *= 2.0;
        //w += 1_m;
        //w += h1;

        Height h;
        h += Height(1_cm);
        //h += 1_cm;
      //h += w;
        //h += w.value();
      //h = w;

      //const auto xxx = w + h;
        //const auto xxx = w + 1_m;

        const Centimetres mmm = 1_m + 1_km + 1_cm;
    }
    {
        using Length   = Tagged<Millimetres, class _length>;
        using Height   = Tagged<Millimetres, class _height>;
        using Width    = Tagged<Millimetres, class _width>;
        using Position = QuantityPoint<Length>;

        const Height h(0);

        Length x;
        Position y;
        y += x;
        x += Length(h);
        //x += h;
        //y += y;
        x += x;
        //x += y;
        const auto t = Millimetres(x);
        const auto z = x + y;
        const auto w = y + x;
        x = y - y;

        //Width W;
        //Height H;
        //W = H;
        //W += H;

        //Width W;
        // W += Height(1_cm);
        // W = x;
        //W = 1_cm;
        //W = x.value();

        constexpr Height h2(20.0);
        constexpr Width w2(30.0);
        constexpr auto p = h2 * w2;
        constexpr auto p2 = SquareCentimetres(p);
        assert(p2 == 6_cm2);
    }
    {
        using Duration = Seconds; // QuantityT<Seconds, class _duration>;
        using TimePoint = QuantityPoint<Duration>;
    }
#if 1
    {
        using ReinfCont = Tagged<decltype(1_cm2/1_m), class AreaPerLength>;

        ReinfCont r(3.1415_cm2 / 1_m);
        // r = 3.1415_cm2 / 1_m;
        r += r;
        r += r / 2.0;

        constexpr auto x = ReinfCont(3.1415_cm2 / 1_m);
//      constexpr auto y = x * 1_m;
        constexpr auto z = x * 1_m;
    }
#endif
    {
        //using Rainfall = Tagged<decltype(1_L/1_m2), class LitresPerSquareMeter>;

        //const Rainfall r(2.71828_L / 1_m2);
    }
    {
        //constexpr Steradians sr = 1_rad * 1_rad;
        //constexpr Radians rad = sr / 1_rad;
    }
    {
        //constexpr SquareCentimetres x = 1_cm * 1_cm;
        //constexpr Centimetres y = x / 1_cm;
    }
}

static void test3()
{
    {
        constexpr auto t01 = 1_cm * 1_cm;
        constexpr auto t02 = 1_cm2;
        const auto xxx = fma(1_cm, 1_cm, 1_cm2);
    }
}