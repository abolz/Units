#include "Unit.h"
#include "UnitLiterals.h"
#include "UnitMath.h"

using namespace uom;
using namespace uom::literals;

struct LocalHeight1Zero {
    static constexpr double value = 200.0;
    [[nodiscard]] constexpr double operator()() const noexcept {
        return value;
    }
};

#if 0
struct LocalHeight2Zero {
    static constexpr double value = 100.0;
    [[nodiscard]] constexpr double operator()() const noexcept {
        return value;
    }
};
#else
using LocalHeight2Zero = std::integral_constant<int32_t, 100>;
#endif

static void test1()
{
    using GlobalHeight = Absolute<Meters>;
    using LocalHeight1 = Absolute<Meters, LocalHeight1Zero>;
    using LocalHeight2 = Absolute<Meters, LocalHeight2Zero>;

    constexpr LocalHeight1 h1(123.0);
    static_assert(h1.count_internal() == 123.0);
    static_assert(h1.count<GlobalHeight>() == 323.0);

    constexpr LocalHeight2 h2 = convert_to<LocalHeight2>(h1);
    static_assert(h2.count_internal() == 223.0);
    static_assert(h2.count<GlobalHeight>() == 323.0);

    constexpr GlobalHeight hg = convert_to<GlobalHeight>(h1);
    static_assert(hg.count_internal() == 323.0);
    static_assert(hg.count<LocalHeight1>() == 123.0);
    static_assert(hg.count<LocalHeight2>() == 223.0);
}