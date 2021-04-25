#include "doctest.h"

#include "Unit.h"
#include "UnitLiterals.h"
#include "UnitMath.h"

using namespace uom;
using namespace uom::literals;

inline constexpr Kilometers EarthRadius = 6371_km;

using Latitude
    = TaggedQuantity<Degrees, class _latitude>;
using Longitude
    = TaggedQuantity<Degrees, class _longitude>;
using SphericalDistance
    = decltype(1_km * 1_rad);

struct Position {
    Latitude lat;
    Longitude lon;
};

static SphericalDistance spherical_distance_cosine(Position pos1, Position pos2)
{
    const Radians lat1 = convert_to<Radians>(pos1.lat);
    const Radians lat2 = convert_to<Radians>(pos2.lat);
    const Radians lon1 = convert_to<Radians>(pos1.lon);
    const Radians lon2 = convert_to<Radians>(pos2.lon);

    // Spherical law of cosines
    const auto central_angle = acos(sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon2 - lon1));

    return EarthRadius * central_angle;
}

static SphericalDistance spherical_distance_haversine(Position pos1, Position pos2)
{
    const Radians lat1 = convert_to<Radians>(pos1.lat);
    const Radians lat2 = convert_to<Radians>(pos2.lat);
    const Radians lon1 = convert_to<Radians>(pos1.lon);
    const Radians lon2 = convert_to<Radians>(pos2.lon);

    // Haversine formula
    const auto sin_lat = sin((lat2 - lat1) / 2);
    const auto sin_lon = sin((lon2 - lon1) / 2);
    const auto central_angle = 2 * asin(sqrt(sin_lat * sin_lat + cos(lat1) * cos(lat2) * sin_lon * sin_lon));

    return EarthRadius * central_angle;
}

TEST_CASE("Spherical distance - 1")
{
    const Latitude  lat1(Degrees( 50.94130));
    const Longitude lon1(Degrees(  6.95828));
    const Latitude  lat2(Degrees( 52.51664));
    const Longitude lon2(Degrees( 13.37760));

    const Position pos1{lat1, lon1}; // Koelner Dom
    const Position pos2{lat2, lon2}; // Brandenbuger Tor

    const SphericalDistance dist_1 = round(spherical_distance_cosine(pos1, pos2));
    const SphericalDistance dist_2 = round(spherical_distance_haversine(pos1, pos2));

    const double km_1 = count_as<Kilometers>(dist_1);
    const double km_2 = count_as<Kilometers>(dist_2);
    CHECK(km_1 == 475);
    CHECK(km_2 == 475);
}

TEST_CASE("Spherical distance - 2")
{
    const Latitude  lat1(Degrees( 59.33797));
    const Longitude lon1(Degrees( 18.07269));
    const Latitude  lat2(Degrees(-32.84555));
    const Longitude lon2(Degrees( 18.36146));

    const Position pos1{lat1, lon1}; // Stockholm
    const Position pos2{lat2, lon2}; // Kapstadt

    const SphericalDistance dist_1 = round(spherical_distance_cosine(pos1, pos2));
    const SphericalDistance dist_2 = round(spherical_distance_haversine(pos1, pos2));

    const double km_1 = count_as<Kilometers>(dist_1);
    const double km_2 = count_as<Kilometers>(dist_2);
    CHECK(km_1 == 10250);
    CHECK(km_2 == 10250);
}

struct DMS
{
    Degrees    d;
    ArcMinutes m;
    ArcSeconds s;
};

static DMS DmsFromDegrees(Degrees deg)
{
    const auto d = trunc(deg);
    const auto m = trunc(ArcMinutes(deg - d));
    return {d, m, deg - d - m};
}

TEST_CASE("DMS")
{
    {
        const Degrees degs(50.94130);

        const DMS dms = DmsFromDegrees(degs);
        CHECK(dms.d == Degrees(50));
        CHECK(dms.m == ArcMinutes(56));
        CHECK(dms.s == ArcSeconds(28.679999999993697));

        const ArcSeconds secs = dms.d + dms.m + dms.s;
        CHECK(secs == degs);

        const Degrees degs2(secs);
        CHECK(degs2 == degs);
    }
    {
        const Degrees degs( 6.95828);

        const DMS dms = DmsFromDegrees(degs);
        CHECK(dms.d == Degrees(6));
        CHECK(dms.m == ArcMinutes(57));
        CHECK(dms.s == ArcSeconds(29.808000000000874));

        const ArcSeconds secs = dms.d + dms.m + dms.s;
        CHECK(secs == degs);

        const Degrees degs2(secs);
        CHECK(degs2 == degs);
    }
}
