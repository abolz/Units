#include "Unit.h"
#include "UnitLiterals.h"
#include "UnitMath.h"

using namespace uom;
using namespace uom::literals;

// See:
// https://github.com/mpusz/units/tree/master/example

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

static SphericalDistance spherical_distance(Position pos1, Position pos2)
{
    const Radians lat1 = cast<Radians>(pos1.lat);
    const Radians lat2 = cast<Radians>(pos2.lat);
    const Radians lon1 = cast<Radians>(pos1.lon);
    const Radians lon2 = cast<Radians>(pos2.lon);

    if (1)
    {
        // Spherical law of cosines
        const auto central_angle = acos(sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon2 - lon1));
        return EarthRadius * central_angle;
    }
    else
    {
        // The haversine formula
        const auto sin_lat = sin(lat2 - lat1) / 2;
        const auto sin_lon = sin(lon2 - lon1) / 2;
        const auto central_angle = 2 * asin(sqrt(sin_lat * sin_lat + cos(lat1) * cos(lat2) * sin_lon * sin_lon));
        return EarthRadius * central_angle;
    }
}
