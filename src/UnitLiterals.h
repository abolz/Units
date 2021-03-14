// Copyright 2021 Alexander Bolz
//
// Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at https://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include "Unit.h"

namespace uom::literals {

//--------------------------------------------------------------------------------------------------
// One

[[nodiscard]] constexpr auto operator""_q(long double x) noexcept
{
    return Dimensionless(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_q(unsigned long long x) noexcept
{
    return Dimensionless(static_cast<double>(x));
}

//--------------------------------------------------------------------------------------------------
// Length

[[nodiscard]] constexpr auto operator""_mm(long double x) noexcept
{
    return Millimetres(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_mm(unsigned long long x) noexcept
{
    return Millimetres(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_cm(long double x) noexcept
{
    return Centimetres(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_cm(unsigned long long x) noexcept
{
    return Centimetres(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_dm(long double x) noexcept
{
    return Decimetres(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_dm(unsigned long long x) noexcept
{
    return Decimetres(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_m(long double x) noexcept
{
    return Metres(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_m(unsigned long long x) noexcept
{
    return Metres(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_km(long double x) noexcept
{
    return Kilometres(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_km(unsigned long long x) noexcept
{
    return Kilometres(static_cast<double>(x));
}

[[nodiscard]] constexpr auto operator""_in(long double x) noexcept
{
    return Inches(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_in(unsigned long long x) noexcept
{
    return Inches(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_ft(long double x) noexcept
{
    return Feet(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_ft(unsigned long long x) noexcept
{
    return Feet(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_yd(long double x) noexcept
{
    return Yards(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_yd(unsigned long long x) noexcept
{
    return Yards(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_mi(long double x) noexcept
{
    return Miles(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_mi(unsigned long long x) noexcept
{
    return Miles(static_cast<double>(x));
}

//--------------------------------------------------------------------------------------------------
// Area

[[nodiscard]] constexpr auto operator""_mm2(long double x) noexcept
{
    return SquareMillimetres(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_mm2(unsigned long long x) noexcept
{
    return SquareMillimetres(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_cm2(long double x) noexcept
{
    return SquareCentimetres(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_cm2(unsigned long long x) noexcept
{
    return SquareCentimetres(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_m2(long double x) noexcept
{
    return SquareMetres(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_m2(unsigned long long x) noexcept
{
    return SquareMetres(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_km2(long double x) noexcept
{
    return SquareKilometres(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_km2(unsigned long long x) noexcept
{
    return SquareKilometres(static_cast<double>(x));
}

//--------------------------------------------------------------------------------------------------
// Volume

[[nodiscard]] constexpr auto operator""_mm3(long double x) noexcept
{
    return CubicMillimetres(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_mm3(unsigned long long x) noexcept
{
    return CubicMillimetres(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_cm3(long double x) noexcept
{
    return CubicCentimetres(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_cm3(unsigned long long x) noexcept
{
    return CubicCentimetres(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_dm3(long double x) noexcept
{
    return CubicDecimetres(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_dm3(unsigned long long x) noexcept
{
    return CubicDecimetres(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_m3(long double x) noexcept
{
    return CubicMetres(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_m3(unsigned long long x) noexcept
{
    return CubicMetres(static_cast<double>(x));
}

//[[nodiscard]] constexpr auto operator""_L(long double x) noexcept
//{
//    return Litres(static_cast<double>(x));
//}
//[[nodiscard]] constexpr auto operator""_L(unsigned long long x) noexcept
//{
//    return Litres(static_cast<double>(x));
//}

//--------------------------------------------------------------------------------------------------
// Time

[[nodiscard]] constexpr auto operator""_ms(long double x) noexcept
{
    return Milliseconds(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_ms(unsigned long long x) noexcept
{
    return Milliseconds(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_s(long double x) noexcept
{
    return Seconds(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_s(unsigned long long x) noexcept
{
    return Seconds(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_min(long double x) noexcept
{
    return Minutes(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_min(unsigned long long x) noexcept
{
    return Minutes(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_h(long double x) noexcept
{
    return Hours(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_h(unsigned long long x) noexcept
{
    return Hours(static_cast<double>(x));
}

//--------------------------------------------------------------------------------------------------
// Frequency

[[nodiscard]] constexpr auto operator""_Hz(long double x) noexcept
{
    return Hertz(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_Hz(unsigned long long x) noexcept
{
    return Hertz(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_kHz(long double x) noexcept
{
    return Kilohertz(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_kHz(unsigned long long x) noexcept
{
    return Kilohertz(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_MHz(long double x) noexcept
{
    return Megahertz(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_MHz(unsigned long long x) noexcept
{
    return Megahertz(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_GHz(long double x) noexcept
{
    return Gigahertz(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_GHz(unsigned long long x) noexcept
{
    return Gigahertz(static_cast<double>(x));
}

//--------------------------------------------------------------------------------------------------
// Mass

[[nodiscard]] constexpr auto operator""_g(long double x) noexcept
{
    return Grams(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_g(unsigned long long x) noexcept
{
    return Grams(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_kg(long double x) noexcept
{
    return Kilograms(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_kg(unsigned long long x) noexcept
{
    return Kilograms(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_t(long double x) noexcept
{
    return Tonnes(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_t(unsigned long long x) noexcept
{
    return Tonnes(static_cast<double>(x));
}

//--------------------------------------------------------------------------------------------------
// Velocity

[[nodiscard]] constexpr auto operator""_m_per_s(long double x) noexcept
{
    return MetresPerSecond(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_m_per_s(unsigned long long x) noexcept
{
    return MetresPerSecond(static_cast<double>(x));
}

[[nodiscard]] constexpr auto operator""_km_per_h(long double x) noexcept
{
    return KilometresPerHour(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_km_per_h(unsigned long long x) noexcept
{
    return KilometresPerHour(static_cast<double>(x));
}

//--------------------------------------------------------------------------------------------------
// Temperature

[[nodiscard]] constexpr auto operator""_K(long double x) noexcept
{
    return Kelvin(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_K(unsigned long long x) noexcept
{
    return Kelvin(static_cast<double>(x));
}

//--------------------------------------------------------------------------------------------------
// Plane angle

[[nodiscard]] constexpr auto operator""_rad(long double x) noexcept
{
    return Radians(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_rad(unsigned long long x) noexcept
{
    return Radians(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_deg(long double x) noexcept
{
    return Degrees(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_deg(unsigned long long x) noexcept
{
    return Degrees(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_gon(long double x) noexcept
{
    return Gons(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_gon(unsigned long long x) noexcept
{
    return Gons(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_rev(long double x) noexcept
{
    return Revolutions(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_rev(unsigned long long x) noexcept
{
    return Revolutions(static_cast<double>(x));
}

//--------------------------------------------------------------------------------------------------
// Solid angle

[[nodiscard]] constexpr auto operator""_sr(long double x) noexcept
{
    return Steradians(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_sr(unsigned long long x) noexcept
{
    return Steradians(static_cast<double>(x));
}

//--------------------------------------------------------------------------------------------------
// Force

[[nodiscard]] constexpr auto operator""_N(long double x) noexcept
{
    return Newtons(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_N(unsigned long long x) noexcept
{
    return Newtons(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_kN(long double x) noexcept
{
    return Kilonewtons(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_kN(unsigned long long x) noexcept
{
    return Kilonewtons(static_cast<double>(x));
}

//--------------------------------------------------------------------------------------------------
// Energy

[[nodiscard]] constexpr auto operator""_J(long double x) noexcept
{
    return Joules(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_J(unsigned long long x) noexcept
{
    return Joules(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_kJ(long double x) noexcept
{
    return Kilojoules(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_kJ(unsigned long long x) noexcept
{
    return Kilojoules(static_cast<double>(x));
}

//[[nodiscard]] constexpr auto operator""_Nm(long double x) noexcept
//{
//    return Joules(static_cast<double>(x));
//}
//[[nodiscard]] constexpr auto operator""_Nm(unsigned long long x) noexcept
//{
//    return Joules(static_cast<double>(x));
//}

//--------------------------------------------------------------------------------------------------
// Torque

//[[nodiscard]] constexpr auto operator""_J_per_rad(long double x) noexcept
//{
//    return NewtonMetres(static_cast<double>(x));
//}
//[[nodiscard]] constexpr auto operator""_J_per_rad(unsigned long long x) noexcept
//{
//    return NewtonMetres(static_cast<double>(x));
//}

[[nodiscard]] constexpr auto operator""_Nm_per_rad(long double x) noexcept
{
    return NewtonMetres(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_Nm_per_rad(unsigned long long x) noexcept
{
    return NewtonMetres(static_cast<double>(x));
}

//--------------------------------------------------------------------------------------------------
// Data

[[nodiscard]] constexpr auto operator""_b(long double x) noexcept
{
    return Bits(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_b(unsigned long long x) noexcept
{
    return Bits(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_B(long double x) noexcept
{
    return Bytes(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_B(unsigned long long x) noexcept
{
    return Bytes(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_kB(long double x) noexcept
{
    return Kilobytes(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_kB(unsigned long long x) noexcept
{
    return Kilobytes(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_MB(long double x) noexcept
{
    return Megabytes(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_MB(unsigned long long x) noexcept
{
    return Megabytes(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_GB(long double x) noexcept
{
    return Gigabytes(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_GB(unsigned long long x) noexcept
{
    return Gigabytes(static_cast<double>(x));
}

//--------------------------------------------------------------------------------------------------
// Area per Length

[[nodiscard]] constexpr auto operator""_cm2_per_m(long double x) noexcept
{
    return SquareCentimetresPerMetre(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_cm2_per_m(unsigned long long x) noexcept
{
    return SquareCentimetresPerMetre(static_cast<double>(x));
}

[[nodiscard]] constexpr auto operator""_m2_per_m(long double x) noexcept
{
    return SquareMetresPerMetre(static_cast<double>(x));
}
[[nodiscard]] constexpr auto operator""_m2_per_m(unsigned long long x) noexcept
{
    return SquareMetresPerMetre(static_cast<double>(x));
}

} // namespace uom::literals
