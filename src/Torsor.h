// Copyright Alexander Bolz 2019
// Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at https://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <cassert>
#include <type_traits>

namespace sc {

template <typename T>
class Torsor
{
    //static_assert(!IsTorsor<T>::value, "");

    T point_;

public:
    constexpr Torsor() noexcept = default;
    constexpr Torsor(const Torsor&) noexcept = default;
    constexpr Torsor& operator=(const Torsor&) noexcept = default;

    template <typename U,
        std::enable_if_t< !std::is_same<std::remove_cv_t<U>, Torsor>::value && std::is_constructible< T, U&& >::value, int > = 0>
    constexpr explicit Torsor(U&& init) noexcept
        : point_(static_cast<U&&>(init))
    {
    }

    // constexpr T& point() noexcept {
    //     return point_;
    // }

    constexpr const T& point() const noexcept {
        return point_;
    }

    //----------------------------------------------------------------------------------------------
    // Arithmetic

    // Torsor + Element -> Torsor
    template <typename U>
    constexpr friend auto operator+(const Torsor& lhs, const U& rhs) noexcept
        -> Torsor<decltype(lhs.point() + rhs)>
    {
        return Torsor<decltype(lhs.point() + rhs)>(lhs.point() + rhs);
    }

    // Torsor - Element -> Torsor
    template <typename T2>
    constexpr friend auto operator-(const Torsor& lhs, const Torsor<T2>& rhs) noexcept
        -> decltype(( lhs.point() - rhs.point() ))
    {
        return lhs.point() - rhs.point();
    }

    // Torsor - Torsor -> Element
    template <typename U>
    constexpr friend auto operator-(const Torsor& lhs, const U& rhs) noexcept
        -> Torsor<decltype(lhs.point() - rhs)>
    {
        return Torsor<decltype(lhs.point() - rhs)>(lhs.point() - rhs);
    }

    // Torsor += Element
    template <typename U>
    constexpr friend auto operator+=(Torsor& lhs, const U& rhs) noexcept
        -> decltype(( static_cast<void>(lhs.point_ += rhs), lhs ))
    {
        lhs.point_ += rhs;
        return lhs;
    }

    // Torsor -= Element
    template <typename U>
    constexpr friend auto operator-=(Torsor& lhs, const U& rhs) noexcept
        -> decltype(( static_cast<void>(lhs.point_ += rhs), lhs ))
    {
        lhs.point_ -= rhs;
        return lhs;
    }

    //----------------------------------------------------------------------------------------------
    // Comparisons

    template <typename T2,
        std::enable_if_t< std::is_convertible<decltype(std::declval<T>() == std::declval<T2>()), bool>::value, int > = 0>
    constexpr friend bool operator==(const Torsor& lhs, const Torsor<T2>& rhs) noexcept
    {
        return lhs.point() == rhs.point();
    }

    template <typename T2,
        std::enable_if_t< std::is_convertible<decltype(std::declval<T>() != std::declval<T2>()), bool>::value, int > = 0>
    constexpr friend bool operator!=(const Torsor& lhs, const Torsor<T2>& rhs) noexcept
    {
        return lhs.point() != rhs.point();
    }

    template <typename T2,
        std::enable_if_t< std::is_convertible<decltype(std::declval<T>() < std::declval<T2>()), bool>::value, int > = 0>
    constexpr friend bool operator<(const Torsor& lhs, const Torsor<T2>& rhs) noexcept
    {
        return lhs.point() < rhs.point();
    }

    template <typename T2,
        std::enable_if_t< std::is_convertible<decltype(std::declval<T>() > std::declval<T2>()), bool>::value, int > = 0>
    constexpr friend bool operator>(const Torsor& lhs, const Torsor<T2>& rhs) noexcept
    {
        return lhs.point() > rhs.point();
    }

    template <typename T2,
        std::enable_if_t< std::is_convertible<decltype(std::declval<T>() <= std::declval<T2>()), bool>::value, int > = 0>
    constexpr friend bool operator<=(const Torsor& lhs, const Torsor<T2>& rhs) noexcept
    {
        return lhs.point() <= rhs.point();
    }

    template <typename T2,
        std::enable_if_t< std::is_convertible<decltype(std::declval<T>() >= std::declval<T2>()), bool>::value, int > = 0>
    constexpr friend bool operator>=(const Torsor& lhs, const Torsor<T2>& rhs) noexcept
    {
        return lhs.point() >= rhs.point();
    }
};

template <typename T> using Absolute = Torsor<T>;

} // namespace sc
