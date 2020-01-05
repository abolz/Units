// Copyright Alexander Bolz 2019
// Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at https://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <cassert>
#include <type_traits>

namespace sc {

//==================================================================================================
// Torsor
//==================================================================================================

template <typename T>
class Torsor;

template <typename T>
struct IsTorsor : std::false_type {};

template <typename T>
struct IsTorsor<Torsor<T>> : std::true_type {};

//--------------------------------------------------------------------------------------------------
//
//--------------------------------------------------------------------------------------------------

template <typename T>
class Torsor final
{
    template <typename T2> friend class Torsor;

    static_assert(!IsTorsor<T>::value, "T must not be a Torsor");
    // TODO:
    // Check more conditions on type T here...

    T point_;

public:
    constexpr Torsor() noexcept = default;
    constexpr Torsor(const Torsor&) noexcept = default;
    constexpr Torsor& operator=(const Torsor&) noexcept = default;

    template <typename U,
        std::enable_if_t< !IsTorsor<U>::value && std::is_constructible< T, U&& >::value, int > = 0>
    constexpr explicit Torsor(U&& init) noexcept
        : point_(static_cast<U&&>(init))
    {
    }

    template <typename T2,
        std::enable_if_t< std::is_constructible< T, T2 >::value, int > = 0>
    constexpr explicit Torsor(const Torsor<T2>& init) noexcept
        : point_(static_cast<T>(init.point_))
    {
    }

    //constexpr T& point() noexcept {
    //    return point_;
    //}

    //constexpr const T& point() const noexcept {
    //    return point_;
    //}

    //----------------------------------------------------------------------------------------------
    // Arithmetic

    // Torsor + Element -> Torsor
    template <typename U,
        std::enable_if_t< !IsTorsor<U>::value, int > = 0>
    constexpr friend auto operator+(const Torsor& lhs, const U& rhs) noexcept
        -> Torsor<decltype(lhs.point_ + rhs)>
    {
        return Torsor<decltype(lhs.point_ + rhs)>(lhs.point_ + rhs);
    }

    // Element + Torsor -> Torsor
    template <typename U,
        std::enable_if_t< !IsTorsor<U>::value, int > = 0>
    constexpr friend auto operator+(const U& lhs, const Torsor& rhs) noexcept
        -> Torsor<decltype(lhs + rhs.point_)>
    {
        return Torsor<decltype(lhs + rhs.point_)>(lhs + rhs.point_);
    }

    // Torsor - Torsor -> Element
    template <typename T2>
    constexpr friend auto operator-(const Torsor& lhs, const Torsor<T2>& rhs) noexcept
        -> decltype(( lhs.point_ - rhs.point_ ))
    {
        return lhs.point_ - rhs.point_;
    }

    // Torsor - Element -> Torsor
    template <typename U,
        std::enable_if_t< !IsTorsor<U>::value, int > = 0>
    constexpr friend auto operator-(const Torsor& lhs, const U& rhs) noexcept
        -> Torsor<decltype(lhs.point_ - rhs)>
    {
        return Torsor<decltype(lhs.point_ - rhs)>(lhs.point_ - rhs);
    }

    // Torsor + Element -> Torsor
    template <typename U,
        std::enable_if_t< !IsTorsor<U>::value, int > = 0>
    constexpr friend auto operator+=(Torsor& lhs, const U& rhs) noexcept
        -> decltype(( static_cast<void>(lhs.point_ += rhs), lhs ))
    {
        lhs.point_ += rhs;
        return lhs;
    }

    // Torsor - Element -> Torsor
    template <typename U,
        std::enable_if_t< !IsTorsor<U>::value, int > = 0>
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
        return lhs.point_ == rhs.point_;
    }

    template <typename T2,
        std::enable_if_t< std::is_convertible<decltype(std::declval<T>() != std::declval<T2>()), bool>::value, int > = 0>
    constexpr friend bool operator!=(const Torsor& lhs, const Torsor<T2>& rhs) noexcept
    {
        return lhs.point_ != rhs.point_;
    }

    template <typename T2,
        std::enable_if_t< std::is_convertible<decltype(std::declval<T>() < std::declval<T2>()), bool>::value, int > = 0>
    constexpr friend bool operator<(const Torsor& lhs, const Torsor<T2>& rhs) noexcept
    {
        return lhs.point_ < rhs.point_;
    }

    template <typename T2,
        std::enable_if_t< std::is_convertible<decltype(std::declval<T>() > std::declval<T2>()), bool>::value, int > = 0>
    constexpr friend bool operator>(const Torsor& lhs, const Torsor<T2>& rhs) noexcept
    {
        return lhs.point_ > rhs.point_;
    }

    template <typename T2,
        std::enable_if_t< std::is_convertible<decltype(std::declval<T>() <= std::declval<T2>()), bool>::value, int > = 0>
    constexpr friend bool operator<=(const Torsor& lhs, const Torsor<T2>& rhs) noexcept
    {
        return lhs.point_ <= rhs.point_;
    }

    template <typename T2,
        std::enable_if_t< std::is_convertible<decltype(std::declval<T>() >= std::declval<T2>()), bool>::value, int > = 0>
    constexpr friend bool operator>=(const Torsor& lhs, const Torsor<T2>& rhs) noexcept
    {
        return lhs.point_ >= rhs.point_;
    }
};

template <typename T> using Absolute = Torsor<T>;

} // namespace sc
