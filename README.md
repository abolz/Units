Units of measure
====

A small, _simple to use_, (C++17) compile-time SI-unit library.

Quantities
----

The primary class in this library is the `Quantity<U>` template, where `U` denotes a `Unit`.
The underlying value type currently always is `double` (but this will hopefully be changed
in the near future). The library provides a few typedef's and user-defined literals ready
for use. E.g.
```c++
Metres        length    = 1_m;
Kilograms     mass      = 1_kg;
Seconds       time      = 1_s;
SquareMetres  area      = length * length;
Hertz         frequency = 1 / time;
Metres        mean      = sqrt(area);
Dimensionless ratio     = 1_m / 2_m;
```

Strong units
----

In some contexts, different quanities with the same physical unit are not interchangeable.
To support such quanties, all predefined `Quantity`s can be (re-)tagged. E.g.:
```c++
using Width  = Tagged<Millimetres, struct _width>;
using Height = Tagged<Millimetres, struct _height>;

Width  w(1.0);      // 1 mm
Height h(2.0_cm);   // 2 cm

#if 0
w = h; // will not compile
w + h; // will not compile
#endif

const auto area = w * h; // (some kind of "area")
w = area / h; // works
h = area / w; // works
const auto value = area.count<SquareCentimetres>();

const auto a = SquareMillimetres(area); // explicit cast works
```

The library differentiates "simple" and "complex" (or "tagged") units.

Quantity conversions/arithmetic
----

Implicit conversions are supported when the source quantity is a multiple of the target quantity. E.g. `Centimetres` are implicitly convertible to `Millimetres` - but not vice versa:
```c++
void f(Millimetres mm);

f(1_cm);
```
works fine.

Angles
----

Absolute quantities
----

are not correctly supported yet.