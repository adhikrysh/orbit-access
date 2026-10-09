# orbit-access

a leo satellite gives a ground station about six to eight minutes a pass, a few passes a day. miss the window and you wait hours. the earth doesn't help: it's fat at the equator, so J2 drags inclined orbits around, and the station is spinning under the satellite at ~465 m/s. get either wrong and the dish is pointing at a cloud.

this propagates an orbit with J2, puts a WGS84 station on a spinning earth, and finds every window above the elevation mask.

```text
state at t0 -> adaptive J2 propagation -> WGS84 station + earth spin
            -> elevation above mask -> bracket -> bisect -> pass
```

default run (circular, 7,000 km, 51.6°, station at 1.35° N 103.82° E) gives four passes above 10° in a day: 390.278, 449.958, 455.508 and 404.475 s. model outputs, not a forecast for anything real.

![One-day comparison](docs/orbit-comparison.png)

## how

dynamics are central gravity plus J2 in [dynamics.cpp](src/dynamics.cpp), integrated with Dormand-Prince 5(4) using separate tolerances for km and km/s (1e-8 km, 1e-11 km/s). it either reaches the requested time or throws, it never quietly hands back half a trajectory.

with J2 on, energy and z angular momentum are conserved but the full angular momentum vector isn't, which matters because it decides which conservation test is even valid ([numerics.md](docs/numerics.md)).

the station uses the WGS84 ellipsoid normal as up, not the line from earth's centre, rotated with an explicit Greenwich angle. passes are found by sampling every 10 s, then bisecting each crossing to 1e-5 s, re-propagating from the same endpoint every time instead of interpolating. a pass shorter than the sample step can still sneak through between samples, and that's stated, not hidden.

## tests

150 elliptic two-body runs against an independent Kepler solver (under 0.5 m and 0.5 mm/s), conservation over long runs, J2 acceleration against a numerical potential gradient, node drift within 1% of secular theory, analytic rise and set times, and the same pass list when the sample step drops from 31 s to 7 s. [ci](.github/workflows/check.yml) runs it under asan and ubsan.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build
ctest --test-dir build --output-on-failure
./build/orbit-access passes 86400 1.3521 103.8198 10 j2
```

no drag, no moon, no sun, no UTC or polar motion, no refraction. good for poking at the numerics, not for booking a real antenna.
