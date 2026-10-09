# orbit-access

a leo pass over a ground station lasts six to eight minutes in this model. four a day, about half an hour of contact in total. predict the window late and the pass is gone; the next one is hours away.

the earth makes this annoying. it's oblate, so J2 swings inclined orbit planes around the pole, and the station is moving under the satellite at about 465 m/s. ignore either one and the dish points at empty sky.

orbit-access propagates a satellite in an inertial frame and computes its elevation above a ground station on a rotating, oblate earth. it reports intervals above the station's elevation mask.

```text
state at t0 (km, km/s)
          |
adaptive J2 propagation
          |
satellite position at t
          |
WGS84 station + earth spin
          |
elevation above mask
          |
bracket -> bisection -> pass
```

the default CLI scenario starts with a circular orbit at a 7,000 km geocentric radius and 51.6 degrees inclination. it reports four above-10-degree windows for a one-day run at 1.3521 degrees north, 103.8198 degrees east: 390.278 s, 449.958 s, 455.508 s, and 404.475 s long. those are outputs of this stated model, not predictions for a named spacecraft.

![One-day comparison](docs/orbit-comparison.png)

## a model with stated conventions

the state is position in kilometres and velocity in kilometres per second. time is elapsed seconds. vectors live in a right-handed inertial frame whose positive z axis is earth's spin axis. library angles are radians; the CLI converts its latitude, longitude, and mask arguments from degrees. [orbit.hpp](include/orbit/orbit.hpp) defines these units and conventions.

the force model combines central gravity with J2. its potential per unit mass is:

```text
U = -mu/r [1 - J2/2 (R/r)^2
            (3z^2/r^2 - 1)]
```

units are km, km/s and km²/s² throughout. [dynamics.cpp](src/dynamics.cpp) evaluates the matching acceleration directly. setting J2 to zero leaves the two-body model. with J2 enabled, the orbit plane precesses; energy and z angular momentum remain conserved, while the full angular-momentum vector does not. [numerics.md](docs/numerics.md) makes that distinction explicit because it determines which conservation check is valid.

## propagation

[dynamics.cpp](src/dynamics.cpp) uses Dormand-Prince 5(4). position and velocity have separate absolute tolerances because kilometres and kilometres per second have different scales: 1e-8 km, 1e-11 km/s, with 1e-11 relative tolerance. accepted steps can grow to 60 s; rejected steps shrink. the propagator either reaches the requested time or throws when its minimum step or work budget is exhausted. it does not return a partial trajectory as a complete one.

## station geometry

[access.cpp](src/access.cpp) uses WGS84 geodetic coordinates for the station. the station's local vertical is the ellipsoid normal, which generally differs from its centre-to-station direction. the code rotates that position and normal into the inertial frame with an explicit Greenwich angle and a constant 7.2921150e-5 rad/s earth rotation rate. elevation is then the angle between the station-to-satellite line of sight and that local vertical.

## pass boundaries

the pass finder samples elevation every 10 s by default. when a sample pair changes visibility, it brackets the crossing and bisects it, propagating from the same left endpoint at each midpoint instead of interpolating elevation. the default root tolerance is 1e-5 s. pass records retain `clipped_at_start` and `clipped_at_end` when the search begins or ends inside a visible window. a crossing precision does not prove complete detection: a shorter grazing pass can fall entirely between samples, and a tangency without an above-mask interval is intentionally omitted.

## numerical checks

the test executable passes seven named checks locally. [test_orbit.cpp](tests/test_orbit.cpp) compares 150 elliptic, forward and backward two-body propagations with an independent Kepler-equation solver, requiring position error below 0.5 m and velocity error below 0.5 mm/s. it also checks long-run conservation, the J2 acceleration against a centred potential gradient, and node-drift sign and scale against first-order secular theory within a 1% envelope.

the access tests use an analytic circular equatorial case to check rise and set times, then reduce the sample step from 31 s to 7 s and require the same pass list. they also check WGS84 equatorial and polar geometry, overhead and antipodal elevation, clipped intervals, tangencies, invalid input, adaptive-step rejection, and exhausted work budgets. the [checked-in CI workflow](.github/workflows/check.yml) runs these tests on Ubuntu with address and undefined-behaviour sanitizers and on macOS in release mode.

## build and reproduce

requires C++20 and CMake 3.20+. it has no runtime dependencies.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 3
ctest --test-dir build --output-on-failure

./build/orbit-access propagate 86400 60 j2 > trajectory.csv
./build/orbit-access passes 86400 1.3521 103.8198 10 j2
```

the first command writes one day of states at 60-second intervals. the second searches passes above a 10-degree mask. append `x y z vx vy vz` to provide another initial state, in kilometres and kilometres per second; use `two-body` to disable J2. reproduce the figure with Matplotlib and NumPy using `python examples/plot.py`.

time zero uses a caller-defined earth orientation; the CLI uses zero Greenwich angle and station height. it has no UTC conversion, precession, nutation, polar motion, drag, lunar or solar gravity, atmospheric refraction, impact detection, or event-completeness guarantee. it is appropriate for inspecting numerical and geometric assumptions, not operational contact prediction.
