# orbit-access

Propagate a satellite's position and velocity, then find when a ground station can see it. The library supports two-body gravity and Earth's J2 term, which accounts for the largest effect of Earth's equatorial bulge.

The useful comparison is how much a simple circular orbit changes when J2 is included. Both models start from the same state. After a day, their orbital planes and positions differ even though each model conserves its own mechanical energy.

![One-day comparison](docs/orbit-comparison.png)

## Build and run

Requires C++20 and CMake 3.20+. No runtime dependencies.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 3
ctest --test-dir build --output-on-failure

./build/orbit-access propagate 86400 60 j2 > trajectory.csv
./build/orbit-access passes 86400 1.3521 103.8198 10 j2
```

The first command outputs one day of states at 60-second intervals. The second finds passes above a 10° elevation mask for a station at the given latitude and longitude. Times are elapsed seconds. This is a synthetic scenario, not a prediction for a named spacecraft.

The default state has a 7,000 km geocentric radius and 51.6° inclination. To supply another state, append `x y z vx vy vz` in kilometres and kilometres per second. `two-body` disables J2. `--help` lists the arguments.

## How it works

Position and velocity are integrated with the Dormand–Prince 5(4) method. Two estimates of each step provide an error estimate. A failed step is retried at a smaller size; a successful step can increase the next size. Position and velocity have separate absolute tolerances because they have different units and scales.

Ground stations use WGS84 geodetic coordinates. The library rotates their positions into the inertial frame and measures elevation against the station's local vertical. A change from below to above the mask brackets a rise time. Bisection with fresh propagation refines that crossing. The same process finds the set time.

The default search samples every 10 seconds. A very short grazing pass can fall between samples and be missed. The root tolerance controls the precision of a detected crossing; it does not guarantee that every crossing was detected. A pass already underway at the start or still visible at the end carries a clipping flag.

## Checks and limits

The tests compare 150 elliptic trajectories, including backward propagation, against a separate Kepler-equation solver. They also check long-run energy and angular momentum, the J2 force against a numerical potential gradient, node drift against first-order secular theory, and pass times against an analytic circular equatorial case.

Time zero has a user-defined Earth orientation. The CLI uses a Greenwich angle of zero, constant Earth rotation, and zero station height. There is no UTC conversion, precession, nutation, drag, lunar/solar gravity, or atmospheric refraction. The propagator rejects force evaluations inside its reference radius; it does not locate impact times. Use this model for controlled experiments, not operational contact predictions.

To repeat the plot, install `matplotlib` and `numpy` in a Python environment, then run `python examples/plot.py`. To run memory and undefined-behavior checks, configure another build with `-DORBIT_SANITIZERS=ON -DCMAKE_BUILD_TYPE=Debug`.

The public API is in [orbit.hpp](include/orbit/orbit.hpp). [Numerical notes](docs/numerics.md) define the model and its checks. Useful references are JPL's [two-body propagation documentation](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/prop2b_c.html) and the [Dormand–Prince method used by SciPy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.RK45.html).
