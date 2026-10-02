# orbit-access

calculate how a satellite moves, then find when a ground station can see it. The library supports two-body gravity and Earth's J2 term, which models the largest effect of Earth's equatorial bulge.

the included comparison shows how a simple circular orbit changes when J2 is included. Both models start from the same state. After one day, their orbital planes and positions differ, even though each model conserves its own mechanical energy.

![One-day comparison](docs/orbit-comparison.png)

## build and run

requires C++20 and CMake 3.20+. It has no runtime dependencies.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 3
ctest --test-dir build --output-on-failure

./build/orbit-access propagate 86400 60 j2 > trajectory.csv
./build/orbit-access passes 86400 1.3521 103.8198 10 j2
```

the first command writes one day of states at 60-second intervals. The second finds passes above a 10° elevation mask for a station at the given latitude and longitude. Times are elapsed seconds. This is a synthetic scenario, not a prediction for a named spacecraft.

the default state has a 7,000 km geocentric radius and 51.6° inclination. To provide another state, append `x y z vx vy vz` in kilometres and kilometres per second. `two-body` disables J2. `--help` lists the arguments.

## model and search

the propagator integrates position and velocity with Dormand–Prince 5(4), an adaptive numerical method. Two estimates for each step give an error estimate. A failed step retries at a smaller size; a successful step can increase the next size. Position and velocity use separate absolute tolerances because their units and scales differ.

ground stations use WGS84 geodetic coordinates. The library rotates a station's position into the inertial frame and measures elevation from the station's local vertical. A transition from below to above the mask brackets a rise time. Bisection with fresh propagation refines that crossing. The same process finds the set time.

the default search samples every 10 seconds. A very short grazing pass can fall between samples and be missed. The root tolerance controls the precision of a detected crossing; it does not guarantee that every crossing was detected. A pass already underway at the start or still visible at the end carries a clipping flag.

## checks and limits

the tests compare 150 elliptic trajectories, including backward propagation, against a separate Kepler-equation solver. They also check long-run energy and angular momentum, the J2 force against a numerical potential gradient, node drift against first-order secular theory, and pass times against an analytic circular equatorial case.

time zero uses a user-defined Earth orientation. The CLI sets the Greenwich angle to zero, assumes constant Earth rotation, and sets station height to zero. It has no UTC conversion, precession, nutation, drag, lunar or solar gravity, or atmospheric refraction. The propagator rejects force evaluations inside its reference radius; it does not find impact times. Use it for controlled experiments, not operational contact predictions.

to recreate the plot, install `matplotlib` and `numpy` in a Python environment, then run `python examples/plot.py`. To run memory and undefined-behavior checks, configure another build with `-DORBIT_SANITIZERS=ON -DCMAKE_BUILD_TYPE=Debug`.

the public API is in [orbit.hpp](include/orbit/orbit.hpp). [Numerical notes](docs/numerics.md) define the model and its checks. Useful references are JPL's [two-body propagation documentation](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/prop2b_c.html) and the [Dormand–Prince method used by SciPy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.RK45.html).
