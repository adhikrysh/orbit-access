# Numerical model

All vectors are in a right-handed inertial frame with Earth's rotation axis along positive z. Distance is kilometres, time is seconds, and angles in the library are radians. A `State` contains position and velocity at an elapsed time; the propagator itself is autonomous and does not interpret calendar dates.

The potential per unit mass is

```text
U = -mu/r * [1 - J2/2 * (R/r)^2 * (3*z^2/r^2 - 1)]
```

Here `mu` is the gravitational parameter, `R` the reference radius, and `r` the distance from Earth's centre. Acceleration is the negative gradient of this potential. The tests evaluate that gradient with centred differences independently of the acceleration implementation.

With J2 set to zero, the force is central. Both specific energy and the full angular momentum vector should remain constant. With J2 enabled, energy and the z component of angular momentum remain constant, but the orbital plane precesses. Testing the full momentum vector for conservation in that model would be a mistake.

The integrator uses the fifth-order state and a fourth-order estimate of its local error. Each component is divided by its own absolute tolerance plus a relative tolerance times the larger endpoint magnitude. The largest scaled error determines acceptance. This is a local error control rule, not a promised bound on accumulated trajectory error.

The default tolerances are 1e-8 km for position, 1e-11 km/s for velocity, and 1e-11 relative. Maximum step size is 60 seconds. These values were checked against elliptic trajectories and the included conservation tests; they are not universal choices for every force model or orbit. If the minimum step or work budget is reached, the call fails rather than reporting an unfinished trajectory as complete.

The station model uses the WGS84 ellipsoid and geodetic latitude. Its surface normal is not generally the same as the direction from Earth's centre. Earth rotation is a constant rotation about z. The initial Greenwich angle is explicit in the library API, so an application can state which orientation it assumes.

Pass detection samples elevation and refines sign-changing intervals. It deliberately makes no claim about an interval that never brackets a sign change. Halving the sampling interval and comparing the pass list is one practical check; more demanding applications would need event-completeness bounds, an impact model, and a higher-fidelity Earth-orientation model.
