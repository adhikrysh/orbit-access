#include "orbit/orbit.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

namespace orbit {
Vec3 Vec3::operator+(Vec3 b) const {
    return {x + b.x, y + b.y, z + b.z};
}
Vec3 Vec3::operator-(Vec3 b) const {
    return {x - b.x, y - b.y, z - b.z};
}
Vec3 Vec3::operator*(double s) const {
    return {x * s, y * s, z * s};
}
Vec3 Vec3::operator/(double s) const {
    return {x / s, y / s, z / s};
}
double dot(Vec3 a, Vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
Vec3 cross(Vec3 a, Vec3 b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
double norm(Vec3 a) {
    return std::hypot(a.x, a.y, a.z);
}
bool finite(Vec3 a) {
    return std::isfinite(a.x) && std::isfinite(a.y) && std::isfinite(a.z);
}

namespace {
void validate(Gravity g) {
    if (!(std::isfinite(g.mu) && g.mu > 0 && std::isfinite(g.radius) && g.radius > 0 &&
          std::isfinite(g.j2) && std::abs(g.j2) < 0.1)) {
        throw std::invalid_argument("gravity requires finite positive mu/radius and |J2| < 0.1");
    }
}
Vec3 force(Vec3 p, Gravity g) {
    const double r = norm(p);
    if (!finite(p) || !std::isfinite(r) || r <= g.radius) {
        throw std::domain_error("state is nonfinite or inside the gravity reference radius");
    }
    const double z2 = (p.z / r) * (p.z / r);
    const double central = -g.mu / (r * r * r);
    const double j = 1.5 * g.j2 * g.mu / (r * r * r) * (g.radius / r) * (g.radius / r);
    return {central * p.x + j * p.x * (5 * z2 - 1), central * p.y + j * p.y * (5 * z2 - 1),
            central * p.z + j * p.z * (5 * z2 - 3)};
}
using Six = std::array<double, 6>;
Six pack(State s) {
    return {s.position.x, s.position.y, s.position.z, s.velocity.x, s.velocity.y, s.velocity.z};
}
State unpack(const Six &y) {
    return {{y[0], y[1], y[2]}, {y[3], y[4], y[5]}};
}
Six derivative(const Six &y, Gravity gravity) {
    const Vec3 a = force({y[0], y[1], y[2]}, gravity);
    return {y[3], y[4], y[5], a.x, a.y, a.z};
}
// Dormand-Prince 5(4). The fifth-order state advances the solution; its difference
// from the fourth-order state estimates local error for step-size control.
constexpr double stages[7][7] = {
    {},
    {1.0 / 5},
    {3.0 / 40, 9.0 / 40},
    {44.0 / 45, -56.0 / 15, 32.0 / 9},
    {19372.0 / 6561, -25360.0 / 2187, 64448.0 / 6561, -212.0 / 729},
    {9017.0 / 3168, -355.0 / 33, 46732.0 / 5247, 49.0 / 176, -5103.0 / 18656},
    {35.0 / 384, 0, 500.0 / 1113, 125.0 / 192, -2187.0 / 6784, 11.0 / 84}};
constexpr double low[7] = {5179.0 / 57600, 0,       7571.0 / 16695, 393.0 / 640, -92097.0 / 339200,
                           187.0 / 2100,   1.0 / 40};
} // namespace

Vec3 acceleration(Vec3 position, Gravity g) {
    validate(g);
    return force(position, g);
}

double specific_energy(State state, Gravity g) {
    validate(g);
    (void)force(state.position, g);
    if (!finite(state.velocity))
        throw std::invalid_argument("nonfinite velocity");
    const double r = norm(state.position);
    const double z = state.position.z / r;
    const double potential =
        -g.mu / r * (1 - 0.5 * g.j2 * (g.radius / r) * (g.radius / r) * (3 * z * z - 1));
    return 0.5 * dot(state.velocity, state.velocity) + potential;
}

Propagator::Propagator(Gravity gravity, Tolerances tolerances)
    : gravity_(gravity), tolerances_(tolerances) {
    validate(gravity);
    const auto &o = tolerances_;
    for (double value :
         {o.position_abs, o.velocity_abs, o.relative, o.initial_step, o.min_step, o.max_step}) {
        if (!(std::isfinite(value) && value > 0)) {
            throw std::invalid_argument(
                "integration tolerances and step sizes must be finite and positive");
        }
    }
    if (o.min_step > o.max_step || o.max_attempts == 0) {
        throw std::invalid_argument("invalid integration step bounds or work budget");
    }
}

Propagation Propagator::advance(State initial, double dt) const {
    if (!finite(initial.velocity) || !std::isfinite(dt)) {
        throw std::invalid_argument("state and elapsed time must be finite");
    }
    (void)force(initial.position, gravity_);
    if (dt == 0)
        return {initial, {}};
    const auto &o = tolerances_;
    Six y = pack(initial);
    const double direction = std::copysign(1.0, dt);
    double h = direction * std::clamp(o.initial_step, o.min_step, o.max_step);
    double time = 0;
    IntegrationStats stats;
    for (std::size_t attempt = 0; attempt < o.max_attempts; ++attempt) {
        h = direction * std::min(std::abs(h), std::abs(dt - time));
        if (time + h == time)
            throw std::runtime_error("integration clock cannot advance");
        std::array<Six, 7> k{};
        Six fifth{};
        for (std::size_t stage = 0; stage < 7; ++stage) {
            Six trial = y;
            for (std::size_t component = 0; component < 6; ++component) {
                for (std::size_t earlier = 0; earlier < stage; ++earlier) {
                    trial[component] += h * stages[stage][earlier] * k[earlier][component];
                }
            }
            k[stage] = derivative(trial, gravity_);
            ++stats.force_evaluations;
            if (stage == 6)
                fifth = trial;
        }
        double error = 0;
        for (std::size_t i = 0; i < 6; ++i) {
            double fourth = y[i];
            for (std::size_t stage = 0; stage < 7; ++stage)
                fourth += h * low[stage] * k[stage][i];
            const double absolute = i < 3 ? o.position_abs : o.velocity_abs;
            const double scale =
                absolute + o.relative * std::max(std::abs(y[i]), std::abs(fifth[i]));
            error = std::max(error, std::abs(fifth[i] - fourth) / scale);
            if (!std::isfinite(fifth[i]) || !std::isfinite(fourth)) {
                throw std::runtime_error("nonfinite integration state");
            }
        }
        const double factor = error == 0 ? 5.0 : std::clamp(0.9 * std::pow(error, -0.2), 0.2, 5.0);
        if (error <= 1) {
            time += h;
            y = fifth;
            ++stats.accepted;
            if (direction * (dt - time) <= 0)
                return {unpack(y), stats};
        } else {
            ++stats.rejected;
            // A final fractional step may be shorter than min_step. Accept it
            // only if its error passes; never override a failed error estimate.
            if (std::abs(h) <= o.min_step) {
                throw std::runtime_error("requested accuracy requires a smaller minimum step");
            }
        }
        h = direction * std::clamp(std::abs(h) * factor, o.min_step, o.max_step);
    }
    throw std::runtime_error("integration work budget exhausted before the requested time");
}
} // namespace orbit
