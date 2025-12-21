#include "orbit/orbit.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace orbit {
namespace {
constexpr double pi = std::numbers::pi;
constexpr double wgs84_a = 6378.137;
constexpr double wgs84_f = 1.0 / 298.257223563;
void validate(Station s, EarthRotation r) {
    if (!(std::isfinite(s.latitude) && std::abs(s.latitude) <= pi / 2 &&
          std::isfinite(s.longitude) && std::abs(s.longitude) <= pi && std::isfinite(s.height) &&
          s.height >= 0 && std::isfinite(r.angle_at_zero) && std::isfinite(r.rate))) {
        throw std::invalid_argument("invalid station coordinates or Earth rotation model");
    }
}
Vec3 rotate(Vec3 v, double angle) {
    return {std::cos(angle) * v.x - std::sin(angle) * v.y,
            std::sin(angle) * v.x + std::cos(angle) * v.y, v.z};
}
double earth_angle(double t, EarthRotation r) {
    if (!std::isfinite(t))
        throw std::invalid_argument("nonfinite observation time");
    const double angle = r.angle_at_zero + r.rate * t;
    if (!std::isfinite(angle))
        throw std::invalid_argument("Earth rotation angle overflow");
    return std::remainder(angle, 2 * pi);
}
} // namespace

Vec3 station_position(Station s, double t, EarthRotation r) {
    validate(s, r);
    const double e2 = wgs84_f * (2 - wgs84_f);
    const double sin_lat = std::sin(s.latitude);
    const double n = wgs84_a / std::sqrt(1 - e2 * sin_lat * sin_lat);
    return rotate({(n + s.height) * std::cos(s.latitude) * std::cos(s.longitude),
                   (n + s.height) * std::cos(s.latitude) * std::sin(s.longitude),
                   (n * (1 - e2) + s.height) * sin_lat},
                  earth_angle(t, r));
}

double elevation(Vec3 satellite, Station s, double t, EarthRotation r) {
    const Vec3 line = satellite - station_position(s, t, r);
    const double range = norm(line);
    if (!finite(satellite) || !(range > 0) || !std::isfinite(range)) {
        throw std::invalid_argument("invalid satellite-to-station line of sight");
    }
    const Vec3 up = rotate({std::cos(s.latitude) * std::cos(s.longitude),
                            std::cos(s.latitude) * std::sin(s.longitude), std::sin(s.latitude)},
                           earth_angle(t, r));
    return std::asin(std::clamp(dot(line / range, up), -1.0, 1.0));
}

std::vector<Pass> find_passes(const Propagator &propagator, State initial, Station station,
                              AccessOptions o, EarthRotation rotation) {
    validate(station, rotation);
    if (!(std::isfinite(o.duration) && o.duration > 0 && std::isfinite(o.sample_step) &&
          o.sample_step > 0 && std::isfinite(o.root_tolerance) && o.root_tolerance > 0 &&
          o.root_tolerance < o.sample_step && std::isfinite(o.minimum_elevation) &&
          o.minimum_elevation >= 0 && o.minimum_elevation < pi / 2 && o.max_samples > 0)) {
        throw std::invalid_argument("invalid pass-search options");
    }
    if (std::ceil(o.duration / o.sample_step) > static_cast<double>(o.max_samples)) {
        throw std::invalid_argument("pass search exceeds its sample budget");
    }
    State state = propagator.advance(initial, 0).state;
    double time = 0;
    bool visible = elevation(state.position, station, time, rotation) > o.minimum_elevation;
    Pass open{0, 0, visible, false};
    std::vector<Pass> passes;
    while (time < o.duration) {
        const double next_time = std::min(o.duration, time + o.sample_step);
        if (next_time == time)
            throw std::runtime_error("pass-search clock cannot advance");
        const State next = propagator.advance(state, next_time - time).state;
        const bool next_visible =
            elevation(next.position, station, next_time, rotation) > o.minimum_elevation;
        if (next_visible != visible) {
            double left = time;
            double right = next_time;
            // Propagate from the same left sample at every midpoint. Updating
            // the bracket must not silently change the epoch of the base state.
            for (unsigned iteration = 0; right - left > o.root_tolerance; ++iteration) {
                if (iteration == 100)
                    throw std::runtime_error("pass-root refinement exhausted");
                const double middle = left + (right - left) / 2;
                if (middle == left || middle == right)
                    throw std::runtime_error("pass-root precision exhausted");
                const State trial = propagator.advance(state, middle - time).state;
                const bool mid_visible =
                    elevation(trial.position, station, middle, rotation) > o.minimum_elevation;
                if (mid_visible == visible)
                    left = middle;
                else
                    right = middle;
            }
            const double crossing = left + (right - left) / 2;
            if (next_visible)
                open = {crossing, 0, false, false};
            else {
                open.set = crossing;
                if (open.set > open.rise)
                    passes.push_back(open);
            }
        }
        visible = next_visible;
        state = next;
        time = next_time;
    }
    if (visible) {
        open.set = o.duration;
        open.clipped_at_end = true;
        if (open.set > open.rise)
            passes.push_back(open);
    }
    return passes;
}
} // namespace orbit
