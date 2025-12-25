#include "orbit/orbit.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <string>

namespace {
constexpr double radians = std::numbers::pi / 180;
double number(const char *text) {
    const std::string value(text);
    std::size_t consumed{};
    const double result = std::stod(value, &consumed);
    if (consumed != value.size() || !std::isfinite(result)) {
        throw std::invalid_argument("expected a finite number: " + value);
    }
    return result;
}
orbit::Gravity model(const std::string &name) {
    orbit::Gravity gravity;
    if (name == "two-body")
        gravity.j2 = 0;
    else if (name != "j2")
        throw std::invalid_argument("model must be two-body or j2");
    return gravity;
}
orbit::State initial_state(char **args, int count) {
    if (count == 0) {
        const double speed = std::sqrt(orbit::Gravity{}.mu / 7000.0);
        return {{7000, 0, 0},
                {0, speed * std::cos(51.6 * radians), speed * std::sin(51.6 * radians)}};
    }
    if (count != 6)
        throw std::invalid_argument("state requires x y z vx vy vz");
    return {{number(args[0]), number(args[1]), number(args[2])},
            {number(args[3]), number(args[4]), number(args[5])}};
}
void row(double time, orbit::State state, orbit::Gravity gravity) {
    const auto &p = state.position;
    const auto &v = state.velocity;
    std::cout << time << ',' << p.x << ',' << p.y << ',' << p.z << ',' << v.x << ',' << v.y << ','
              << v.z << ',' << orbit::specific_energy(state, gravity) << ',' << orbit::cross(p, v).z
              << '\n';
}
constexpr const char *usage =
    "Usage:\n"
    "  orbit-access propagate SECONDS SAMPLE_SECONDS two-body|j2 [x y z vx vy vz]\n"
    "  orbit-access passes SECONDS LAT_DEG LON_DEG MASK_DEG two-body|j2 [x y z vx vy vz]\n"
    "Positions: km. Velocities: km/s. Output: CSV.\n"
    "Default state: 7000-km radius, 51.6-degree inclined circular two-body orbit.\n"
    "Station height and Greenwich angle at t=0 are zero; time is elapsed seconds, not UTC.\n";
} // namespace

int main(int argc, char **argv) {
    if (argc == 2 && std::string(argv[1]) == "--help") {
        std::cout << usage;
        return 0;
    }
    try {
        std::cout << std::setprecision(15);
        if (argc >= 5 && std::string(argv[1]) == "propagate") {
            const double duration = number(argv[2]);
            const double sample = number(argv[3]);
            const auto gravity = model(argv[4]);
            const orbit::Propagator propagator(gravity);
            auto state = initial_state(argv + 5, argc - 5);
            if (!(duration > 0 && sample > 0) || std::ceil(duration / sample) > 1'000'000) {
                throw std::invalid_argument(
                    "positive duration/sample required, at most one million rows");
            }
            (void)propagator.advance(state, 0);
            std::cout << "time_s,x_km,y_km,z_km,vx_km_s,vy_km_s,vz_km_s,energy_km2_s2,hz_km2_s\n";
            row(0, state, gravity);
            for (double time = 0; time < duration;) {
                const double next = std::min(duration, time + sample);
                if (next == time)
                    throw std::runtime_error("sample interval cannot advance clock");
                state = propagator.advance(state, next - time).state;
                row(next, state, gravity);
                time = next;
            }
        } else if (argc >= 7 && std::string(argv[1]) == "passes") {
            orbit::AccessOptions options;
            options.duration = number(argv[2]);
            options.minimum_elevation = number(argv[5]) * radians;
            const orbit::Station station{number(argv[3]) * radians, number(argv[4]) * radians, 0};
            const orbit::Propagator propagator(model(argv[6]));
            const auto state = initial_state(argv + 7, argc - 7);
            const auto passes = orbit::find_passes(propagator, state, station, options);
            std::cout << "rise_s,set_s,duration_s,clipped_at_start,clipped_at_end\n";
            for (const auto &pass : passes) {
                std::cout << pass.rise << ',' << pass.set << ',' << pass.set - pass.rise << ','
                          << pass.clipped_at_start << ',' << pass.clipped_at_end << '\n';
            }
        } else {
            std::cerr << usage;
            return 2;
        }
        if (!std::cout)
            throw std::runtime_error("failed writing output");
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "orbit-access: " << error.what() << '\n';
        return 1;
    }
}
