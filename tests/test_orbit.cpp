#include "orbit/orbit.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <numbers>
#include <random>
#include <stdexcept>
#include <string>

namespace {
using namespace orbit;
constexpr double pi = std::numbers::pi;
constexpr double mu = 398600.4418;
const Gravity two_body{mu, 6378.137, 0};

void require(bool condition, const std::string &reason) {
    if (!condition)
        throw std::runtime_error(reason);
}
void near(double actual, double expected, double tolerance, const std::string &reason) {
    if (!std::isfinite(actual) || std::abs(actual - expected) > tolerance) {
        throw std::runtime_error(reason + ": actual=" + std::to_string(actual) +
                                 " expected=" + std::to_string(expected));
    }
}
template <class F> void rejects(F action) {
    bool threw = false;
    try {
        action();
    } catch (const std::exception &) {
        threw = true;
    }
    require(threw, "invalid request was accepted");
}
State circle(double radius = 7000, double inclination = 0) {
    const double v = std::sqrt(mu / radius);
    return {{radius, 0, 0}, {0, v * std::cos(inclination), v * std::sin(inclination)}};
}

// This reference solves Kepler's algebraic equation, not the numerical ODE.
// Starting at periapsis makes the independent initial condition unambiguous.
State kepler(double a, double eccentricity, double dt) {
    const double n = std::sqrt(mu / (a * a * a));
    const double mean = std::remainder(n * dt, 2 * pi);
    double anomaly = mean;
    for (unsigned i = 0; i < 100; ++i) {
        const double residual = anomaly - eccentricity * std::sin(anomaly) - mean;
        const double step = residual / (1 - eccentricity * std::cos(anomaly));
        anomaly -= step;
        if (std::abs(step) < 1e-14)
            break;
        require(i < 99, "reference Kepler solver failed to converge");
    }
    const double b = std::sqrt(1 - eccentricity * eccentricity);
    const double denominator = 1 - eccentricity * std::cos(anomaly);
    return {
        {a * (std::cos(anomaly) - eccentricity), a * b * std::sin(anomaly), 0},
        {-a * n * std::sin(anomaly) / denominator, a * n * b * std::cos(anomaly) / denominator, 0}};
}

void kepler_ensemble() {
    std::mt19937_64 rng(90210);
    std::uniform_real_distribution<double> unit(0, 1);
    const Propagator propagator(two_body);
    for (int i = 0; i < 150; ++i) {
        const double e = 0.75 * unit(rng);
        const double a = (6600 + 5000 * unit(rng)) / (1 - e);
        const double period = 2 * pi * std::sqrt(a * a * a / mu);
        const double dt = (4 * unit(rng) - 2) * period;
        const State initial = kepler(a, e, 0);
        const State expected = kepler(a, e, dt);
        const auto result = propagator.advance(initial, dt);
        require(norm(result.state.position - expected.position) < 5e-4,
                "elliptic propagation differs from Kepler by more than 0.5 m");
        require(norm(result.state.velocity - expected.velocity) < 5e-7,
                "elliptic velocity differs from Kepler by more than 0.5 mm/s");
        require(result.stats.force_evaluations ==
                    7 * (result.stats.accepted + result.stats.rejected),
                "incorrect integration counters");
    }
}

void conserved_quantities() {
    const State initial = circle(7000, 0.9);
    const double period = 2 * pi * std::sqrt(7000.0 * 7000 * 7000 / mu);
    const State result = Propagator(two_body).advance(initial, 100 * period).state;
    near(specific_energy(result, two_body), specific_energy(initial, two_body), 2e-7,
         "two-body energy drift");
    require(norm(cross(result.position, result.velocity) -
                 cross(initial.position, initial.velocity)) < 2e-4,
            "two-body angular momentum drift");
    const Gravity j2;
    const State perturbed = Propagator(j2).advance(initial, 20 * period).state;
    near(specific_energy(perturbed, j2), specific_energy(initial, j2), 1e-7,
         "J2 potential energy drift");
    near(cross(perturbed.position, perturbed.velocity).z,
         cross(initial.position, initial.velocity).z, 1e-4,
         "axisymmetric J2 must conserve z momentum");
}

void potential_gradient() {
    const Gravity g;
    for (Vec3 p : {Vec3{7000, 1500, 800}, Vec3{-10000, 12000, -5000}, Vec3{0, 0, 8000}}) {
        const double epsilon = 0.01;
        Vec3 gradient;
        gradient.x = (specific_energy({p + Vec3{epsilon, 0, 0}, {}}, g) -
                      specific_energy({p - Vec3{epsilon, 0, 0}, {}}, g)) /
                     (2 * epsilon);
        gradient.y = (specific_energy({p + Vec3{0, epsilon, 0}, {}}, g) -
                      specific_energy({p - Vec3{0, epsilon, 0}, {}}, g)) /
                     (2 * epsilon);
        gradient.z = (specific_energy({p + Vec3{0, 0, epsilon}, {}}, g) -
                      specific_energy({p - Vec3{0, 0, epsilon}, {}}, g)) /
                     (2 * epsilon);
        require(norm(acceleration(p, g) + gradient) < 2e-12,
                "J2 force is not minus potential gradient");
    }
}

void j2_node_drift() {
    const double radius = 7000;
    const double inclination = 0.9;
    const double n = std::sqrt(mu / (radius * radius * radius));
    const double period = 2 * pi / n;
    const Gravity g;
    const double expected_rate =
        -1.5 * g.j2 * n * std::pow(g.radius / radius, 2) * std::cos(inclination);
    auto state = circle(radius, inclination);
    const Propagator propagator(g);
    double sum_t = 0, sum_angle = 0, sum_tt = 0, sum_t_angle = 0;
    constexpr int samples = 200;
    for (int i = 1; i <= samples; ++i) {
        state = propagator.advance(state, period / 10).state;
        const Vec3 node = cross({0, 0, 1}, cross(state.position, state.velocity));
        const double angle = std::atan2(node.y, node.x);
        const double t = static_cast<double>(i) * period / 10;
        sum_t += t;
        sum_angle += angle;
        sum_tt += t * t;
        sum_t_angle += t * angle;
    }
    const double rate =
        (samples * sum_t_angle - sum_t * sum_angle) / (samples * sum_tt - sum_t * sum_t);
    // First-order mean-element theory differs slightly from these osculating
    // initial elements. A 1% envelope checks sign and scale without claiming
    // that the secular approximation is an exact trajectory oracle.
    require(std::abs((rate - expected_rate) / expected_rate) < 0.01,
            "J2 node drift differs from secular theory");
}

void analytic_access() {
    constexpr double radius = 7000;
    constexpr double station_radius = 6378.137;
    const double n = std::sqrt(mu / (radius * radius * radius));
    const EarthRotation rotation;
    const double relative_rate = n - rotation.rate;
    AccessOptions options;
    options.duration = 1.1 * 2 * pi / relative_rate;
    options.sample_step = 31;
    options.minimum_elevation = 0.1;
    options.root_tolerance = 1e-6;
    const auto passes = find_passes(Propagator(two_body), circle(radius), {}, options, rotation);
    // Circular equatorial geometry gives the central angle at a given elevation
    // directly, independently of the event finder and coordinate transforms.
    const double half_angle =
        std::acos(station_radius / radius * std::cos(options.minimum_elevation)) -
        options.minimum_elevation;
    const double half_duration = half_angle / relative_rate;
    const double revisit = 2 * pi / relative_rate;
    require(passes.size() == 2, "analytic equatorial case should have two visible intervals");
    require(passes[0].clipped_at_start && !passes[0].clipped_at_end, "initial clipping flags");
    near(passes[0].rise, 0, 0, "initial rise");
    near(passes[0].set, half_duration, 2e-5, "analytic first set");
    near(passes[1].rise, revisit - half_duration, 2e-5, "analytic next rise");
    near(passes[1].set, revisit + half_duration, 2e-5, "analytic next set");
    options.sample_step = 7;
    const auto finer = find_passes(Propagator(two_body), circle(radius), {}, options, rotation);
    require(finer.size() == passes.size(), "pass count changes with smaller step");
    for (std::size_t i = 0; i < passes.size(); ++i) {
        near(finer[i].rise, passes[i].rise, 2e-5, "rise convergence");
        near(finer[i].set, passes[i].set, 2e-5, "set convergence");
    }
}

void geometry_and_clipping() {
    near(norm(station_position({}, 0)), 6378.137, 1e-10, "WGS84 equator");
    near(station_position({pi / 2, 0, 0}, 0).z, 6356.752314245, 1e-8, "WGS84 pole");
    near(elevation({7000, 0, 0}, {}, 0), pi / 2, 1e-15, "overhead elevation");
    near(elevation({-7000, 0, 0}, {}, 0), -pi / 2, 1e-15, "antipodal elevation");
    AccessOptions options;
    options.duration = 1;
    auto passes = find_passes(Propagator(two_body), circle(), {}, options);
    require(passes.size() == 1 && passes[0].clipped_at_start && passes[0].clipped_at_end,
            "observation inside one pass must retain both clipping flags");
    require(find_passes(Propagator(two_body), circle(), {pi / 2, 0, 0}, options).empty(),
            "polar station sees equatorial spacecraft");
    const Station tangent_station{0.1, 0, 0};
    const EarthRotation fixed_earth{0, 0};
    options.minimum_elevation = elevation(circle().position, tangent_station, 0, fixed_earth);
    require(
        find_passes(Propagator(two_body), circle(), tangent_station, options, fixed_earth).empty(),
        "touching the elevation mask must not create a positive-duration pass");
}

void rejected_steps_and_limits() {
    Tolerances loose;
    loose.relative = 1e-5;
    loose.position_abs = 1e-3;
    loose.velocity_abs = 1e-6;
    loose.initial_step = 1500;
    loose.max_step = 1500;
    const auto initial = kepler(15000, 0.55, 0);
    const auto reference = kepler(15000, 0.55, 8000);
    const auto approximate = Propagator(two_body, loose).advance(initial, 8000);
    const auto precise = Propagator(two_body).advance(initial, 8000);
    require(approximate.stats.rejected > 0, "large initial step should be rejected");
    require(norm(precise.state.position - reference.position) <
                norm(approximate.state.position - reference.position) * 0.01,
            "tightening tolerance should improve accuracy");
    Tolerances capped;
    capped.max_attempts = 1;
    rejects([&] { (void)Propagator(two_body, capped).advance(circle(), 1000); });
    rejects([] { (void)Propagator({-1, 6378, 0}); });
    rejects([] { (void)Propagator(two_body).advance({{}, {}}, 1); });
    rejects([] {
        (void)Propagator(two_body).advance(circle(), std::numeric_limits<double>::infinity());
    });
    rejects([] { (void)station_position({pi, 0, 0}, 0); });
    rejects([] {
        AccessOptions o;
        o.sample_step = 0;
        (void)find_passes(Propagator(two_body), circle(), {}, o);
    });
    rejects([] {
        AccessOptions o;
        o.max_samples = 1;
        (void)find_passes(Propagator(two_body), circle(), {}, o);
    });
    const auto no_time = Propagator(two_body).advance(circle(), 0);
    require(no_time.stats.accepted == 0 && norm(no_time.state.position - circle().position) == 0,
            "zero-time propagation");
}
} // namespace

int main() {
    unsigned failed = 0;
    auto test = [&](const char *name, const std::function<void()> &body) {
        try {
            body();
            std::cout << "PASS " << name << '\n';
        } catch (const std::exception &e) {
            ++failed;
            std::cerr << "FAIL " << name << ": " << e.what() << '\n';
        }
    };
    test("150 independent elliptic Kepler comparisons", kepler_ensemble);
    test("long-run conserved quantities", conserved_quantities);
    test("J2 potential-gradient check", potential_gradient);
    test("J2 secular node drift", j2_node_drift);
    test("analytic equatorial visibility and event convergence", analytic_access);
    test("WGS84 geometry and clipped windows", geometry_and_clipping);
    test("adaptive rejection, input validation, and work budgets", rejected_steps_and_limits);
    return failed == 0 ? 0 : 1;
}
