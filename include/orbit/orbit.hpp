#pragma once

#include <cstddef>
#include <vector>

namespace orbit {

// Distances are kilometres, elapsed time is seconds, and angles are radians.
// The inertial frame is right-handed, with Earth's spin along +z.
struct Vec3 {
    double x{}, y{}, z{};
    Vec3 operator+(Vec3 b) const;
    Vec3 operator-(Vec3 b) const;
    Vec3 operator*(double scalar) const;
    Vec3 operator/(double scalar) const;
};
double dot(Vec3 a, Vec3 b);
Vec3 cross(Vec3 a, Vec3 b);
double norm(Vec3 a);
bool finite(Vec3 a);

struct State {
    Vec3 position; // km in the inertial frame
    Vec3 velocity; // km/s in the same frame
};

struct Gravity {
    double mu{398600.4418};   // km^3/s^2
    double radius{6378.137};  // km; also the forbidden interior boundary
    double j2{1.08262668e-3}; // dimensionless; set to zero for two-body motion
};

Vec3 acceleration(Vec3 position, Gravity gravity = {});
double specific_energy(State state, Gravity gravity = {}); // km^2/s^2

struct Tolerances {
    double position_abs{1e-8};  // km, separate from the velocity scale
    double velocity_abs{1e-11}; // km/s
    double relative{1e-11};
    double initial_step{10.0};
    double min_step{1e-7};
    double max_step{60.0};
    std::size_t max_attempts{1'000'000};
};

struct IntegrationStats {
    std::size_t accepted{};
    std::size_t rejected{};
    std::size_t force_evaluations{};
};
struct Propagation {
    State state;
    IntegrationStats stats;
};

class Propagator {
  public:
    explicit Propagator(Gravity gravity = {}, Tolerances tolerances = {});
    // dt may be negative. Throws on invalid input, an interior force evaluation,
    // underflow, or exhausted work budget. It never returns a partial trajectory
    // as though it had reached the requested time.
    [[nodiscard]] Propagation advance(State initial, double dt) const;
    [[nodiscard]] Gravity gravity() const {
        return gravity_;
    }

  private:
    Gravity gravity_;
    Tolerances tolerances_;
};

struct Station {
    double latitude{};  // geodetic radians, [-pi/2, pi/2]
    double longitude{}; // radians, [-pi, pi]
    double height{};    // km above the WGS84 ellipsoid, >= 0
};
struct EarthRotation {
    double angle_at_zero{};    // Greenwich inertial angle at elapsed t=0
    double rate{7.2921150e-5}; // rad/s; constant spin, no polar motion
};

// WGS84 station geometry; the gravity reference radius is a separate parameter.
Vec3 station_position(Station station, double time, EarthRotation rotation = {});
double elevation(Vec3 satellite, Station station, double time, EarthRotation rotation = {});

struct Pass {
    double rise{};
    double set{};
    bool clipped_at_start{};
    bool clipped_at_end{};
};
struct AccessOptions {
    double duration{86400.0};
    double sample_step{10.0};
    double root_tolerance{1e-5};
    double minimum_elevation{0.17453292519943295}; // 10 degrees
    std::size_t max_samples{1'000'000};
};

// initial is the state at time zero. Sign-changing rise/set brackets are refined
// with bisection and propagation, not linear interpolation of elevation.
// A pass shorter than sample_step may be missed; root precision is not a claim
// of complete detection. Tangencies without an above-mask interval are omitted.
std::vector<Pass> find_passes(const Propagator &propagator, State initial, Station station,
                              AccessOptions options = {}, EarthRotation rotation = {});

} // namespace orbit
