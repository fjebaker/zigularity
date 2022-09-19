const std = @import("std");

const pow = std.math.pow;

// Bardeen et al. (1972) Eq. (2.10)
fn bl_T(comptime T: type, E: T, L: T, r: T, a: T) T {
    return E * (pow(T, r, 2) + pow(T, a, 2)) - L * a;
}

// Bardeen et al. (1972) Eq. (2.10), special case of mu = 0
// Effective potential in r
fn bl_Vr(comptime T: type, E: T, L: T, M: T, Q: T, r: T, a: T) T {
    const t = bl_T(T, E, L, r, a);
    const delta = bl_Delta(T, M, r, a);
    return pow(T, t, 2) - delta * (pow(T, L - a * E, 2) + Q);
}

test "bl_Vr" {
    const vr = bl_Vr(f64, 1.0, -9.96194698091745, 1.0, 0.7520465799607337, 14.0, 0.998);
    try std.testing.expectApproxEqAbs(vr, 22396.380196731767, 1e-5);
}

// Effective potential in theta
fn bl_Vtheta(comptime T: type, E: T, L: T, Q: T, a: T, theta: T) T {
    const cos_theta_2 = pow(T, @cos(theta), 2);
    const sin_theta = @sin(theta);
    return Q + cos_theta_2 * (pow(T, a * E, 2) - pow(T, L / sin_theta, 2));
}

test "bl_Vtheta" {
    const vtheta = bl_Vtheta(
        f64,
        1.0,
        -9.96194698091745,
        0.7520465799607337,
        0.998,
        std.math.degreesToRadians(f64, 85.0),
    );
    try std.testing.expectApproxEqAbs(vtheta, -4.444766776856568e-11, 1e-5);
}

// Bardeen et al. (1972) Eq. (2.9)
// Differential terms in t, r, theta, phi with respect to integration parameter s
fn dt_ds(comptime T: type, E: T, L: T, M: T, r: T, a: T, theta: T) T {
    const t = bl_T(T, E, L, r, a);
    const delta = bl_Delta(T, M, r, a);
    const sin_theta = @sin(theta);

    return -a * (a * E * pow(T, sin_theta, 2) - L) +
        (pow(T, r, 2) + pow(T, a, 2) * t / delta);
}
fn dr_ds(comptime T: type, E: T, L: T, M: T, Q: T, r: T, a: T) T {
    const potential = bl_Vr(T, E, L, M, Q, r, a);
    return @sqrt(@fabs(potential));
}
fn dtheta_ds(comptime T: type, E: T, L: T, Q: T, a: T, theta: T) T {
    const potential = bl_Vtheta(T, E, L, Q, a, theta);
    return @sqrt(@fabs(potential));
}
fn dphi_ds(comptime T: type, E: T, L: T, M: T, r: T, a: T, theta: T) T {
    const t = bl_T(T, E, L, r, a);
    const delta = bl_Delta(T, M, r, a);
    const sin_theta = @sin(theta);

    return (L / pow(T, sin_theta, 2)) - (a * E) + a * t / delta;
}

// Fanton et al. (1997), Eq. (76)
fn fa_S(comptime T: type, theta: T, alpha: T, omega: T) T {
    return 1 + alpha * omega * @sin(theta);
}

pub fn ObserverAngles(comptime T: type) type {
    return struct { angle1: T, angle2: T };
}

// Fanton et al. (1997), Eq. (74-75)
// Map impact parameters to observer angles
fn observerAngles(comptime T: type, M: T, r: T, a: T, theta: T, alpha: T, beta: T) ObserverAngles(T) {
    const sigma = bl_Sigma(T, r, a, theta);
    const sin_theta = @sin(theta);
    const A = bl_A(T, M, r, a, theta);
    const delta = bl_Delta(T, M, r, a);
    const omega = 2.0 * a * r / A;
    const S = fa_S(T, theta, alpha, omega);
    return _implObserverAngles(T, sigma, sin_theta, A, delta, S, r, a, alpha, beta);
}
fn _implObserverAngles(
    comptime T: type,
    sigma: T,
    sin_theta: T,
    A: T,
    delta: T,
    S: T,
    r: T,
    a: T,
    alpha: T,
    beta: T,
) ObserverAngles(T) {
    const denom1 = pow(T, beta, 2) +
        pow(T, alpha + a * sin_theta, 2) +
        (A * pow(T, S, 2) - pow(T, (pow(T, r, 2) + pow(T, a, 2) + a * alpha * sin_theta), 2)) / delta;

    const angle1 = ((alpha * sigma) / @sqrt(A)) / @sqrt(denom1);
    const angle2 = -(alpha * sigma * @sqrt(delta)) / (S * A * angle1);

    return .{ .angle1 = angle1, .angle2 = angle2 };
}

test "observerAngles" {
    const obs_angles = observerAngles(f64, 1.0, 1000.0, 0.998, std.math.degreesToRadians(f64, 85), 10.0, 0.0);
    try std.testing.expectApproxEqAbs(obs_angles.angle1, 0.999999999, 1e-5);
    try std.testing.expectApproxEqAbs(obs_angles.angle2, -0.00998998, 1e-5);
}

pub fn ConstantsOfMotion(comptime T: type) type {
    return struct { angular_momentum: T, carters_q: T };
}

// Fanton et al. Eq. (69-73)
fn constantsOfMotion(
    comptime T: type,
    M: T,
    r: T,
    a: T,
    theta: T,
    obs_angles: ObserverAngles(T),
) ConstantsOfMotion(T) {
    const sin_theta = @sin(theta);
    const sigma = bl_Sigma(T, r, a, theta);
    const A = bl_A(T, M, r, a, theta);
    const delta = bl_Delta(T, M, r, a);
    const omega = 2.0 * a * r / A;

    const upsilon1 = sin_theta * obs_angles.angle1 * obs_angles.angle2;
    const upsilon2 = (sigma * @sqrt(delta) / A) + omega * upsilon1;

    const angmom = upsilon1 / upsilon2;
    const p = (pow(T, r, 2) + pow(T, a, 2) - a * angmom);
    const carters_q = (pow(T, p, 2) / delta) -
        pow(T, angmom - a, 2) -
        ((pow(T, sigma, 2) / A) * pow(T, @cos(std.math.asin(obs_angles.angle2)) / upsilon2, 2));
    return .{ .angular_momentum = angmom, .carters_q = carters_q };
}

test "constantsOfMotion" {
    const lq = constantsOfMotion(
        f64,
        1.0,
        1000.0,
        0.998,
        std.math.degreesToRadians(f64, 85.0),
        .{ .angle1 = 0.9999999999988463, .angle2 = -0.00998998981155132 },
    );
    try std.testing.expectApproxEqAbs(lq.angular_momentum, -9.96194698091745, 1e-5);
    try std.testing.expectApproxEqAbs(lq.carters_q, 0.7520465799607337, 1e-5);
}

// Below are all metric parameters
// Bardeen et al. (1972) Eq. (2.3)
fn bl_Sigma(comptime T: type, r: T, a: T, theta: T) T {
    return pow(T, r, 2) + pow(T, a * @cos(theta), 2);
}

test "bl_Sigma" {
    const sigma = bl_Sigma(f64, 1.0, 17.0, 0.998);
    try std.testing.expectApproxEqAbs(sigma, 85.8928356952155, 1e-5);
}

fn bl_Delta(comptime T: type, M: T, r: T, a: T) T {
    return pow(T, r, 2) - 2 * M * r + pow(T, a, 2);
}

test "bl_Delta" {
    const delta = bl_Delta(f64, 1.0, 17.0, 0.998);
    try std.testing.expectApproxEqAbs(delta, 255.996004, 1e-5);
}

fn bl_A(comptime T: type, M: T, r: T, a: T, theta: T) T {
    const delta = bl_Delta(T, M, r, a);
    const sin_theta = @sin(theta);

    return pow(T, pow(T, r, 2) + pow(T, a, 2), 2) - pow(T, a, 2) * delta * pow(T, sin_theta, 2);
}

test "bl_A" {
    const a = bl_A(f64, 1.0, 17.0, 0.998, std.math.degreesToRadians(f64, 85.0));
    try std.testing.expectApproxEqAbs(a, 83844.6460987296, 1e-5);
}

fn innerHoziron(comptime T: type, M: T, a: T) T {
    return M + @sqrt(pow(T, M, 2) - pow(T, a, 2));
}

test "innerHorizon" {
    const hor = innerHoziron(f64, 1.0, 0.998);
    try std.testing.expectApproxEqAbs(hor, 1.0632139225171164, 1e-5);
}

pub fn BoyerLindquist(comptime T: type) type {
    return struct {
        const Self = @This();

        mass: T,
        spin: T,
        angular_momentum: T,
        carters_q: T,
        // no reason why you wouldn't want this to be 1.0
        energy: T = 1.0,
        positive_radial_potential: bool = false,
        positive_angular_potential: bool = true,
        inner_horizon: T,

        pub fn init(mass: T, spin: T, angmom: T, carters_q: T) Self {
            const inner_horizon = innerHoziron(T, mass, spin);
            return .{
                .mass = mass,
                .spin = spin,
                .angular_momentum = angmom,
                .carters_q = carters_q,
                .inner_horizon = inner_horizon,
            };
        }

        pub fn initImpactParameters(
            mass: T,
            spin: T,
            alpha: T,
            beta: T,
            obs_r: T,
            obs_theta: T,
        ) Self {
            const obs_angles = observerAngles(T, mass, obs_r, spin, obs_theta, alpha, beta);
            const lq = constantsOfMotion(T, mass, obs_r, spin, obs_theta, obs_angles);
            var self = Self.init(mass, spin, lq.angular_momentum, lq.carters_q);
            if (beta >= 0) {
                self.positive_angular_potential = true;
            } else {
                self.positive_angular_potential = false;
            }
            return self;
        }

        pub fn fourVelocity(self: *Self, u: *const [4]T) [4]T {
            // bindings purely so i don't make a mistake
            // whilst writing out these equations
            const E = self.energy;
            const Q = self.carters_q;
            const L = self.angular_momentum;
            const M = self.mass;
            const a = self.spin;
            const r = u[1];
            const theta = u[2];

            const sigma = bl_Sigma(T, r, a, theta);

            const vt = dt_ds(T, E, L, M, r, a, theta) / sigma;
            const vr = dr_ds(T, E, L, M, Q, r, a) / sigma;
            const vtheta = dtheta_ds(T, E, L, Q, a, theta) / sigma;
            const vphi = dphi_ds(T, E, L, M, r, a, theta) / sigma;

            return .{
                vt,
                if (self.positive_radial_potential) vr else -vr,
                if (self.positive_angular_potential) vtheta else -vtheta,
                vphi,
            };
        }

        pub fn radial_potential(self: *Self, u: *const [4]T) T {
            return bl_Vr(
                T,
                self.energy,
                self.angular_momentum,
                self.mass,
                self.carters_q,
                u[1],
                self.spin,
            );
        }

        pub fn angular_potential(self: *Self, u: *const [4]T) T {
            return bl_Vtheta(
                T,
                self.energy,
                self.angular_momentum,
                self.carters_q,
                self.spin,
                u[2],
            );
        }
    };
}
