const std = @import("std");
const zigode = @import("zigode");
const bl = @import("./boyer-lindquist.zig");
const render = @import("./render.zig");

// type aliases
pub const BoyerLindquist = render.BoyerLindquist;
pub const Solver = render.Solver;
pub const Integrator = render.Integrator;

const disc_inner_radius = 10.0;
const disc_outer_radius = 50.0;

pub fn problem(du: *[4]f64, u: *const [4]f64, _: f64, p: *BoyerLindquist) void {
    const v: [4]f64 = p.fourVelocity(u);
    for (du) |*a, i| {
        a.* = v[i];
    }
}

pub fn callback(solver: *Solver, u: *const [4]f64, _: f64, p: *BoyerLindquist) void {
    // check we're not about to hit the event horizon or at effective infinity
    if ((p.inner_horizon * 1.05) > u[1] or u[1] > 1100) {
        solver.terminate();
    } else {
        // update the signs of the potentials
        if (p.radial_potential(u) < 0.0) {
            p.positive_radial_potential = !p.positive_radial_potential;
        }
        if (p.angular_potential(u) < 0.0) {
            p.positive_angular_potential = !p.positive_angular_potential;
        }

        // check intersection with disc
        if (u[1] < disc_outer_radius and u[1] > disc_inner_radius) {
            const uprev = solver.uprev;
            const s1 = u[1] * @cos(u[2]);
            const s2 = uprev[1] * @cos(uprev[2]);
            if (s1 * s2 < 0.0) {
                solver.terminate();
            }
        }
    }
}

pub fn run(allocator: std.mem.Allocator) !void {
    try render.run(allocator, problem, callback, .{}, "disc.data");
}
