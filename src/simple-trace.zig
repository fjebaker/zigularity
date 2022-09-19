const std = @import("std");
const zigode = @import("zigode");
const bl = @import("./boyer-lindquist.zig");

// type aliases
pub const BoyerLindquist = bl.BoyerLindquist(f64);
pub const Solver = zigode.Solver(f64, 4, BoyerLindquist);
pub const Integrator = zigode.AdaptiveTsit5(f64, 4, BoyerLindquist);

pub fn callback(solver: *Solver, u: *const [4]f64, _: f64, p: *BoyerLindquist) void {
    // check we're not about to hit the event horizon or at effective infinity
    if ((p.inner_horizon * 1.05) > u[1]) {
        solver.terminate();
    } else {
        // update the signs of the potentials
        if (p.radial_potential(u) < 0.0) {
            p.positive_radial_potential = !p.positive_radial_potential;
        }
        if (p.angular_potential(u) < 0.0) {
            p.positive_angular_potential = !p.positive_angular_potential;
        }
    }
}

pub fn problem(du: *[4]f64, u: *const [4]f64, _: f64, p: *BoyerLindquist) void {
    const v: [4]f64 = p.fourVelocity(u);
    for (du) |*a, i| {
        a.* = v[i];
    }
}

pub fn run(allocator: std.mem.Allocator) !void {
    // observer position
    const obs_r = 1000.0;
    const obs_theta = std.math.degreesToRadians(f64, 85.0);

    const alpha = 7.0;
    const beta = 0.0;

    // black hole parameters
    const mass = 1.0;
    const spin = 0.998;

    // initial position four vector
    const u: [4]f64 = .{ 0.0, obs_r, obs_theta, 0.0 };
    const p = BoyerLindquist.initImpactParameters(mass, spin, alpha, beta, obs_r, obs_theta);
    var prob = Integrator.init(problem, p);
    var solver = prob.solver(allocator);
    var sol = try solver.solve(
        u,
        0.0,
        2000.0,
        .{ .save = true, .dt = 1.0, .max_iters = 10_000, .callback = callback },
    );
    defer sol.deinit();

    const stdout = std.io.getStdErr();
    try sol.printInfo(stdout);

    var file = try std.fs.cwd().openFile("adaptive-trace.txt", .{ .mode = std.fs.File.OpenMode.write_only });
    defer file.close();

    for (sol.u) |*v| {
        try file.writer().print("{e}\n", .{v.*});
    }
}
