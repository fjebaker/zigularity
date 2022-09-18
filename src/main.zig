const std = @import("std");
const zigode = @import("zigode");
const bl = @import("./boyer-lindquist.zig");

// type aliases
const BoyerLindquist = bl.BoyerLindquist(f64);
const Solver = zigode.Solver(f64, 4, BoyerLindquist);
const Integrator = zigode.Tsit5(f64, 4, BoyerLindquist);

pub fn callback(solver: *Solver, u: *const [4]f64, _: f64, p: *BoyerLindquist) void {
    // check we're not about to hit the event horizon
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

var gpa = std.heap.GeneralPurposeAllocator(.{}){};
pub fn main() !void {
    defer _ = gpa.deinit();
    // observer position
    const obs_r = 1000.0;
    const obs_theta = std.math.degreesToRadians(f64, 85.0);
    const bl1 = bl.BoyerLindquist(f64).initImpactParameters(1.0, 0.998, 7.0, 0.0, obs_r, obs_theta);

    std.debug.print("L:{e} Q:{e}\n", .{ bl1.angular_momentum, bl1.carters_q });

    var prob = Integrator.init(problem, bl1);
    const u: [4]f64 = .{ 0.0, obs_r, obs_theta, 0.0 };

    var solver = prob.solver(gpa.allocator());

    var sol = try solver.solve(u, 0.0, 2000.0, .{ .save = true, .dt = 0.06, .max_iters = 30_000, .callback = callback });
    defer sol.deinit();

    const stdout = std.io.getStdIn();
    try sol.printInfo(stdout);

    var file = try std.fs.cwd().openFile("out.txt", .{ .mode = std.fs.File.OpenMode.write_only });
    defer file.close();

    for (sol.u[1..sol.index]) |*v| {
        try file.writer().print("{e}\n", .{v.*});
    }
}
