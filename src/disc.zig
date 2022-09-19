const std = @import("std");
const zigode = @import("zigode");
const bl = @import("./boyer-lindquist.zig");

// type aliases
const BoyerLindquist = bl.BoyerLindquist(f64);
const Solver = zigode.Solver(f64, 4, BoyerLindquist);
const Integrator = zigode.AdaptiveTsit5(f64, 4, BoyerLindquist);

const disc_inner_radius = 10;
const disc_outer_radius = 50;

fn callback(solver: *Solver, u: *const [4]f64, _: f64, p: *BoyerLindquist) void {
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

fn problem(du: *[4]f64, u: *const [4]f64, _: f64, p: *BoyerLindquist) void {
    const v: [4]f64 = p.fourVelocity(u);
    for (du) |*a, i| {
        a.* = v[i];
    }
}

pub fn run(allocator: std.mem.Allocator) !void {
    // observer position
    const obs_r = 1000.0;
    const obs_theta = std.math.degreesToRadians(f64, 85.0);
    // image properties
    const height = 200;
    const width = 400;

    // black hole parameters
    const mass = 1.0;
    const spin = 0.998;

    // allocate output image
    var image: [][]f64 = try allocator.alloc([]f64, width);
    defer allocator.free(image);
    for (image) |*row| {
        row.* = try allocator.alloc(f64, height);
    }
    defer for (image) |row| {
        allocator.free(row);
    };

    // initial position four vector
    const u: [4]f64 = .{ 0.0, obs_r, obs_theta, 0.0 };

    // define impact parameter ranges
    const square = 30.0;
    const min_alpha = -2 * square;
    const min_beta = -square;
    const max_alpha = 2 * square;
    const max_beta = square;

    const alpha_step = (max_alpha - min_alpha) / @as(f64, width);
    // add slight offset to avoid central coordinate singularity
    const beta_step = (max_beta - min_beta + 1e-2) / @as(f64, height);

    var alpha: f64 = min_alpha;
    var i: usize = 0;
    while (alpha <= max_alpha) : (alpha += alpha_step) {
        var j: usize = 0;
        std.debug.print("{d}\n", .{i});

        var beta: f64 = min_beta;
        while (beta <= max_beta) : (beta += beta_step) {
            // std.debug.print("alpha {e}, beta {e}\n", .{alpha, beta});
            const p = BoyerLindquist.initImpactParameters(mass, spin, alpha, beta, obs_r, obs_theta);
            var prob = Integrator.init(problem, p);
            prob.dtmax = 3.0;
            var solver = prob.solver(allocator);
            var sol = try solver.solve(
                u,
                0.0,
                2000.0,
                .{ .save = false, .dt = 1.01, .max_iters = 30_000, .callback = callback },
            );
            defer sol.deinit();

            // read out data
            image[i][j] = sol.u[0][0];

            j += 1;
            if (j >= height) break;
        }
        i += 1;
        if (i >= width) break;
    }

    var file = try std.fs.cwd().createFile("disc.data", .{});
    defer file.close();

    var writer = file.writer();
    for (image) |row| {
        for (row) |v| {
            try writer.print("{e},", .{v});
        }
        try writer.print("\n", .{});
    }
}
