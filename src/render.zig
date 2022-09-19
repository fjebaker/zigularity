const std = @import("std");
const zigode = @import("zigode");
const bl = @import("./boyer-lindquist.zig");

// type aliases
pub const BoyerLindquist = bl.BoyerLindquist(f64);
pub const Solver = zigode.Solver(f64, 4, BoyerLindquist);
pub const Integrator = zigode.AdaptiveTsit5(f64, 4, BoyerLindquist);

pub const ImageConfig = struct {
    observer_distance: f64 = 1000.0,
    observer_inclination: f64 = 85.0, // in degrees
    image_width: usize = 400,
    image_height: usize = 200,
    dtmax: f64 = 3.0,
    mass: f64 = 1.0,
    spin: f64 = 0.998,
    // impact parameters
    min_alpha: f64 = -60.0,
    max_alpha: f64 = 60.0,
    min_beta: f64 = -30.0,
    max_beta: f64 = 30.0,
};

pub fn run(
    allocator: std.mem.Allocator,
    comptime problem: fn (du: *[4]f64, u: *const [4]f64, t: f64, p: *BoyerLindquist) void,
    comptime callback: fn (ptr: *Solver, u: *const [4]f64, t: f64, p: *BoyerLindquist) void,
    config: ImageConfig,
    outfile: []const u8,
) !void {
    // observer position
    const obs_r = config.observer_distance;
    const obs_theta = std.math.degreesToRadians(f64, config.observer_inclination);
    // image properties
    const height = config.image_height;
    const width = config.image_width;

    // black hole parameters
    const mass = config.mass;
    const spin = config.spin;

    // initial position four vector
    const u: [4]f64 = .{ 0.0, obs_r, obs_theta, 0.0 };

    // define impact parameter ranges
    const min_alpha = config.min_alpha;
    const min_beta = config.min_beta;
    const max_alpha = config.max_alpha;
    const max_beta = config.max_beta;

    const alpha_step = (max_alpha - min_alpha) / @intToFloat(f64, width);
    // add slight offset to avoid central coordinate singularity
    const beta_step = (max_beta - min_beta + 1e-2) / @intToFloat(f64, height);

    // allocate output image
    var image: [][]f64 = try allocator.alloc([]f64, width);
    defer allocator.free(image);
    for (image) |*row| {
        row.* = try allocator.alloc(f64, height);
    }
    defer for (image) |row| {
        allocator.free(row);
    };

    var beta: f64 = min_beta;
    var i: usize = 0;

    std.debug.print("Starting render:\n", .{});
    while (beta <= max_beta) : (beta += beta_step) {
        // print which row we just finished
        std.debug.print("{d:>5}/{d}\n", .{ i + 1, height });

        var alpha: f64 = min_alpha;
        var j: usize = 0;
        while (alpha <= max_alpha) : (alpha += alpha_step) {
            // create parameters for this photon
            const p = BoyerLindquist.initImpactParameters(mass, spin, alpha, beta, obs_r, obs_theta);
            // initialise the ODE problem
            var prob = Integrator.init(problem, p);
            // set dtmax to something small so we don't miss the disc
            prob.dtmax = config.dtmax;
            var solver = prob.solver(allocator);
            var sol = try solver.solve(
                u,
                0.0,
                2000.0,
                .{ .save = false, .dt = 1.01, .max_iters = 30_000, .callback = callback },
            );
            defer sol.deinit();

            // read out data
            image[j][i] = sol.u[0][0];

            j += 1;
            // safety catch
            if (j >= width) break;
        }
        i += 1;
        // safety catch
        if (i >= height) break;
    }
    std.debug.print("Done: writing to file {s}\n", .{outfile});
    // save to file
    var file = try std.fs.cwd().createFile(outfile, .{});
    defer file.close();

    var writer = file.writer();
    for (image) |row| {
        for (row) |v| {
            try writer.print("{e},", .{v});
        }
        try writer.print("\n", .{});
    }
}
