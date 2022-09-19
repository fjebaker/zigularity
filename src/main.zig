const std = @import("std");
const zigode = @import("zigode");
const bl = @import("./boyer-lindquist.zig");

// example problems
const shadow = @import("./shadow.zig");
const disc = @import("./disc.zig");
const trace = @import("./simple-trace.zig");

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    // try shadow.run(allocator);
    // try trace.run(allocator);
    try disc.run(allocator);
}
