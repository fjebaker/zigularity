# Zigularity

_Itsy Bitsy Teeny Weeny Black Hole Visualisation Machiney_

A little long weekend project: relativistic ray tracing in Zig, visualizing thin accretion discs around a Kerr black hole with a handful of configurable parameters:

- observer position
- disc inner and outer radius
- black hole mass and spin

```bash
git clone https://github.com/fjebaker/zigularity
cd zigularity
zigmod fetch
zig build -Drelease-fast
```

I write relativistic ray tracing codes as part of my research, feel free to [take a look](https://github.com/astro-group-bristol/Gradus.jl) c: