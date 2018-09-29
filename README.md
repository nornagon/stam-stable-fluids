# Real-Time Fluid Dynamics for Games, Jos Stam 2003 [JS port]

This is a direct port from the C code presented in Stam's 2003 paper, with a small convenience wrapper.

Example:
```js
const {Fluid} = require('stam-stable-fluids')
const size = 128
const f = new Fluid(size)
for (let y = 0; y < size; y++)
  for (let x = 0; x < size/2; x++)
    f.addDensity(x, y, 10)
for (let y = 0; y < size; y++)
  f.addForce(10, y, 10, 0)
f.step(1/60)

const [vx, vy] = f.velocityAt(0, 0)
const d = f.densityAt(0, 0)
```

## API

#### `new Fluid(size, viscosity = 0, diffusion = 0)`

Make a new fluid simulation with width and height equal to `size`. If you want a goopier fluid, like honey, set `viscosity` to a small positive number, like `0.0001` (larger numbers will probably make the fluid so "thick" it will barely move). `diffusion` will make the "stuff" in the fluid spread out over time. The default for both viscosity and diffusion, if not provided, is 0.

#### `Fluid.step(dt)`

Step the simulation by `dt` (e.g. `1/60`).

#### `Fluid.velocityAt(x, y)`

Returns `[vx, vy]` representing the velocity of the fluid at (_x_, _y_).

#### `Fluid.densityAt(x, y)`

Returns the density of "stuff" at (_x_, _y_).

#### `Fluid.addForce(x, y, fx, fy)`

Adds an external force to the fluid, e.g. from a UI or other system. The force field will be cleared after every `step()`.

#### `Fluid.addDensity(x, y, d)`

Adds an external density flow to the fluid. `d` is the rate of flow. The density flow field will be cleared after every `step()`.
