1-D Hydrodyanmic model
---

This is a channel model, which implements depth-averaged momentum and continuity equations.
The channel is represented as an equispaced grid (dx). The surface evelation is forced at either
end of the channel, using the etaspec function (see config). Hydrodynamic models are time-stepping,
but the delta time (dt in this case) must not be less than the time it take a parcel of water moving
at the maximum velocity to traverse a cell dimension (dx). To select an appropriate dt, consider
the [courant number](https://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition)

__Copyright (C) 2020, SoftWaring Solutions ATF The Miss Trust__
