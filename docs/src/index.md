# TrackingLoopFilters.jl

Loop filters for GNSS tracking algorithms.

## Installation

```julia
pkg> add TrackingLoopFilters
```

## Quick Start

```julia
using TrackingLoopFilters
using Unitful: Hz, s

# Create a loop filter
lf = ThirdOrderBilinearLF()

# Process discriminator output
output, next_lf = filter_loop(lf, δθ, Δt, bandwidth)
```

## API Reference

```@autodocs
Modules = [TrackingLoopFilters]
```
