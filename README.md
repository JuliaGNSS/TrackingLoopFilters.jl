[![Build Status](https://travis-ci.org/JuliaGNSS/TrackingLoopFilters.jl.svg?branch=master)](https://travis-ci.org/JuliaGNSS/TrackingLoopFilters.jl)
[![Coverage Status](https://coveralls.io/repos/github/JuliaGNSS/TrackingLoopFilters.jl/badge.svg?branch=master)](https://coveralls.io/github/JuliaGNSS/TrackingLoopFilters.jl?branch=master)

# TrackingLoopFilters
This implements multiple loop filters for GNSS tracking algorithms.

## Features

* First, second and third order loop filters
* Boxcar and bilinear loop filters

## Getting started

Install:
```julia
julia> ]
pkg> add TrackingLoopFilters
```

## Usage

```julia
using TrackingLoopFilters
carrier_loop_filter = ThirdOrderBilinearLF()
output, next_carrier_loop_filter = filter_loop(carrier_loop_filter, discriminator_output, Δt, bandwidth)
```
