using Documenter
using TrackingLoopFilters

makedocs(
    sitename = "TrackingLoopFilters.jl",
    modules = [TrackingLoopFilters],
    pages = [
        "Home" => "index.md",
    ],
)

deploydocs(
    repo = "github.com/JuliaGNSS/TrackingLoopFilters.jl.git",
)
