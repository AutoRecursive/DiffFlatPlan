module Visualization

using Plots
using ..BezierCurve: bezier_curve

export plot_bezier_curve

function plot_bezier_curve(control_points::AbstractMatrix{T}, num_points::Integer=100) where {T<:Real}
    curve_points = bezier_curve(control_points, num_points)
    plt = plot(curve_points[1,:], curve_points[2,:], label="Bezier Curve")
    scatter!(control_points[1,:], control_points[2,:], label="Control Points")
    return plt
end

end # module