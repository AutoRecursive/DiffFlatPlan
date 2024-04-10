module BezierCurves

include("BezierCurve.jl")
include("visualization.jl")

using .BezierCurve
using .Visualization

export bezier_curve, evaluate_bezier, plot_bezier_curve, bezier_value, evaluate_bezier_derivative

end # module