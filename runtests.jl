using Test
using Plots

include("BezierCurves/BezierCurves.jl")

@testset "BezierCurve2Pts" begin
    # 简单的测试
    control_points1 = [0.0 1.0 2.0; 0.0 1.0 0.0]
    curve_points1 = BezierCurves.bezier_curve(control_points1, 100)
    @test size(curve_points1) == (2, 100)
    plt1 = BezierCurves.plot_bezier_curve(control_points1)
    @test plt1 isa Plots.Plot

    # 更复杂的测试，使用 5 个控制点
    control_points2 = [0.0 1.0 2.0 3.0 4.0; 0.0 2.0 -1.0 1.5 0.0]
    curve_points2 = BezierCurves.bezier_curve(control_points2, 200)
    @test size(curve_points2) == (2, 200)
    plt2 = BezierCurves.plot_bezier_curve(control_points2)
    @test plt2 isa Plots.Plot
    
    plt = plot(plt1, plt2, layout=(1, 2), size=(800, 400))

    display(plt)
end

