module BezierCurve

export bezier_curve, evaluate_bezier, bezier_value, evaluate_bezier_derivative

# 使用矩阵形式计算Bezier曲线上的点
function evaluate_bezier(t::Real, control_points::AbstractMatrix{T}) where {T<:Real}
    n = size(control_points, 2) - 1
    u = 1 - t
    bc = [binomial(n, i) for i in 0:n]
    terms = [bc[i+1] * t^i * u^(n-i) for i in 0:n]
    return control_points * terms
end

# 生成Bezier曲线上的点
function bezier_curve(control_points::AbstractMatrix{T}, num_points::Integer) where {T<:Real}
    t_values = range(0, 1, length=num_points)
    return hcat([evaluate_bezier(t, control_points) for t in t_values]...)
end

# Helper functions
function bezier_value(t, control_points)
    n = size(control_points, 2) - 1
    sum(binomial(n, i) * (1 - t)^(n - i) * t^i * control_points[:, i+1] for i in 0:n)
end

# Evaluate the derivative of the Bezier curve at parameter t
function evaluate_bezier_derivative(t, control_points)
    n = size(control_points, 2) - 1
    sum(
        binomial(n - 1, i) * (1 - t)^(n - 1 - i) * t^i * (control_points[:, i+2] - control_points[:, i+1])
        for i in 0:(n-1)
    ) * n
end
end # module