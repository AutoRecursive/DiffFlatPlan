using JuMP
using Ipopt
using Plots
include("BezierCurves/BezierCurves.jl")

# Bicycle model parameters
const L = 2.0  # wheelbase (m)
const v_desired = 11.1  # desired velocity (m/s)
const v_initial = 0.0  # initial velocity (m/s)
const y_desired = 3.75  # initial velocity (m/s)

const total_time = 10.0  # total time (s)
const dt = 0.1  # time interval between samples (s)
const ε = 1e-8  # Small positive constant to avoid division by zero

# Helper functions
function bezier_value(t, control_points)
	n = size(control_points, 2) - 1
	sum(binomial(n, i) * (1 - t)^(n - i) * t^i * control_points[:, i+1] for i in 0:n)
end

function bezier_derivative(t, control_points)
	n = size(control_points, 2) - 1
	sum(
		binomial(n - 1, i) * (1 - t)^(n - 1 - i) * t^i * (control_points[:, i+2] - control_points[:, i+1])
		for i in 0:(n-1)
	) * n
end

function optimize_trajectory(num_control_points)
    model = Model(Ipopt.Optimizer)

    # Decision variables
    @variable(model, xc[1:num_control_points])
    @variable(model, yc[1:num_control_points])

    # Objective function
    num_samples = Int(total_time / dt) + 1
    @NLobjective(model, Min,
        sum(
            (bezier_derivative(i / num_samples, [xc'; yc'])[1] - v_desired)^2 +
            (bezier_value(i / num_samples, [xc'; yc'])[2] - 3.75)^2 +
            (bezier_derivative(i / num_samples, [bezier_derivative(i / num_samples, [xc'; yc'])[1]; yc'])[1])^2 +
            (bezier_derivative(i / num_samples, [bezier_derivative(i / num_samples, [bezier_derivative(i / num_samples, [xc'; yc'])[1]; yc'])[1]; yc'])[1])^2 +
            (atan(clamp(bezier_derivative(i / num_samples, [xc'; yc'])[2] / (bezier_derivative(i / num_samples, [xc'; yc'])[1] + ε), -10, 10)))^2
            for i in 1:num_samples
        )
    )

    # Constraints
    @constraint(model, xc[1] == 0.0)
    @constraint(model, yc[1] == 0.0)
    # @constraint(model, yc[end] == 3.75)  # Fix the y-coordinate of the last control point to 3.75
    @constraint(model, bezier_derivative(0.0, [xc'; yc'])[1] == v_initial)

    # Constraint on initial heading
    @NLconstraint(model, atan(bezier_derivative(0.0, [xc'; yc'])[2] / (bezier_derivative(0.0, [xc'; yc'])[1] + ε)) == 0.0)

    # Constraint on y-velocity
    for i in 1:num_samples
        @NLconstraint(model, bezier_derivative(i / num_samples, [xc'; yc'])[2] <= 1.0)
        @NLconstraint(model, bezier_derivative(i / num_samples, [xc'; yc'])[2] >= -1.0)
    end

    # Optimize
    optimize!(model)

    # Return optimized control points
    return value.(xc), value.(yc)
end

# Run optimization
num_control_points = 5
x_opt, y_opt = optimize_trajectory(num_control_points)

# Generate points for the optimized trajectory
num_samples = Int(total_time / dt) + 1
t = range(0, total_time, length = num_samples)
x_traj = [bezier_value(ti / total_time, [x_opt'; y_opt'])[1] for ti in t]
y_traj = [bezier_value(ti / total_time, [x_opt'; y_opt'])[2] for ti in t]
v_traj = [bezier_derivative(ti / total_time, [x_opt'; y_opt'])[1] for ti in t]

# Plot y(x)
plt_yx = plot(x_traj, y_traj, xlabel = "x (m)", ylabel = "y (m)", title = "y(x)", legend = false)

# Plot x(t)
plt_xt = plot(t, x_traj, xlabel = "t (s)", ylabel = "x (m)", title = "x(t)", legend = false)

# Plot y(t)
plt_yt = plot(t, y_traj, xlabel = "t (s)", ylabel = "y (m)", title = "y(t)", legend = false)

# Plot v(t)
plt_vt = plot(t, v_traj, xlabel = "t (s)", ylabel = "v (m/s)", title = "v(t)", legend = false)

# Calculate heading (in degrees) at each sample point
heading_traj = [atan(bezier_derivative(ti / total_time, [x_opt'; y_opt'])[2] / bezier_derivative(ti / total_time, [x_opt'; y_opt'])[1]) * 180 / π for ti in t]

# Plot heading(t)
plt_heading = plot(t, heading_traj, xlabel = "t (s)", ylabel = "heading (deg)", title = "heading(t)", legend = false)

# Combine the plots into a 2x3 layout
plt = plot(plt_yx, plt_xt, plt_yt, plt_vt, plt_heading, layout = (2, 3), size = (1200, 600))

# Display the combined plot
display(plt)
