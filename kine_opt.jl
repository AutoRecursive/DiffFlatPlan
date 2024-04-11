using JuMP
using Ipopt
using Plots

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

function bicycle_model(x, y, v, ψ, a, δ, dt)
	dx = v * cos(ψ)
	dy = v * sin(ψ)
	dv = a
	dψ = v * tan(δ) / L

	x_next = x + dx * dt
	y_next = y + dy * dt
	v_next = v + dv * dt
	ψ_next = ψ + dψ * dt

	return x_next, y_next, v_next, ψ_next
end

function optimize_trajectory(num_control_points)
	model = Model(Ipopt.Optimizer)

	# Decision variables
	@variable(model, ac[1:num_control_points])
	@variable(model, δc[1:num_control_points])

	# 初始化轨迹变量
	num_samples = Int(total_time / dt) + 1
	@variable(model, x_traj[1:num_samples])
	@variable(model, y_traj[1:num_samples])
	@variable(model, v_traj[1:num_samples])
	@variable(model, ψ_traj[1:num_samples])

	# 设置初始条件
	@constraint(model, x_traj[1] == 0.0)
	@constraint(model, y_traj[1] == 0.0)
	@constraint(model, v_traj[1] == v_initial)
	@constraint(model, ψ_traj[1] == 0.0)

	# 使用自行车模型更新轨迹
	for i in 2:num_samples
		@NLconstraint(model, x_traj[i] == bicycle_model(x_traj[i-1], y_traj[i-1], v_traj[i-1], ψ_traj[i-1],
			bezier_value((i - 1) / (num_samples - 1), [ac'; δc'])[1],
			bezier_value((i - 1) / (num_samples - 1), [ac'; δc'])[2],
			dt)[1])
		@NLconstraint(model, y_traj[i] == bicycle_model(x_traj[i-1], y_traj[i-1], v_traj[i-1], ψ_traj[i-1],
			bezier_value((i - 1) / (num_samples - 1), [ac'; δc'])[1],
			bezier_value((i - 1) / (num_samples - 1), [ac'; δc'])[2],
			dt)[2])
		@NLconstraint(model, v_traj[i] == bicycle_model(x_traj[i-1], y_traj[i-1], v_traj[i-1], ψ_traj[i-1],
			bezier_value((i - 1) / (num_samples - 1), [ac'; δc'])[1],
			bezier_value((i - 1) / (num_samples - 1), [ac'; δc'])[2],
			dt)[3])
		@NLconstraint(model, ψ_traj[i] == bicycle_model(x_traj[i-1], y_traj[i-1], v_traj[i-1], ψ_traj[i-1],
			bezier_value((i - 1) / (num_samples - 1), [ac'; δc'])[1],
			bezier_value((i - 1) / (num_samples - 1), [ac'; δc'])[2],
			dt)[4])
	end

	# 目标函数
	@NLobjective(model, Min,
		sum(
            ifelse((i - 1) * dt > 3, (y_traj[i] - y_desired)^2, 0) +
			(v_traj[i] - v_desired)^2 +
            100 * (ψ_traj[i] - 0)^2 +
			1 * (bezier_value(i / num_samples, [ac'; δc'])[1])^2 +
			1 * (bezier_value(i / num_samples, [ac'; δc'])[2])^2
			for i in 1:num_samples
		)
	)

	# 收紧约束条件
	for i in 1:num_samples
		@NLconstraint(model, bezier_value(i / num_samples, [ac'; δc'])[1] <= 1.0)
		@NLconstraint(model, bezier_value(i / num_samples, [ac'; δc'])[1] >= -3.0)
		@NLconstraint(model, bezier_value(i / num_samples, [ac'; δc'])[2] <= π / 3)
		@NLconstraint(model, bezier_value(i / num_samples, [ac'; δc'])[2] >= -π / 3)
	end

	# 提供初始猜测值
	set_start_value.(ac, 0.1 * ones(num_control_points))
	set_start_value.(δc, 0.1 * ones(num_control_points))

	# 优化
	optimize!(model)

	# 返回优化后的控制点和轨迹
	return value.(ac), value.(δc), value.(x_traj), value.(y_traj), value.(v_traj), value.(ψ_traj)
end
# Run optimization
num_control_points = 5
a_opt, δ_opt, x_traj_opt, y_traj_opt, v_traj_opt, ψ_traj_opt = optimize_trajectory(num_control_points)
print(a_opt, δ_opt)

# Generate points for the optimized trajectory
num_samples = Int(total_time / dt) + 1
t = range(0, total_time, length = num_samples)
x_traj = zeros(num_samples)
y_traj = zeros(num_samples)
v_traj = zeros(num_samples)
ψ_traj = zeros(num_samples)

x_traj[1] = 0.0
y_traj[1] = 0.0
v_traj[1] = v_initial
ψ_traj[1] = 0.0

for i in 2:num_samples
	a = bezier_value((i - 1) / (num_samples - 1), [a_opt'; δ_opt'])[1]
	δ = bezier_value((i - 1) / (num_samples - 1), [a_opt'; δ_opt'])[2]
	x_traj[i], y_traj[i], v_traj[i], ψ_traj[i] = bicycle_model(x_traj[i-1], y_traj[i-1], v_traj[i-1], ψ_traj[i-1], a, δ, dt)
end

# Plot x(t)
plt_xt = plot(t, x_traj, xlabel = "t (s)", ylabel = "x (m)", title = "x(t)", legend = false)

# Plot y(t)  
plt_yt = plot(t, y_traj, xlabel = "t (s)", ylabel = "y (m)", title = "y(t)", legend = false)

# Plot v(t)
plt_vt = plot(t, v_traj, xlabel = "t (s)", ylabel = "v (m/s)", title = "v(t)", legend = false)

# Plot ψ(t)
plt_ψt = plot(t, ψ_traj .* 180 / π, xlabel = "t (s)", ylabel = "ψ (deg)", title = "ψ(t)", legend = false)

# Plot y(x)
plt_yx = plot(x_traj, y_traj, xlabel = "x (m)", ylabel = "y (m)", title = "y(x)", aspect_ratio = :equal, legend = false)

# Combine the plots into a 2x3 layout  
plt = plot(plt_xt, plt_yt, plt_vt, plt_ψt, plt_yx, layout = (2, 3), size = (1200, 600))

# Display the combined plot
display(plt)