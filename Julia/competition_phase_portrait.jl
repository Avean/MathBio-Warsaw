using Plots

# -----------------------------------------------------------------------------
# Parameters that can be changed for the slide
# -----------------------------------------------------------------------------

const r_u = 1.00
const r_v = 0.85
const alpha = 0.65
const beta = 0.75

const u0 = 1.35
const v0 = 1.38

const u_max = 1.55
const v_max = 1.55
const t_max = 14.0
const dt = 0.01

# Every nonzero vector in the phase portrait has this length
const arrow_length = 0.085

const output_path = normpath(joinpath(
    @__DIR__,
    "..",
    "Latex",
    "figures",
    "local_dynamics",
    "vector_field.png",
))


# -----------------------------------------------------------------------------
# Classical nondimensional two-species competition model
#
#   u_t = r_u u (1 - u - alpha v)
#   v_t = r_v v (1 - v - beta u)
# -----------------------------------------------------------------------------

function competition_rhs(u, v)
    du = r_u * u * (1.0 - u - alpha * v)
    dv = r_v * v * (1.0 - v - beta * u)
    return du, dv
end


# -----------------------------------------------------------------------------
# Fixed-step fourth-order Runge--Kutta method
# -----------------------------------------------------------------------------

function rk4_trajectory(u_initial, v_initial, final_time, time_step)
    number_of_steps = Int(ceil(final_time / time_step))
    u = Vector{Float64}(undef, number_of_steps + 1)
    v = Vector{Float64}(undef, number_of_steps + 1)

    u[1] = u_initial
    v[1] = v_initial

    for n in 1:number_of_steps
        un = u[n]
        vn = v[n]

        k1u, k1v = competition_rhs(un, vn)
        k2u, k2v = competition_rhs(
            un + 0.5 * time_step * k1u,
            vn + 0.5 * time_step * k1v,
        )
        k3u, k3v = competition_rhs(
            un + 0.5 * time_step * k2u,
            vn + 0.5 * time_step * k2v,
        )
        k4u, k4v = competition_rhs(
            un + time_step * k3u,
            vn + time_step * k3v,
        )

        u[n + 1] = un + time_step * (k1u + 2k2u + 2k3u + k4u) / 6
        v[n + 1] = vn + time_step * (k1v + 2k2v + 2k3v + k4v) / 6
    end

    return u, v
end


# -----------------------------------------------------------------------------
# Normalized vector field: every visible arrow has the same length
# -----------------------------------------------------------------------------

u_grid = range(0.07, u_max - 0.07, length=17)
v_grid = range(0.07, v_max - 0.07, length=17)

arrow_u = Float64[]
arrow_v = Float64[]
arrow_du = Float64[]
arrow_dv = Float64[]

for v in v_grid, u in u_grid
    du, dv = competition_rhs(u, v)
    speed = hypot(du, dv)

    if speed > 1.0e-10
        push!(arrow_u, u)
        push!(arrow_v, v)
        push!(arrow_du, arrow_length * du / speed)
        push!(arrow_dv, arrow_length * dv / speed)
    end
end


# -----------------------------------------------------------------------------
# Trajectory and figure
# -----------------------------------------------------------------------------

trajectory_u, trajectory_v = rk4_trajectory(u0, v0, t_max, dt)

hydra_navy = colorant"#020817"
hydra_line = colorant"#263B59"
hydra_white = colorant"#F6F8FC"
hydra_muted = colorant"#AAB7CD"
hydra_cyan = colorant"#43D5FF"
hydra_gold = colorant"#FFBE4A"

phase_portrait = plot(
    xlims=(0.0, u_max),
    ylims=(0.0, v_max),
    aspect_ratio=:equal,
    size=(1000, 760),
    dpi=180,
    background_color=hydra_navy,
    background_color_inside=hydra_navy,
    foreground_color_subplot=hydra_muted,
    foreground_color_axis=hydra_muted,
    foreground_color_border=hydra_muted,
    foreground_color_text=hydra_muted,
    grid=false,
    legend=false,
    framestyle=:axes,
    xlabel="u",
    ylabel="v",
    guidefontcolor=hydra_white,
    guidefontsize=16,
    tickfontcolor=hydra_muted,
    tickfontsize=11,
    margin=5Plots.mm,
)

quiver!(
    phase_portrait,
    arrow_u,
    arrow_v,
    quiver=(arrow_du, arrow_dv),
    color=hydra_line,
    linewidth=1.15,
    arrow=:closed,
)

plot!(
    phase_portrait,
    trajectory_u,
    trajectory_v,
    color=hydra_gold,
    linewidth=3.4,
)

scatter!(
    phase_portrait,
    [u0],
    [v0],
    color=hydra_cyan,
    markerstrokecolor=hydra_white,
    markerstrokewidth=1.2,
    markersize=7.5,
)

annotate!(
    phase_portrait,
    u0 - 0.03,
    v0 + 0.09,
    text("(u₀, v₀)", 12, hydra_white, :right),
)

mkpath(dirname(output_path))
savefig(phase_portrait, output_path)

println("Saved phase portrait to: " * output_path)
