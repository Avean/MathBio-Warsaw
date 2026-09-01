using Plots

# -----------------------------------------------------------------------------
# Output files required by the Lyapunov-stability slide
# -----------------------------------------------------------------------------

const output_directory = normpath(joinpath(
    @__DIR__,
    "..",
    "Latex",
    "figures",
    "B_reaction_diffusion",
))

const stable_output = joinpath(output_directory, "stability_converging.png")
const unstable_output = joinpath(output_directory, "stability_diverging.png")


# -----------------------------------------------------------------------------
# Colours used by the LaTeX presentation
# -----------------------------------------------------------------------------

const hydra_navy = colorant"#020817"
const hydra_line = colorant"#263B59"
const hydra_white = colorant"#F6F8FC"
const hydra_muted = colorant"#AAB7CD"
const hydra_cyan = colorant"#43D5FF"
const hydra_gold = colorant"#FFBE4A"


# -----------------------------------------------------------------------------
# Generic fourth-order Runge--Kutta solver for a planar autonomous system
# -----------------------------------------------------------------------------

function rk4_trajectory(rhs, initial_state, final_time; dt=0.01)
    number_of_steps = Int(ceil(final_time / dt))
    u = Vector{Float64}(undef, number_of_steps + 1)
    v = Vector{Float64}(undef, number_of_steps + 1)

    u[1], v[1] = initial_state

    for n in 1:number_of_steps
        un = u[n]
        vn = v[n]

        k1u, k1v = rhs(un, vn)
        k2u, k2v = rhs(un + 0.5dt * k1u, vn + 0.5dt * k1v)
        k3u, k3v = rhs(un + 0.5dt * k2u, vn + 0.5dt * k2v)
        k4u, k4v = rhs(un + dt * k3u, vn + dt * k3v)

        u[n + 1] = un + dt * (k1u + 2k2u + 2k3u + k4u) / 6
        v[n + 1] = vn + dt * (k1v + 2k2v + 2k3v + k4v) / 6
    end

    return u, v
end


# -----------------------------------------------------------------------------
# Equal-length arrows for a planar vector field
# -----------------------------------------------------------------------------

function normalized_vector_field(rhs, u_limits, v_limits;
                                 grid_size=16, arrow_length=0.075)
    u_margin = 0.045 * (u_limits[2] - u_limits[1])
    v_margin = 0.045 * (v_limits[2] - v_limits[1])

    u_grid = range(
        u_limits[1] + u_margin,
        u_limits[2] - u_margin,
        length=grid_size,
    )
    v_grid = range(
        v_limits[1] + v_margin,
        v_limits[2] - v_margin,
        length=grid_size,
    )

    arrow_u = Float64[]
    arrow_v = Float64[]
    arrow_du = Float64[]
    arrow_dv = Float64[]

    for v in v_grid, u in u_grid
        du, dv = rhs(u, v)
        speed = hypot(du, dv)

        if speed > 1.0e-10
            push!(arrow_u, u)
            push!(arrow_v, v)
            push!(arrow_du, arrow_length * du / speed)
            push!(arrow_dv, arrow_length * dv / speed)
        end
    end

    @assert all(
        isapprox(hypot(du, dv), arrow_length; atol=1.0e-12)
        for (du, dv) in zip(arrow_du, arrow_dv)
    )

    return arrow_u, arrow_v, arrow_du, arrow_dv
end


function base_phase_portrait(u_limits, v_limits)
    return plot(
        xlims=u_limits,
        ylims=v_limits,
        aspect_ratio=:equal,
        size=(900, 720),
        dpi=180,
        background_color=hydra_navy,
        background_color_inside=hydra_navy,
        foreground_color_subplot=hydra_muted,
        foreground_color_axis=hydra_muted,
        foreground_color_border=hydra_muted,
        foreground_color_text=hydra_muted,
        guidefontcolor=hydra_white,
        guidefontsize=16,
        tickfontcolor=hydra_muted,
        tickfontsize=10,
        xlabel="u",
        ylabel="v",
        grid=false,
        legend=false,
        framestyle=:axes,
        margin=5Plots.mm,
    )
end


function add_vector_field!(phase_portrait, rhs, u_limits, v_limits;
                           arrow_length=0.075)
    arrow_u, arrow_v, arrow_du, arrow_dv = normalized_vector_field(
        rhs,
        u_limits,
        v_limits;
        arrow_length=arrow_length,
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
end


function add_two_trajectories!(phase_portrait, rhs, initial_a, initial_b,
                               final_time; dt=0.01,
                               label_offsets=((0.0, 0.0), (0.0, 0.0)))
    u_a, v_a = rk4_trajectory(rhs, initial_a, final_time; dt=dt)
    u_b, v_b = rk4_trajectory(rhs, initial_b, final_time; dt=dt)

    plot!(phase_portrait, u_a, v_a, color=hydra_cyan, linewidth=3.2)
    plot!(phase_portrait, u_b, v_b, color=hydra_gold, linewidth=3.2)

    scatter!(
        phase_portrait,
        [initial_a[1]],
        [initial_a[2]],
        color=hydra_cyan,
        markerstrokecolor=hydra_white,
        markerstrokewidth=1.5,
        markersize=9.5,
    )
    scatter!(
        phase_portrait,
        [initial_b[1]],
        [initial_b[2]],
        color=hydra_gold,
        markerstrokecolor=hydra_white,
        markerstrokewidth=1.5,
        markersize=9.5,
        markershape=:diamond,
    )

    annotate!(
        phase_portrait,
        initial_a[1] + label_offsets[1][1],
        initial_a[2] + label_offsets[1][2],
        text("(u₀, v₀)", 26, hydra_cyan),
    )
    annotate!(
        phase_portrait,
        initial_b[1] + label_offsets[2][1],
        initial_b[2] + label_offsets[2][2],
        text("(ũ₀, ṽ₀)", 26, hydra_gold),
    )

    initial_distance = hypot(
        initial_a[1] - initial_b[1],
        initial_a[2] - initial_b[2],
    )
    final_distance = hypot(u_a[end] - u_b[end], v_a[end] - v_b[end])

    return initial_distance, final_distance
end


# -----------------------------------------------------------------------------
# Example 1: stable coexistence in the competition model
# -----------------------------------------------------------------------------

const stable_r_u = 1.00
const stable_r_v = 0.85
const stable_alpha = 0.65
const stable_beta = 0.75

stable_rhs(u, v) = (
    stable_r_u * u * (1.0 - u - stable_alpha * v),
    stable_r_v * v * (1.0 - v - stable_beta * u),
)

stable_initial_a = (1.30, 0.20)
stable_initial_b = (1.26, 0.24)
stable_u_limits = (0.0, 1.55)
stable_v_limits = (0.0, 1.55)

stable_portrait = base_phase_portrait(stable_u_limits, stable_v_limits)
add_vector_field!(
    stable_portrait,
    stable_rhs,
    stable_u_limits,
    stable_v_limits;
    arrow_length=0.078,
)

stable_initial_distance, stable_final_distance = add_two_trajectories!(
    stable_portrait,
    stable_rhs,
    stable_initial_a,
    stable_initial_b,
    50.0;
    label_offsets=((0.08, -0.05), (-0.08, 0.06)),
)

@assert stable_final_distance < 0.02 * stable_initial_distance


# -----------------------------------------------------------------------------
# Example 2: nearby states separate across an unstable threshold
#
#   u_t = 2.2 u (1-u) (u-1/2)
#   v_t = -0.65 (v-1/4)
#
# The equilibrium (1/2,1/4) is unstable in the u direction
# -----------------------------------------------------------------------------

unstable_rhs(u, v) = (
    2.2 * u * (1.0 - u) * (u - 0.5),
    -0.65 * (v - 0.25),
)

unstable_initial_a = (0.48, 1.05)
unstable_initial_b = (0.52, 1.05)
unstable_u_limits = (-0.03, 1.03)
unstable_v_limits = (0.15, 1.21)

unstable_portrait = base_phase_portrait(unstable_u_limits, unstable_v_limits)
add_vector_field!(
    unstable_portrait,
    unstable_rhs,
    unstable_u_limits,
    unstable_v_limits;
    arrow_length=0.055,
)

unstable_initial_distance, unstable_final_distance = add_two_trajectories!(
    unstable_portrait,
    unstable_rhs,
    unstable_initial_a,
    unstable_initial_b,
    12.0;
    label_offsets=((-0.075, -0.05), (0.075, 0.05)),
)

@assert unstable_final_distance > 10.0 * unstable_initial_distance


# -----------------------------------------------------------------------------
# Save both figures
# -----------------------------------------------------------------------------

mkpath(output_directory)
savefig(stable_portrait, stable_output)
savefig(unstable_portrait, unstable_output)

println("Stable example")
println("  initial distance = ", stable_initial_distance)
println("  final distance   = ", stable_final_distance)
println("  saved to         = ", stable_output)

println("Unstable example")
println("  initial distance = ", unstable_initial_distance)
println("  final distance   = ", unstable_final_distance)
println("  saved to         = ", unstable_output)
