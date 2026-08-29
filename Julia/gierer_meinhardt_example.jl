using LinearAlgebra
using Printf
using Plots

# -----------------------------------------------------------------------------
# Concrete parameter set from Models/GiererMainhardt/GiererMainhardt.jl
# -----------------------------------------------------------------------------

const D_u = 1.0e-2
const D_v = 1.0
const a = 1.5
const b = 2.0
const mu_u = 0.5
const mu_v = 1.0
const K = 1.0

# The application uses [0, 1] with homogeneous Neumann boundary conditions
const domain_length = 1.0
const large_domain_length = 10.0

const output_directory = normpath(joinpath(
    @__DIR__,
    "..",
    "Latex",
    "figures",
    "gierer_meinhardt",
))

const nullclines_output = joinpath(
    output_directory,
    "gm_example_nullclines.png",
)

const phase_portrait_output = joinpath(
    output_directory,
    "gm_example_phase_portrait.png",
)

const mode_growth_output = joinpath(
    output_directory,
    "gm_example_mode_growth_rates.png",
)

const large_domain_mode_growth_output = joinpath(
    output_directory,
    "gm_example_mode_growth_rates_L10.png",
)


# -----------------------------------------------------------------------------
# Presentation colours
# -----------------------------------------------------------------------------

const hydra_navy = colorant"#020817"
const hydra_panel = colorant"#0B1930"
const hydra_line = colorant"#263B59"
const hydra_white = colorant"#F6F8FC"
const hydra_muted = colorant"#AAB7CD"
const hydra_cyan = colorant"#43D5FF"
const hydra_gold = colorant"#FFBE4A"
const hydra_green = colorant"#42D6A4"
const hydra_red = colorant"#FF6B6B"


# -----------------------------------------------------------------------------
# Local Gierer--Meinhardt dynamics
#
#   u_t = a u^2 / (K + v) - mu_u u
#   v_t = b u^2 - mu_v v
# -----------------------------------------------------------------------------

function gm_rhs(u, v)
    du = a * u^2 / (K + v) - mu_u * u
    dv = b * u^2 - mu_v * v
    return du, dv
end

function gm_jacobian(u, v)
    f_u = 2a * u / (K + v) - mu_u
    f_v = -a * u^2 / (K + v)^2
    g_u = 2b * u
    g_v = -mu_v

    return [f_u f_v; g_u g_v]
end


# -----------------------------------------------------------------------------
# Steady states and their local eigenvalues
# -----------------------------------------------------------------------------

const steady_state_discriminant =
    (a / mu_u)^2 - 4b * K / mu_v

@assert steady_state_discriminant > 0

const u_minus = mu_v / (2b) * (
    a / mu_u - sqrt(steady_state_discriminant)
)
const u_plus = mu_v / (2b) * (
    a / mu_u + sqrt(steady_state_discriminant)
)
const v_minus = b / mu_v * u_minus^2
const v_plus = b / mu_v * u_plus^2

const steady_states = (
    (name = "zero steady state", point = (0.0, 0.0)),
    (name = "lower positive steady state", point = (u_minus, v_minus)),
    (name = "upper positive steady state", point = (u_plus, v_plus)),
)


# -----------------------------------------------------------------------------
# Turing band and growth rates on [0, 1] with Neumann conditions
# -----------------------------------------------------------------------------

const J_plus = gm_jacobian(u_plus, v_plus)
const determinant_J_plus = det(J_plus)
const A_GM = D_v * mu_u - D_u * mu_v
const turing_discriminant =
    A_GM^2 - 4D_u * D_v * determinant_J_plus

@assert turing_discriminant > 0

const mu_band_minus = (
    A_GM - sqrt(turing_discriminant)
) / (2D_u * D_v)
const mu_band_plus = (
    A_GM + sqrt(turing_discriminant)
) / (2D_u * D_v)

function mode_matrix(k; length=domain_length)
    mu_k = (k * pi / length)^2
    diffusion_matrix = [D_u 0.0; 0.0 D_v]
    return mu_k, J_plus - mu_k * diffusion_matrix
end


# -----------------------------------------------------------------------------
# Fourth-order Runge--Kutta trajectories
# -----------------------------------------------------------------------------

function rk4_trajectory(initial_state, final_time; dt=0.005)
    number_of_steps = Int(ceil(final_time / dt))
    u = Vector{Float64}(undef, number_of_steps + 1)
    v = Vector{Float64}(undef, number_of_steps + 1)

    u[1], v[1] = initial_state

    for n in 1:number_of_steps
        un = u[n]
        vn = v[n]

        k1u, k1v = gm_rhs(un, vn)
        k2u, k2v = gm_rhs(
            un + 0.5dt * k1u,
            vn + 0.5dt * k1v,
        )
        k3u, k3v = gm_rhs(
            un + 0.5dt * k2u,
            vn + 0.5dt * k2v,
        )
        k4u, k4v = gm_rhs(
            un + dt * k3u,
            vn + dt * k3v,
        )

        u[n + 1] = un + dt * (k1u + 2k2u + 2k3u + k4u) / 6
        v[n + 1] = vn + dt * (k1v + 2k2v + 2k3v + k4v) / 6
    end

    return u, v
end


# -----------------------------------------------------------------------------
# Equal-length arrows in normalized plot coordinates
# -----------------------------------------------------------------------------

function normalized_vector_field(u_limits, v_limits;
                                 u_count=18,
                                 v_count=18,
                                 arrow_fraction=0.030)
    u_span = u_limits[2] - u_limits[1]
    v_span = v_limits[2] - v_limits[1]

    u_grid = range(
        u_limits[1] + 0.055u_span,
        u_limits[2] - 0.055u_span,
        length=u_count,
    )
    v_grid = range(
        v_limits[1] + 0.055v_span,
        v_limits[2] - 0.055v_span,
        length=v_count,
    )

    arrow_u = Float64[]
    arrow_v = Float64[]
    arrow_du = Float64[]
    arrow_dv = Float64[]

    for v in v_grid, u in u_grid
        du, dv = gm_rhs(u, v)

        normalized_du = du / u_span
        normalized_dv = dv / v_span
        normalized_speed = hypot(normalized_du, normalized_dv)

        if normalized_speed > 1.0e-10
            push!(arrow_u, u)
            push!(arrow_v, v)
            push!(
                arrow_du,
                arrow_fraction * u_span * normalized_du / normalized_speed,
            )
            push!(
                arrow_dv,
                arrow_fraction * v_span * normalized_dv / normalized_speed,
            )
        end
    end

    normalized_lengths = hypot.(arrow_du ./ u_span, arrow_dv ./ v_span)
    @assert all(
        isapprox(length_value, arrow_fraction; atol=1.0e-12)
        for length_value in normalized_lengths
    )

    return arrow_u, arrow_v, arrow_du, arrow_dv
end


# -----------------------------------------------------------------------------
# Shared plot styling
# -----------------------------------------------------------------------------

const u_limits = (-0.035, 1.35)
const v_limits = (-0.070, 2.70)

function base_plot(; show_legend=false)
    return plot(
        xlims=u_limits,
        ylims=v_limits,
        size=(1100, 780),
        dpi=180,
        background_color=hydra_navy,
        background_color_inside=hydra_navy,
        foreground_color_subplot=hydra_muted,
        foreground_color_axis=hydra_muted,
        foreground_color_border=hydra_muted,
        foreground_color_text=hydra_muted,
        guidefontcolor=hydra_white,
        guidefontsize=17,
        tickfontcolor=hydra_muted,
        tickfontsize=11,
        legend=show_legend ? :topleft : false,
        legendfontcolor=hydra_white,
        legendfontsize=10,
        legend_background_color=hydra_panel,
        legend_foreground_color=hydra_line,
        xlabel="u",
        ylabel="v",
        grid=false,
        framestyle=:axes,
        margin=6Plots.mm,
    )
end

function add_nullclines!(figure; labels=true, alpha_value=1.0)
    f_u_values = range(K * mu_u / a, u_limits[2], length=500)
    f_v_values = @. (a / mu_u) * f_u_values - K

    g_u_values = range(0.0, u_limits[2], length=500)
    g_v_values = @. (b / mu_v) * g_u_values^2

    plot!(
        figure,
        [0.0, 0.0],
        [0.0, v_limits[2]],
        color=hydra_gold,
        linewidth=2.6,
        alpha=alpha_value,
        label="",
    )
    plot!(
        figure,
        f_u_values,
        f_v_values,
        color=hydra_gold,
        linewidth=3.0,
        alpha=alpha_value,
        label=labels ? "f(u,v) = 0" : "",
    )
    plot!(
        figure,
        g_u_values,
        g_v_values,
        color=hydra_cyan,
        linewidth=3.0,
        alpha=alpha_value,
        label=labels ? "g(u,v) = 0" : "",
    )
end

function add_steady_states!(figure; labels=true, annotations=false)
    scatter!(
        figure,
        [0.0],
        [0.0],
        color=hydra_cyan,
        markerstrokecolor=hydra_white,
        markerstrokewidth=1.4,
        markersize=8.5,
        label=labels ? "(0, 0)  stable node" : "",
    )
    scatter!(
        figure,
        [u_minus],
        [v_minus],
        color=hydra_red,
        markerstrokecolor=hydra_white,
        markerstrokewidth=1.4,
        markersize=8.5,
        markershape=:diamond,
        label=labels ? "(0.5, 0.5)  saddle" : "",
    )
    scatter!(
        figure,
        [u_plus],
        [v_plus],
        color=hydra_green,
        markerstrokecolor=hydra_white,
        markerstrokewidth=1.4,
        markersize=9.0,
        markershape=:star5,
        label=labels ? "(1, 2)  stable spiral" : "",
    )

    if annotations
        annotate!(
            figure,
            0.07,
            0.10,
            text("(0, 0)", 11, hydra_cyan, :left),
        )
        annotate!(
            figure,
            u_minus - 0.03,
            v_minus + 0.13,
            text("(0.5, 0.5)", 11, hydra_red, :right),
        )
        annotate!(
            figure,
            u_plus + 0.035,
            v_plus + 0.12,
            text("(1, 2)", 11, hydra_green, :left),
        )
    end
end


# -----------------------------------------------------------------------------
# Figure 1: nullclines and classified steady states
# -----------------------------------------------------------------------------

nullclines_figure = base_plot(show_legend=true)
add_nullclines!(nullclines_figure; labels=true, alpha_value=1.0)
add_steady_states!(nullclines_figure; labels=true, annotations=false)


# -----------------------------------------------------------------------------
# Figure 2: phase portrait with equal-length arrows and trajectories
# -----------------------------------------------------------------------------

phase_portrait = base_plot(show_legend=false)

arrow_u, arrow_v, arrow_du, arrow_dv = normalized_vector_field(
    (0.0, u_limits[2]),
    (0.0, v_limits[2]);
    arrow_fraction=0.027,
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

add_nullclines!(phase_portrait; labels=false, alpha_value=0.62)

zero_basin_initial_states = (
    (0.22, 0.85),
    (0.38, 1.55),
)
positive_basin_initial_states = (
    (0.78, 0.25),
    (1.25, 2.55),
)

for initial_state in zero_basin_initial_states
    trajectory_u, trajectory_v = rk4_trajectory(initial_state, 30.0)
    plot!(
        phase_portrait,
        trajectory_u,
        trajectory_v,
        color=hydra_cyan,
        linewidth=2.8,
        alpha=0.92,
        label="",
    )
    scatter!(
        phase_portrait,
        [initial_state[1]],
        [initial_state[2]],
        color=hydra_cyan,
        markerstrokecolor=hydra_white,
        markerstrokewidth=1.0,
        markersize=5.8,
        label="",
    )

    @assert hypot(trajectory_u[end], trajectory_v[end]) < 0.08
end

for initial_state in positive_basin_initial_states
    trajectory_u, trajectory_v = rk4_trajectory(initial_state, 30.0)
    plot!(
        phase_portrait,
        trajectory_u,
        trajectory_v,
        color=hydra_gold,
        linewidth=2.8,
        alpha=0.92,
        label="",
    )
    scatter!(
        phase_portrait,
        [initial_state[1]],
        [initial_state[2]],
        color=hydra_gold,
        markerstrokecolor=hydra_white,
        markerstrokewidth=1.0,
        markersize=5.8,
        label="",
    )

    @assert hypot(
        trajectory_u[end] - u_plus,
        trajectory_v[end] - v_plus,
    ) < 0.08
end

add_steady_states!(phase_portrait; labels=false, annotations=true)


# -----------------------------------------------------------------------------
# Figure 3: maximum real eigenvalue for the relevant discrete modes
#
# Modes k = 0,...,3 contain both unstable modes and their nearest stable
# neighbours. Higher modes are increasingly strongly damped and would compress
# the positive part of the vertical scale.
# -----------------------------------------------------------------------------

displayed_modes = collect(0:4)
displayed_growth_rates = Float64[]

for k in displayed_modes
    _, matrix_k = mode_matrix(k)
    push!(
        displayed_growth_rates,
        maximum(real.(eigvals(matrix_k))),
    )
end

mode_growth_figure = plot(
    xlims=(-0.35, 4.35),
    ylims=(-1.20, 0.40),
    size=(1100, 700),
    dpi=180,
    background_color=hydra_navy,
    background_color_inside=hydra_navy,
    foreground_color_subplot=hydra_muted,
    foreground_color_axis=hydra_muted,
    foreground_color_border=hydra_muted,
    foreground_color_text=hydra_muted,
    guidefontcolor=hydra_white,
    guidefontsize=17,
    tickfontcolor=hydra_muted,
    tickfontsize=11,
    xlabel="spatial mode k",
    ylabel="unstable eigenvalue",
    xticks=(displayed_modes, string.(displayed_modes)),
    yticks=-1.1:0.1:0.4,
    grid=false,
    legend=false,
    framestyle=:axes,
    margin=6Plots.mm,
)

hspan!(
    mode_growth_figure,
    [0.0, 0.40],
    color=hydra_gold,
    alpha=0.055,
    label="",
)
hspan!(
    mode_growth_figure,
    [-1.2, 0.0],
    color=hydra_cyan,
    alpha=0.035,
    label="",
)
hline!(
    mode_growth_figure,
    [0.0],
    color=hydra_white,
    linewidth=1.1,
    alpha=0.75,
    label="",
)

for (k, growth_rate) in zip(displayed_modes, displayed_growth_rates)
    point_color = growth_rate > 0 ? hydra_gold : hydra_cyan

    plot!(
        mode_growth_figure,
        [k, k],
        [0.0, growth_rate],
        color=point_color,
        linewidth=3.0,
        alpha=0.82,
        label="",
    )
    scatter!(
        mode_growth_figure,
        [k],
        [growth_rate],
        color=point_color,
        markerstrokecolor=hydra_white,
        markerstrokewidth=1.4,
        markersize=8.5,
        label="",
    )

    label_offset = growth_rate > 0 ? 0.032 : -0.036
    vertical_alignment = growth_rate > 0 ? :bottom : :top

    annotate!(
        mode_growth_figure,
        k,
        growth_rate + label_offset,
        text(
            @sprintf("%+.3f", growth_rate),
            11,
            point_color,
            :center,
            vertical_alignment,
        ),
    )
end

annotate!(
    mode_growth_figure,
    2.55,
    0.335,
    text("unstable", 11, hydra_gold, :left),
)
annotate!(
    mode_growth_figure,
    2.0,
    -0.835,
    text("stable", 11, hydra_cyan, :left),
)


# -----------------------------------------------------------------------------
# Figure 4: unstable eigenvalues on the larger interval [0, 10]
#
# The Turing band is unchanged, but the Neumann spectrum is ten times denser in
# wave number and modes k = 2,...,22 become unstable.
# -----------------------------------------------------------------------------

large_domain_modes = collect(0:26)
large_domain_growth_rates = Float64[]

for k in large_domain_modes
    _, matrix_k = mode_matrix(k; length=large_domain_length)
    push!(
        large_domain_growth_rates,
        maximum(real.(eigvals(matrix_k))),
    )
end

large_domain_mode_growth_figure = plot(
    xlims=(-0.8, 26.8),
    ylims=(-0.36, 0.40),
    size=(1100, 700),
    dpi=180,
    background_color=hydra_navy,
    background_color_inside=hydra_navy,
    foreground_color_subplot=hydra_muted,
    foreground_color_axis=hydra_muted,
    foreground_color_border=hydra_muted,
    foreground_color_text=hydra_muted,
    guidefontcolor=hydra_white,
    guidefontsize=17,
    tickfontcolor=hydra_muted,
    tickfontsize=11,
    xlabel="spatial mode k",
    ylabel="unstable eigenvalue",
    xticks=(collect(0:2:26), string.(0:2:26)),
    yticks=-0.3:0.1:0.4,
    grid=false,
    legend=false,
    framestyle=:axes,
    margin=6Plots.mm,
)

hspan!(
    large_domain_mode_growth_figure,
    [0.0, 0.40],
    color=hydra_gold,
    alpha=0.055,
    label="",
)
hspan!(
    large_domain_mode_growth_figure,
    [-0.36, 0.0],
    color=hydra_cyan,
    alpha=0.035,
    label="",
)
hline!(
    large_domain_mode_growth_figure,
    [0.0],
    color=hydra_white,
    linewidth=1.1,
    alpha=0.75,
    label="",
)

for (k, growth_rate) in zip(
    large_domain_modes,
    large_domain_growth_rates,
)
    point_color = growth_rate > 0 ? hydra_gold : hydra_cyan

    plot!(
        large_domain_mode_growth_figure,
        [k, k],
        [0.0, growth_rate],
        color=point_color,
        linewidth=2.1,
        alpha=0.82,
        label="",
    )
    scatter!(
        large_domain_mode_growth_figure,
        [k],
        [growth_rate],
        color=point_color,
        markerstrokecolor=hydra_white,
        markerstrokewidth=1.0,
        markersize=5.6,
        label="",
    )
end

annotate!(
    large_domain_mode_growth_figure,
    8.0,
    large_domain_growth_rates[9] + 0.035,
    text("dominant mode k = 8", 11, hydra_gold, :center, :bottom),
)
annotate!(
    large_domain_mode_growth_figure,
    19.5,
    0.355,
    text("unstable modes k = 2,...,22", 11, hydra_gold, :center),
)


# -----------------------------------------------------------------------------
# Save the requested images
# -----------------------------------------------------------------------------

mkpath(output_directory)
savefig(nullclines_figure, nullclines_output)
savefig(phase_portrait, phase_portrait_output)
savefig(mode_growth_figure, mode_growth_output)
savefig(
    large_domain_mode_growth_figure,
    large_domain_mode_growth_output,
)


# -----------------------------------------------------------------------------
# Numerical values used by the planned example slides
# -----------------------------------------------------------------------------

println("Gierer--Meinhardt example")
@printf(
    "  parameters: Du=%.4f  Dv=%.4f  a=%.4f  b=%.4f  mu_u=%.4f  mu_v=%.4f  K=%.4f\n",
    D_u,
    D_v,
    a,
    b,
    mu_u,
    mu_v,
    K,
)

println("\nSteady states and local eigenvalues")
for steady_state in steady_states
    u_value, v_value = steady_state.point
    eigenvalues = eigvals(gm_jacobian(u_value, v_value))
    @printf(
        "  %-30s (% .6f, % .6f)  eigenvalues = %s\n",
        steady_state.name,
        u_value,
        v_value,
        string(eigenvalues),
    )
end

@printf(
    "\nTuring band: %.9f < mu < %.9f\n",
    mu_band_minus,
    mu_band_plus,
)

println("\nNeumann modes on [0, 1]")
println("  k       mu_k             eigenvalues of M_k                 max Re(lambda)")

unstable_modes = Int[]
for k in 0:6
    mu_k, matrix_k = mode_matrix(k)
    eigenvalues = eigvals(matrix_k)
    growth_rate = maximum(real.(eigenvalues))

    if growth_rate > 0
        push!(unstable_modes, k)
    end

    @printf(
        "  %d   %12.6f   %-38s   % .9f\n",
        k,
        mu_k,
        string(eigenvalues),
        growth_rate,
    )
end


@assert unstable_modes == [1, 2]

large_domain_unstable_modes = Int[]
for k in 0:23
    _, matrix_k = mode_matrix(k; length=large_domain_length)
    growth_rate = maximum(real.(eigvals(matrix_k)))

    if growth_rate > 0
        push!(large_domain_unstable_modes, k)
    end
end

@assert large_domain_unstable_modes == collect(2:22)

println("\nUnstable modes: ", unstable_modes)
println(
    "Unstable modes on [0, 10]: ",
    large_domain_unstable_modes,
)
println("Saved nullclines to:     ", nullclines_output)
println("Saved phase portrait to: ", phase_portrait_output)
println("Saved mode growth plot to: ", mode_growth_output)
println(
    "Saved large-domain growth plot to: ",
    large_domain_mode_growth_output,
)
