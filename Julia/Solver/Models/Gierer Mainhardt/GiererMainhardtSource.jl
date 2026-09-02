# models/GiererMainhardt/GiererMainhardtSource.jl

# ============================================================
# Gierer-Meinhardt system with source-density field
# ============================================================
#
# Model:
#
#     u_t  = Du u_xx  + a s u^2/(v+1) - μu u
#     v_t  = Dv v_xx  + b s u^2       - μv v
#     s_t  = Ds s_xx  + (α u - s)/τ
#
# Homogeneous stationary state:
#
#     s* = α u*
#     v* = (b α / μv) (u*)^3
#
# and u* > 0 solves
#
#     (μu b / μv) u^3 - a u^2 + μu/α = 0.
#
# When two positive equilibria exist, we use the larger one.
#
# ============================================================

τ0 = 1.0


# ------------------------------------------------------------
# Larger positive homogeneous stationary state
# ------------------------------------------------------------

function positive_stationary_state(p)

    # This formula assumes pu = pv = 0
    if p.pu != 0.0 || p.pv != 0.0
        error("positive_stationary_state currently assumes pu = pv = 0.")
    end

    α = p.α

    # Cubic:
    #
    #     A u^3 - a u^2 + C = 0
    #
    # with
    #
    #     A = μu b / μv
    #     C = μu / α

    A = p.μu * p.b / p.μv
    C = p.μu / α

    f(u) = A * u^3 - p.a * u^2 + C

    # Positive local minimum of the cubic
    ucrit = 2.0 * p.a / (3.0 * A)

    if f(ucrit) > 0.0
        error(
            "No positive nonzero stationary state for the current parameters. " *
            "Minimum of stationary cubic is $(f(ucrit))."
        )
    end

    # --------------------------------------------------------
    # We want the LARGER positive root.
    # It lies to the right of ucrit.
    # --------------------------------------------------------

    left = ucrit
    right = max(2.0 * ucrit, 1.0)

    # Find a point where f becomes positive again.
    while f(right) <= 0.0
        right *= 2.0

        if right > 1e6
            error("Could not bracket the upper positive stationary root.")
        end
    end

    # Bisection
    for _ in 1:100
        middle = (left + right) / 2.0

        if f(middle) <= 0.0
            left = middle
        else
            right = middle
        end
    end

    u0 = (left + right) / 2.0

    # Stationary relations
    sd0 = α * u0
    v0 = (p.b * α / p.μv) * u0^3

    return u0, v0, sd0
end


RDModel(
    id = :gierer_meinhardt,

    display_name = "Gierer-Meinhardt system",

    variables = (:u, :v, :sd),

    parameters = (

        τ = τ0,

        Du = 1e-2,
        Dv = 1e0,
        Dsd = 1.5e1 / τ0,

        a = 1.5,
        b = 2.0,

        μu = 0.5,
        μv = 1.0,

        pu = 0.0,
        pv = 0.0,

        # Controls the stationary relation:
        #
        #     sd = α u
        #
        # α appears ONLY in the sd equation.
        α = 6.0,

        ρ0 = 1.0,
        ρ1 = 1.5,
        # ρ1 = 0.5,
    ),

    initial = function (U, x, p)

        # ----------------------------------------------------
        # Compute the upper positive homogeneous equilibrium
        # ----------------------------------------------------

        u0, v0, sd0 = positive_stationary_state(p)

        println("Positive stationary state:")
        println("    u  = ", u0)
        println("    v  = ", v0)
        println("    sd = ", sd0)

        U.u .= u0 .+ 0.01 .* randn(length(x))
        U.v .= v0 .+ 0.01 .* randn(length(x))

        # ----------------------------------------------------
        # Initial sd profile
        # ----------------------------------------------------

        sd_base = @. p.ρ0 - p.ρ1 / 2 + p.ρ1 * x

        H = div(length(x), 2)

        profile = :Constant
        # profile = :FootHead
        # profile = :HeadFoot
        # profile = :FootHeadReverse
        # profile = :HeadFootReverse
        # profile = :ZigZag

        if profile == :FootHead

            U.sd .= sd_base

        elseif profile == :HeadFoot

            U.sd .= [sd_base[(H + 1):end]; sd_base[1:H]]

        elseif profile == :FootHeadReverse

            U.sd .= [
                sd_base[1:H];
                reverse(sd_base[(H + 1):end])
            ]

        elseif profile == :HeadFootReverse

            U.sd .= [
                reverse(sd_base[1:H]);
                sd_base[(H + 1):end]
            ]

        elseif profile == :ZigZag

            x1 = 0.2
            x2 = 0.9

            slope = p.ρ1

            U.sd .= @. ifelse(
                x < x1,
                p.ρ0 + slope * x,
                ifelse(
                    x < x2,
                    p.ρ0 + slope * (x - x1),
                    p.ρ0 + slope * (x - x2),
                )
            )

        elseif profile == :Constant

            U.sd .= sd0 .+ 0.01 .* randn(length(x))

        end

        return nothing
    end,


    reaction = function (F, U, x, p, t)

        @. F.u =
            p.a * U.sd * U.u^2 / (U.v + 1.0) -
            p.μu * U.u +
            p.pu

        @. F.v =
            p.b * U.sd * U.u^2 -
            p.μv * U.v +
            p.pv

        @. F.sd =
            (p.α * U.u - U.sd) / p.τ

        return nothing
    end,


    diffusion = (
        u = :Du,
        v = :Dv,
        sd = :Dsd,
    ),


    latex_equations = (
        raw"\partial_t u = D_u \partial_{xx} u + a s \frac{u^2}{v+1} - \mu_u u",
        raw"\partial_t v = D_v \partial_{xx} v + b s u^2 - \mu_v v",
        raw"\partial_t s = D_s \partial_{xx} s + {\alpha u-s}",
    ),
)