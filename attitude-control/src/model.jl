using Quaternions
using LinearAlgebra
using DifferentialEquations

using Plots, PlotThemes
theme(:rose_pine_dawn)

"""
Body

Represents a rigid body.

Fields:
- `J`: Inertia tensor matrix, kg⋅m²
"""
struct Body
    J::Function
end

"""
Wheels

Represents a reaction wheel assembly for spacecraft attitude control.

Fields:
- `A`: Direction Cosine Matrix (DCM)
- `M`: Pseudo-inverse of the DCM
- `J`: Moment of inertia around wheel gimbal axis, kg⋅m²
- `c`: Damping coefficient, N⋅m⋅s
- `τ`: Maximum wheel torque, N⋅m
- `Ω`: Maximum wheel angular velocity, rad/s
"""
struct Wheels
    A::Matrix{Float64}
    M::Matrix{Float64}
    J::Float64
    c::Float64
    τ::Float64
    Ω::Float64

    function Wheels(A, J, c, τ, Ω)
        new(A, pinv(A), J, c, τ, Ω)
    end
end

"""
Compute control torque

φ: Target orientation quaternion
λ: Current orientation quaternion
ω: Angular velocity vector
"""
function torque(φ, λ, ω)
    γ = 0.006
    β = 0.023

    # Compute orientation error quaternion
    ε = φ \ λ

    -γ * sign(ε.s) * [imag_part(ε)...] - β * ω
end

""" System dynamics """
function system!(du, u, p, t)
    λ = Quaternion(u[1:4]...)
    ω = u[5:7]
    Ω = clamp.(u[8:11], -p.wheels.Ω, p.wheels.Ω)

    τ = torque(p.φ, λ, ω)
    m = clamp.(p.wheels.M * τ, -p.wheels.τ, p.wheels.τ)

    J = p.body.J(t)

    dλ = 0.5 * λ * Quaternion(0, ω...)
    dω = J \ (-ω × (J * ω) + p.wheels.A * m)
    dΩ = (m - p.wheels.c * Ω) / p.wheels.J

    du[1:4] = [dλ.s, imag_part(dλ)...]
    du[5:7] = dω
    du[8:11] = dΩ
end

""" Spacecraft inertia tensor, kg⋅m² """
J = t -> Diagonal([0.050, 0.050, 0.035])

body = Body(J)

""" Wheels skew angle """
α = deg2rad(54.73)

""" Wheels direction cosine matrix """
A = [sin(α) 0 -sin(α) 0; 0 sin(α) 0 -sin(α); cos(α) cos(α) cos(α) cos(α)]

wheels = Wheels(A, 1.15e-4, 1e-5, 4e-3, 600)

""" Target orientation quaternion """
φ = Quaternion(-1 / sqrt(2), -1 / sqrt(6), 1 / sqrt(6), -1 / sqrt(6))

""" Initial orientation quaternion """
λ0 = Quaternion(1.0)

""" Initial spacecraft angular velocity vector, rad/s """
ω0 = [0, 0, 0]

""" Initial wheels angular velocity vector, rad/s """
Ω0 = [0, 0, 0, 0]

u0 = [λ0.s, imag_part(λ0)..., ω0..., Ω0...]

""" System parameters """
params = (φ = φ, wheels = wheels, body = body)

prob = ODEProblem(system!, u0, (0, 80), params)
sol = solve(prob, Tsit5())

# Plot
default(titlefont = 8, legendfontsize = 8, tickfontsize = 8, labelfontsize = 8, legend = :topright)

φ_corrected = sign([φ.s, imag_part(φ)...] · sol[1:4, end]) * [φ.s, imag_part(φ)...]
plot_λ = hline(φ_corrected, line = :dash, color = :gray, label = "")
plot!(sol, idxs = (0, 1:4), tspan = (0, 40), label = ["λ₀" "λ₁" "λ₂" "λ₃"], xlabel = "")

plot_ω = hline([0], color = :gray, label = "")
plot!(sol, idxs = (0, 5:7), tspan = (0, 40), label = ["ω₁" "ω₂" "ω₃"], xlabel = "")

plot_Ω = hline([0], color = :gray, label = "")
plot!(sol, idxs = (0, 8:11), label = ["Ω₁" "Ω₂" "Ω₃" "Ω₄"], xlabel = "")

plot(plot_λ, plot_ω, plot_Ω, layout = (@layout [a b; c]), size = (800, 600))
savefig("results/attitude_control_simulation.pdf")
