using Quaternions
using StaticArrays
using LinearAlgebra
using DifferentialEquations
using Plots, PlotThemes

"""
**Body:** Represents a rigid body.
- `J`: Inertia matrix, kg⋅m² (constant or time-varying)
"""
struct Body
    J::Function
end

"""
**Wheels:** Reaction wheel assembly for spacecraft attitude control.
- `A`: Allocation matrix
- `M`: Pseudo-inverse of the allocation matrix
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

    Wheels(A::AbstractMatrix, J, c, τ, Ω) = new(A, pinv(A), J, c, τ, Ω)
end

""" Convert a quaternion into a 4×1 static vector. """
qvec(q::Quaternion{T}) where {T} = SVector{4, T}(q.s, imag_part(q)...)

# Control law parameters
const γ = 0.006
const β = 0.023

"""
Compute control torque.
- `φ`: Target orientation quaternion
- `λ`: Current orientation quaternion
- `ω`: Angular velocity vector
"""
function torque(φ, λ, ω)
    ε = φ \ λ # Orientation error quaternion
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

    du[1:4]  = qvec(dλ)
    du[5:7]  = dω
    du[8:11] = dΩ
end

""" Spacecraft inertia matrix, kg⋅m² """
J = _ -> Diagonal([0.050, 0.050, 0.035])

body = Body(J)

""" Wheels skew angle """
α = deg2rad(54.73)

""" Wheels allocation matrix """
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

""" Initial state vector """
u0 = [qvec(λ0)..., ω0..., Ω0...]

""" System parameters """
params = (φ = φ, wheels = wheels, body = body)

sys = ODEProblem(system!, u0, (0, 80), params)
sol = solve(sys, Tsit5())

theme(:rose_pine_dawn)
default(titlefont = 8, legendfontsize = 8, tickfontsize = 8, labelfontsize = 8, legend = :topright)

φ_aligned = sign(qvec(φ) · sol[1:4, end]) * qvec(φ) # Align with the solution
plot_λ = hline(φ_aligned, line = :dash, color = :gray, alpha = 0.5, label = "")
plot!(sol, idxs = (0, 1:4), tspan = (0, 40), label = ["λ₀" "λ₁" "λ₂" "λ₃"], xlabel = "")

plot_ω = hline([0], color = :gray, label = "")
plot!(sol, idxs = (0, 5:7), tspan = (0, 40), label = ["ω₁" "ω₂" "ω₃"], xlabel = "")

plot_Ω = hline([0], color = :gray, label = "")
plot!(sol, idxs = (0, 8:11), label = ["Ω₁" "Ω₂" "Ω₃" "Ω₄"], xlabel = "")

plot(plot_λ, plot_ω, plot_Ω, layout = (@layout [a b; c]), size = (800, 600))
savefig("results/attitude_control_simulation.pdf")
