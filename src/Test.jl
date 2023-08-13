using DifferentialEquations
using Plots
using Random, Distributions
using LinearAlgebra

L = 6

function min_symm(x, l) = minimum(x, l-x)

function ParticleNormalForce(rᵢ, rⱼ, R₁, R₂; k=50.)

    r⃗ᵢⱼ = abs.(rᵢ .- rⱼ)

    δ = max(R₁ + R₂ - rᵢⱼ, 0)


    return k*δ .* r̂ᵢⱼ

end

function WallNormalForce(r, R, L; k=300.)

    δ⃗₊ = r .+ R .- L/2
    δ⃗₋ = r .- R .+ L/2

    δ⃗ = -δ⃗₊ .* (δ⃗₊ .> 0) .- δ⃗₋ .* (δ⃗₋ .< 0)

    return k*δ⃗

end


function DEM!(ddu, du, u, p, t)

    L, R = p

    fill!(ddu, 0.)

    for i in axes(u,1)
        rᵢ = u[i, :]
        for j in axes(u,1)
            rⱼ = u[j, :]
            if i != j
                ddu[i, :] .+= ParticleNormalForce(rᵢ, rⱼ, R, R)
            end
        end
        ddu[i, :] .+= WallNormalForce(rᵢ, R, L)
        ddu[i, 2] -= 5
    end

end


tspan = (0.0, 100)

axis = -2:.8:2
grid = Iterators.product(axis, axis) |> collect
N = length(grid)
u0 = zeros(N, 2)

for i in 1:N
    u0[i, :] .= grid[i]
end


du0 = rand(Normal(), (N,2))
# du0 = zeros(N,2)


R = .3

p = (L, R)

prob = SecondOrderODEProblem(DEM!, du0, u0, tspan, p)



sol = solve(prob, VelocityVerlet(), dt=.01)


anim = @animate for step in sol.u
    plot(st=:scatter, ms=2, aspect_ratio=:equal, legend=false)
    positions = step[2,:,:]
    for i in 1:N
        plot!([positions[i,1]], [positions[i,2]], seriestype=:scatter, markersize=20)
    end
    xlims!(-L/2, L/2)
    ylims!(-L/2, L/2)
end

gif(anim, "anim_fps15.gif", fps = 40)


