
using Plots
using DelimitedFiles

include("ProjectedSphericalHarmonics.jl")

# Create output directory
isdir("results") || mkdir("results")

# Rayleigh number
Ra = 1

# Growth rate
γ = 0

# Absorption rate
β = 1

# List of wavenumbers
mspan = collect(0 : 40)

# Preallocate eigenvalues
σₘ = zeros(length(mspan))

# Number of modes
M = 64

println("Discretizing projected spherical harmonics...")

# Discretize
D = disk(M)

# Get points and quadrature weights
ζ, dζ = D.ζ, D.dζ

# Radius
r = abs.(ζ)

println("Computing base state...")

# O(1) buoyancy terms
σ₀ = 2β * ones(length(ζ))
𝒱σ₀ = (β/2) * (r.^2 .- 1)
ℬσ₀ = (β/8) * (r.^4/4 - r.^2 .- 5/4)

# O(1) pressure terms
p₀ = -𝒩inv(γ .+ (Ra / 16) * 𝒱σ₀, D)
𝒮p₀ = 𝒮(p₀, D)

# O(1) normal velocity
U₀ = -∂n(𝒮p₀ .+ (Ra / 16) * ℬσ₀, D)

# Save axisymmetric fields
sol = [real.(ζ)  imag.(ζ)  dζ  real.(p₀)  real.(σ₀)]
writedlm("results/sol.dat", sol, ' ')

println("Beginning main loop...")

# Loop over mode numbers
for (nm, m) in enumerate(mspan)

	# O(ϵ) buoyancy terms
	σ₁ = zeros(length(ζ))
	𝒱σ₁ = zeros(length(ζ))
	δ𝒱σ₀ = β * (abs2.(ζ) .- (m > 0) * (1 / m)) .* ζ.^m

	# O(ϵ) pressure terms
	p₁ = -𝒩inv(δ𝒩(p₀, m, D) + (Ra / 16) * (𝒱σ₁ + δ𝒱σ₀), D)
	
	# O(ϵ) normal velocity
	ψ = 𝒮(p₁, D) .+ δ𝒮(p₀, m, D) .+ (Ra / 16) * (ℬ(σ₁, D) .+ δℬ(σ₀, m, D))
	U₁ = -(m + 1) * U₀ .- ∂n(ψ, D)

	# Get stability coefficient
	σₘ[nm] = real.(U₁[1]) 

	# Print 
	println("(m, σₘ) = ", "($m, ", σₘ[nm], ")")

end

# Save growth rate
writedlm("results/sigma.dat", [mspan σₘ], ' ')
