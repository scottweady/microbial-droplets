
using DelimitedFiles
using ProjectedSphericalHarmonics
using Plots

# Create output directory
isdir("results") || mkdir("results")

# Dimensionless parameters
# Ra, β, γ = 1, 1, 0
Ra, β, γ = 0, 1, 1

# List of wavenumbers
mspan = collect(0 : 8)

# Preallocate eigenvalues
σₘ = zeros(length(mspan))

# Number of modes
M = 64

println("Discretizing projected spherical harmonics...")

# Discretize
D = psh_disk(M)

# Get points and quadrature weights
ζ, dζ = D.ζ, D.dζ

println("Computing base state...")

# O(1) terms
σ₀ = 2β * ones(length(ζ)) #concentration density
f₀ = (γ / 2β) * σ₀ .+ (Ra / 16) * 𝒱(σ₀, D)
p₀ = -𝒩⁻¹(f₀, D) #pressure
ψ₀ = 𝒮(p₀, D) .+ (Ra / 16) * ℬ(σ₀, D)
U₀ = -∂n(ψ₀, D) #normal velocity

# Save axisymmetric fields
sol = [real.(ζ)  imag.(ζ)  dζ  real.(p₀)  real.(σ₀)]
writedlm("results/sol.dat", sol, ' ')

println("Beginning main loop...")

# Loop over mode numbers
for (nm, m) in enumerate(mspan)

	# O(ϵ) terms
	σ₁ = zeros(length(ζ)) #concentration density
    f₁ = δ𝒩(p₀, m, D) .+ (γ / 2β) * σ₁ .+ (Ra / 16) * (𝒱(σ₁, D) .+ δ𝒱(σ₀, m, D))
	p₁ = -𝒩⁻¹(f₁, D) #pressure
	ψ₁ = 𝒮(p₁, D) .+ δ𝒮(p₀, m, D) .+ (Ra / 16) * (ℬ(σ₁, D) .+ δℬ(σ₀, m, D))
	U₁ = -(m + 1) * U₀ .- ∂n(ψ₁, D) #normal velocity

	# Store stability coefficient
	σₘ[nm] = real.(U₁[1]) 

	# Print 
	println("(m, σₘ) = ", "($m, ", σₘ[nm], ")")

end

# Save stability coefficients
writedlm("results/sigma_g$(M).dat", [mspan σₘ], ' ')
