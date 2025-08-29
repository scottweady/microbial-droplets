
using DelimitedFiles

include("ProjectedSphericalHarmonics.jl")

# Create output directory
isdir("results") || mkdir("results")

# Rayleigh number
Ra = 1

# Growth rate
γ = 0

# List of wavenumbers
mspan = collect(0 : 40)

# Preallocate eigenvalues
σₘ = zeros(length(mspan))

# Number of modes
M = 64

println("Discretizing projected spherical harmonics...")

# Discretize
D = disk(M)

# Get points and quadrature weight
ζ, dζ = D.ζ, D.dζ

# Radius
r = abs.(ζ)

println("Computing base state...")

# O(1) concentration terms
σ̂ = 4 / μlm(0, 0)
σ₀ = @. σ̂ / sqrt(1 - r^2)
𝒱σ₀ = @. σ̂ * (log(1 + sqrt(1 - r^2)) - sqrt(1 - r^2))
ℬσ₀ = @. σ̂ * (-sqrt(1 - r^2) * ((1/9) * r^2 + (11/36)) + (1 / 4) * r^2 * (log(r) - 1) + ((1 / 4) * r^2 + 1/6) * atanh(sqrt(1 - r^2)) + (1/6) * log(r))

# O(1) pressure terms
p₀ = -𝒩inv(γ .+ (Ra / 16) * 𝒱σ₀, D)
𝒮p₀ = 𝒮(p₀, D)

# O(1) normal velocity
U₀ = -∂n(𝒮p₀ .+ (Ra / 16) * ℬσ₀, D)

println("Beginning main loop...")

# Loop over mode numbers
for (nm, m) in enumerate(mspan)

	# O(ϵ) concentration terms
	σ₁ = -ζ.^m .* σ₀

	if m == 0
		𝒱σ₁ = 𝒱σ₀ .+ σ̂
	else
		𝒱σ₁ = σ̂ * (1 .- (2 * m + 1) / (2 * m) * (sqrt.(1 .- abs2.(ζ)) .+ Aₘ(ζ, m))) .* ζ.^m
	end

	# O(ϵ) pressure terms
	p₁ = -𝒩inv(δ𝒩(p₀, m, D) + (Ra / 16) * 𝒱σ₁, D)
	
	# O(ϵ) normal velocity
	ψ = 𝒮(p₁, D) .+ δ𝒮(p₀, m, D) .+ (Ra / 16) * (ℬ(σ₁, D) .+ δℬ(σ₀, m, D))
	U₁ = -(m + 1) * U₀ .- ∂n(ψ, D)

	# Get stability coefficient
	σₘ[nm] = real.(U₁[1]) 

	# Print 
	println("(m, σₘ) = ", "($m, ", σₘ[nm], ")")

end

# Save axisymmetric fields
sol = [real.(ζ)  imag.(ζ)  dζ  real.(p₀)  real.(σ₀)]
writedlm("results/sol.dat", sol, ' ')

# Save growth rate
writedlm("results/sigma.dat", [mspan σₘ], ' ')
