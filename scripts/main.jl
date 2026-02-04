
using DelimitedFiles
using ProjectedSphericalHarmonics

# Create output directory
isdir("results") || mkdir("results")

# Dimensionless parameters
Ra, β, γ = 1, 1, 1

# List of wavenumbers
m = collect(0 : 40)

# Number of modes
M = 64

# Discretize
D = disk(M)

# Get points and quadrature weights
ζ, dζ = D.ζ, D.dζ

# Radius
r = abs.(ζ)

# Axisymmetric solution
σ₀ = 2β * ones(length(ζ))

pg = @. γ * λlm(1, 0) * sqrt(1 - r^2)
pb = @. -(β * Ra / 96) * λlm(1, 0) * sqrt(1 - r^2) * (1 +  (4 / 3) * (1 - r^2))

# Save axisymmetric fields
sol_g = [real.(ζ)  imag.(ζ)  dζ  real.(pg)  real.(σ₀)]
writedlm("results/sol_g.dat", sol_g, ' ')
sol_b = [real.(ζ)  imag.(ζ)  dζ  real.(pb)  real.(σ₀)]
writedlm("results/sol_b.dat", sol_b, ' ')

σg = 1/2 .- m .* λlm(1, 0) .* λlm.(m, m) / 4;
σb = (β * Ra / 96) * (m .* λlm(1, 0) .* λlm.(m, m) / 4 .- 9 ./ (2 * (m .- 1) .* (2 * m .+ 1)));
σb[1:2] .= 0.0

# Save growth rate
writedlm("results/sigma_g.dat", [m σg], ' ')
writedlm("results/sigma_b.dat", [m σb], ' ')