
using AssociatedLegendrePolynomials
using FastGaussQuadrature
using FFTW
using LinearAlgebra
using SpecialFunctions

""" Gauss-Legendre quadrature points and weights on the interval dom """
function legpts(N, dom)

  x, dx = gausslegendre(N)
  x = @. dom[1] * (1 - x) / 2 + dom[2] * (1 + x) / 2
  dx = @. dx * (dom[2] - dom[1]) / 2
  return x, dx
  
end

""" Equispaced quadrature points and weights on the interval dom """
function trigpts(N, dom)

  x = range(dom[1], dom[2], N + 1)
  dx = diff(x)
  x = x[1:N]
  return x, dx

end

""" Projected spherical discretization of the unit disk """
function diskpts(Nr, Nθ, rspan=[0, 1], θspan=[0, 2π])

  # Radial grid
  s, ds = legpts(Nr, sqrt.(1 .- rspan.^2))

  # Angular grid
  θ, dθ = trigpts(Nθ, θspan)
  
  # Complex coordinates
  ζ = sqrt.(1 .- s.^2) * exp.(im * θ')
  ζ = reshape(ζ, :, 1)

  # Volume element
  dζ = -(s .* ds) * dθ'
  dζ = reshape(dζ, :, 1)

  return ζ, dζ

end

""" Projected spherical harmonics """ 
function Ylm(M::Int, ζ)

  r, θ = abs.(ζ), angle.(ζ)

  # Compute associated Legendre polynomials
  p = λlm(0 : M, 0 : M, sqrt.(1 .- r.^2)) * sqrt(2)

  # Fill in negative values  
  P = Array{ComplexF64}(undef, length(ζ), M + 1, 2 * M + 1)
  P[:, :, (M + 1) : (2 * M + 1)] = p[:, 1, :, 1 : (M + 1)]
  P[:, :, M : -1 : 1] = p[:, 1, :, 2 : (M + 1)]

  # Preallocate angular part and phase factor
  Z = Array{ComplexF64}(undef, length(ζ), M + 1, 2 * M + 1)

  # Loop over mode numbers and fill in arrays
  for m = -M : M

    nm = (M + 1) + m

    for l = max(abs(m), 0) : M
    
      nl = l + 1
      Z[:, nl, nm] = ϕlm(l, m) * exp.(im * m * θ)

    end
  end

  return P .* Z

end

""" Normal derivative of projected spherical harmonics """
function ∂Ylm∂n(l::Int, m::Int, ζ)

  θ = angle.(ζ)

  if mod(m + l, 2) == 0

    lpm = l + abs(m)
    lmm = l - abs(m)
    tmp = 0.5 * (loggamma(lpm + 1) + loggamma(lmm + 1)) - (loggamma(lpm/2 + 1) + loggamma(lmm/2 + 1)) - l * log(2)

    return ϕlm(l, m) * exp.(im * m * θ) * (-1)^(Int(lpm/2)) * (l + lpm * lmm) * sqrt((2 * l + 1) / 2π) * exp(tmp) 

  end

  return Inf

end

""" Phase factor """
function ϕlm(l::Int, m::Int)
  return m >= 0 ? 1.0 : (-1)^m
end

""" Eigenvalues of projected spherical harmonics """
function μlm(l::Int, m::Int)
  return exp((loggamma((l + m + 1) / 2) + loggamma((l - m + 1) / 2)) - (loggamma((l + m + 2) / 2) + loggamma((l - m + 2) / 2)))
end

""" Single layer for 3D Laplacian """
function 𝒮(u, Ω)

  # Even expansion of u * w
  uwₖ = psh_transform(u .* Ω.w, Ω, kind=:even)

  # Compute ∂Ylm∂n
  fₖ = Ω.S .* uwₖ

  # Evaluate on grid
  return Ω.Y.even * fₖ

end

""" Inverse of 𝒮 """
function 𝒮inv(f, Ω)

  # Even expansion of f
  fₖ = psh_transform(f, Ω, kind=:even)

  # Compute weighted coefficients
  uwₖ = (1 ./ Ω.S) .* fₖ

  # Evaluate on grid
  return (Ω.Y.even * uwₖ) ./ Ω.w

end

""" Hypersingular operator (Δ𝒮) """
function 𝒩(u, Ω)

  # Odd expansion of u
  uₖ = psh_transform(u, Ω, kind=:odd)

  # Compute weighted coefficients
  fwₖ = Ω.N .* uₖ

  # Evaluate on grid
  return (Ω.Y.odd * fwₖ) ./ Ω.w

end


""" Inverse of 𝒩 """
function 𝒩inv(f, Ω)

  # Weighted odd expansion of f 
  fwₖ = psh_transform(f .* Ω.w, Ω, kind=:odd)

  # Compute coefficients
  uₖ = (1 ./ Ω.N) .* fwₖ

  # Evaluate on grid
  return Ω.Y.odd * uₖ

end

""" Bilaplace operator """
function ℬ(u, Ω)

	ζ, dζ = Ω.ζ, Ω.dζ
	ζt = transpose(ζ)
	B = @. (1/8π) * abs(ζ - ζt)^2 * (log(abs(ζ - ζt)) - 1) * dζ';
	B[diagind(B)] .= 0
	return B * u

end


""" Shape derivative of 𝒮 """
function δ𝒮(u, m, Ω)

	ζ = Ω.ζ
	fac = ζ.^0
	arg = ζ.^m .* u
	val = 2 * (m + 1) * 𝒮(arg, Ω)
	
	for _ = 0 : m
		val .+= -𝒮(arg, Ω) .* fac
		fac .*= ζ
		arg ./= ζ
	end

	return val

end

""" Shape derivative of 𝒩 """
function δ𝒩(u, m, Ω)

	ζ = Ω.ζ
	fac = ζ.^0
	arg = ζ.^m .* u
	val = 2 * (m + 1) * 𝒩(arg, Ω)
	
	for _ = 0 : m
		val .+= -3 * 𝒩(arg, Ω) .* fac
		fac .*= ζ
		arg ./= ζ
	end

	return val

end

""" Shape derivative of ℬ """
function δℬ(u, m, Ω)

	ζ = Ω.ζ
	fac = ζ.^0
	arg = ζ.^m .* u
	val = 2 * (m + 1) * ℬ(arg, Ω)

	for _ = 0 : m
		val += 2 * ℬ(arg, Ω) .* fac
		val += (1 / 8π) * fac .* sum(abs2.(ζ .- transpose(ζ)) .* transpose(arg .* dζ), dims=2)
		fac .*= ζ
		arg ./= ζ
	end

	return val

end


""" Integral associated with the concentration perturbation """
function Aₘ(ζ, m; tol=1e-14, maxN=256, N0=16, dN=16)

	# Preallocate
  val = similar(ζ)

	# Loop over grid points
  for (n, ζₙ) in enumerate(ζ)

		# Origin
      if ζₙ == 0
          val[n] = 0
          continue
      end

		# Get limit of integration
    θmax = asin(abs(ζₙ))

		# Initial number of grid points
    N = N0
    prevI = Inf
    I = 0.0

		# Integrate and refine until convergence
    while true

      θ, dθ = legpts(N, [0, θmax])
      I = sum(sin.(θ).^(2 * m + 1) .* dθ)

      if abs(I - prevI) < tol * max(1.0, abs(I))
        break
      elseif N >= maxN
        @warn "Max grid points $maxN reached for ζ = $ζₙ without meeting tolerance $tol"
        break
      end

			# Assign current value as previous and refine
      prevI = I
      N += dN

    end

		# Store value
    val[n] = exp(log(I) - m * log(abs2(ζₙ)))

  end

  return val

end

""" Normal derivative of projected spherical harmonic expansion """
function ∂n(u, Ω)

  # Compute even expansion
  uₖ = psh_transform(u, Ω, kind=:even)

  # Evaluate on boundary
  return Ω.∂Y∂n.even * uₖ

end

""" Discretization of the unit disk """
function disk(M::Int)

  Nr = M + 1
  Nθ = 2 * M + 1

  # Compute interior grid
  ζ, dζ = diskpts(Nr, Nθ)

  # Weight function
  w = sqrt.(1 .- abs2.(ζ))

  # Preallocations
  μ = Array{ComplexF64}(undef, M + 1, 2 * M + 1) #eigenvalues
  odd = falses(M + 1, 2 * M + 1) #odd boolean
  even = falses(M + 1, 2 * M + 1) #even boolean
  modes = Array{Tuple{Int, Int}}(undef, M + 1, 2M + 1) #mode pair

  # Loop over mode numbers and fill in arrays
  for m = -M : M
      
    nm = (M + 1) + m

    for l = max(abs(m), 0) : M

      nl = l + 1

      μ[nl, nm] = μlm(l, m)
      even[nl, nm] = mod(l + m, 2) == 0
      odd[nl, nm] = mod(l + m, 2) == 1
      modes[nl, nm] = (l, m)

    end
  end

  # Evaluate eigenfunctions
  Y = Ylm(M, ζ)
  Y = (even = Y[:, even], odd = Y[:, odd])

  # Create map from mode pair to index
  modeIndex = Dict{Tuple{Int, Int}, Int}()
  for (idx, mode) in enumerate(modes[even])
      modeIndex[mode] = idx
  end

  for (idx, mode) in enumerate(modes[odd])
    modeIndex[mode] = idx
  end

  # Eigenvalues of singular operators
  S = +μ[even] / 4
  N = -1 ./ μ[odd]

  # Construct boundary
  θ, _ = trigpts(Nθ, [0, 2π])
  X = exp.(im * θ)

  # Normal derivatives of eigenfunctions on boundary (only even ones are valid)
  ∂Y∂n = Array{ComplexF64}(undef, length(X), M + 1, 2 * M + 1)

  for m = -M : M

    nm = (M + 1) + m

    for l = max(abs(m), 0) : M

      nl = l + 1

      ∂Y∂n[:, nl, nm] .= ∂Ylm∂n(l, m, X)

    end
  end

  ∂Y∂n = (even = ∂Y∂n[:, even], odd = ∂Y∂n[:, odd])

  # Store
  return (Y = Y, ∂Y∂n = ∂Y∂n, ζ = ζ, dζ = dζ, w = w, S = S, N = N, odd = odd, even = even, modes = modes, modeIndex = modeIndex, M = M)

end


""" Projected spherical harmonics transform using direct integration """
function psh_transform_direct(u, Ω; kind=:even)

  # Get basis functions
  Y = getfield(Ω.Y, kind)

  # Compute coefficients via orthogonality
  uₖ = Y' * (u .* Ω.dζ ./ Ω.w)

  return uₖ

end

""" Projected spherical harmonics transform using FFT (~10x faster for M = 64)"""
function psh_transform_fft(u, Ω; kind=:even)

  # Degree of radial expansion
  M = Ω.M

  # Quadrature weights
  ζ = @view Ω.ζ[1 : (M + 1)]
  dζ = @view Ω.dζ[1 : (M + 1)] 

  # Basis functions
  Y = @view getfield(Ω.Y, kind)[1 : (M + 1), :]

  # Compute transform
  u = reshape(u, M + 1, 2 * M + 1)
  uₖ = fft(u, 2)
  uₖ = Y' * (uₖ .* (dζ ./ sqrt.(1 .- abs2.(ζ))))

  # Get relevant coefficients
  modes = Ω.modes[getfield(Ω, kind)]
  azimuthal_modes = [mod(m, 2 * M + 1) + 1 for (_, m) in modes]
  uₖ = [uₖ[i, j] for (i, j) in enumerate(azimuthal_modes)]

  return uₖ
  
end

# Wrapper for PSH transform
function psh_transform(u, Ω; kind=:even)
  return psh_transform_fft(u, Ω, kind=kind)
end