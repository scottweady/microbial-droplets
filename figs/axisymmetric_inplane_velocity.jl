using LinearAlgebra
using CairoMakie, ColorSchemes
using HCubature
using StaticArrays



set_theme!(fonts = (; regular = "Times", bold = "Times Bold", italic = "Times.kz"))


function Gg_plane(X, z)
    
    nr = norm(X)^2 + z^2
    return  1/(4π)*X*( 1/nr^(3/2) - 3*z^2/nr^(5/2) )   
    end

function Gg_z(X, z)
    
    return  z/(4π)*(norm(X)^2 -2*z^2 )/(norm(X)^2 + z^2)^(5/2)

end

function Gb_plane(X, z; R=1)

    term1 = 1/(128π)*X*( 1 + 2*log(R) + 4 *z/sqrt(norm(X)^2+z^2)  )
    term2 = 1/(128π)*X*( -2*log(-z + sqrt(norm(X)^2+z^2) ) )
    term3 =  -z*(6/(128π))*X/(-z + sqrt(norm(X)^2+z^2))

    return term1 + term2 + term3

end

function Gb_z(X, z; R=1)

    term1 = 1/(32π)*z*( -log(R) + z/sqrt(norm(X)^2+z^2) + 
                        log(-z + sqrt(norm(X)^2+z^2) ) )
 

    return term1

end



pg_func = r-> 4/π*sqrt(1 - r^2)

pb_func = r-> -4/π*sqrt(1 - r^2)/96*(1 + 4/3*(1 - r^2) )
σ_func = r -> 2


rtol = 1e-11
atol = 1e-12

#flow buoyancy


up = X -> hcubature(  rθ -> pb_func(rθ[1])*rθ[1]*Gg_plane( SVector(X[1] - rθ[1]*cos(rθ[2]), X[2] - rθ[1]*sin(rθ[2])), X[3])[1], 
                    (0,0), (1 ,2*π),rtol =rtol, atol = atol)[1]
uσ = X -> hcubature(rθ -> σ_func(rθ[1])*rθ[1] *Gb_plane(SVector(X[1] - rθ[1]*cos(rθ[2]), X[2] - rθ[1]*sin(rθ[2])), X[3] )[1] , 
                    (0,0), (1 ,2*π),rtol =rtol, atol = atol)[1]                  


wp = X -> hcubature(  rθ -> pb_func(rθ[1])*rθ[1]*Gg_z( SVector(X[1] - rθ[1]*cos(rθ[2]), X[2] - rθ[1]*sin(rθ[2])), X[3]), 
                    (0,0), (1 ,2*π),rtol =rtol, atol = atol)[1]
wσ = X -> hcubature(rθ -> σ_func(rθ[1])*rθ[1] *Gb_z(SVector(X[1] - rθ[1]*cos(rθ[2]), X[2] - rθ[1]*sin(rθ[2])), X[3] ) , 
                    (0,0), (1 ,2*π),rtol =rtol, atol = atol)[1] 


flow_buoyancy = (x,z) -> SVector( up(SVector(x, 0, z)) .+ uσ(SVector(x, 0, z)), wp(SVector(x, 0, z)) .+ wσ(SVector(x, 0, z)) )                    



#flow growth

integrand_ug = (X, rθ) -> pg_func(rθ[1])*Gg_plane( SVector(X[1] - rθ[1]*cos(rθ[2]), X[2] - rθ[1]*sin(rθ[2])), X[3])[1]
ug = X -> hcubature(  rθ -> integrand_ug(X, rθ)*rθ[1],
                    (0,0), (1 ,2*π),rtol =rtol, atol = atol)[1]                


wg = X -> hcubature(  rθ -> pg_func(rθ[1])*rθ[1]*Gg_z( SVector(X[1] - rθ[1]*cos(rθ[2]), X[2] - rθ[1]*sin(rθ[2])), X[3]), 
                    (0,0), (1 ,2*π),rtol =rtol, atol = atol)[1]

flow_growth = (x,z) -> SVector( ug(SVector(x, 0, z)), wg(SVector(x, 0, z)) )




x = LinRange(1.000,2.0, 500)[2:end]
y = [0]#-LinRange(0.001,2.0, 200)

color_b = RGBf(0.2, 0.4, 0.8)
color_g = RGBf(0.4, 0.8, 0.4)

u_b = zeros(500)
u_g = zeros(500)

u_b[1] = 0.0
u_g[1] = 1/2


u_b[2:end] .= [flow_buoyancy(x, 0)[1] for x in x]
u_g[2:end] .= [flow_growth(x, 0)[1] for x in x]


x = LinRange(1.000,2.0, 500)

fig_b = Figure(size = (800, 400), fontsize = 22)
ax_b = Axis(fig_b[1, 1] , xlabel = L"r", ylabel = L"\text{radial velocity, }u",
            xgridvisible = false, ygridvisible = false,
            xreversed = false,
            backgroundcolor = :white, aspect = 1)

xlims!(ax_b, 1, 2)
ylims!(ax_b, -0.5, 0.5)


l1 = lines!(ax_b, x, u_b*96, color = color_b, linewidth = 3)
l2 = lines!(ax_b, x, u_g, color = color_g, linewidth = 3)


fig_b



