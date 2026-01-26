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


rtol = 1e-4
atol = 1e-6

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



x = LinRange(-2,2.0, 200)
y = -LinRange(0.001,2.0, 200)

colormap = :ice

Mag_b = [norm(flow_buoyancy(x, y)) for x in x, y in y] 
max_vel = maximum(Mag_b)
Mag_b = Mag_b/maximum(Mag_b)

begin
fig_b = Figure(size = (800, 400), fontsize = 22, fonts = (;regular="CMU Serif"))
ax_b = Axis(fig_b[1, 1] , xlabel = L"x", ylabel = L"z",
            xgridvisible = false, ygridvisible = false,
            xreversed = false,
            xticksvisible = false,
            backgroundcolor = :white)
ax_b_2 = Axis(fig_b[1, 1] , ylabel = L"z",
            xgridvisible = false, ygridvisible = false,
            xreversed = true,
            xticksvisible = false,
            xminorticksvisible = false,
            backgroundcolor = :white)
hidedecorations!(ax_b_2)
fs = heatmap!(ax_b, x, y , Mag_b, colormap = colormap )
fs2 = heatmap!(ax_b_2, x, y , Mag_b*NaN)
# fig_b
# fs.colormap = Reverse(colormap)

# (0.05, x[end]), (y[5], y[end])
# st 
# starts = [Point2f(0,0), Point2f(1,1), Point2f(2,0)]
st1 = streamplot!(ax_b, X -> Point2f((flow_buoyancy(X[1], X[2])/max_vel)...), (0, 2), (y[10], y[end]), colormap = [:white, :white],
    gridsize = (17, 24), arrow_size = 10,  maxsteps = 800, density = 0.9, stepsize = 0.01) 
st2 = streamplot!(ax_b_2, X -> Point2f((flow_buoyancy(X[1], X[2])/max_vel)...), (0, 2), (y[10], y[end]), colormap = [:white, :white],
    gridsize = (17, 24), arrow_size = 10,  maxsteps = 800, density = 0.9, stepsize = 0.01) 
        

lines!(ax_b_2, [-1,1], [0,0], color = RGBf(209/255, 231/255, 206/255), linewidth = 15)
end
# save("flow_buoyancy.png", fig_b)
# fig_b



#flow growth

integrand_ug = (X, rθ) -> pg_func(rθ[1])*Gg_plane( SVector(X[1] - rθ[1]*cos(rθ[2]), X[2] - rθ[1]*sin(rθ[2])), X[3])[1]
ug = X -> hcubature(  rθ -> integrand_ug(X, rθ)*rθ[1],
                    (0,0), (1 ,2*π),rtol =rtol, atol = atol)[1]                


wg = X -> hcubature(  rθ -> pg_func(rθ[1])*rθ[1]*Gg_z( SVector(X[1] - rθ[1]*cos(rθ[2]), X[2] - rθ[1]*sin(rθ[2])), X[3]), 
                    (0,0), (1 ,2*π),rtol =rtol, atol = atol)[1]

flow_growth = (x,z) -> SVector( ug(SVector(x, 0, z)), wg(SVector(x, 0, z)) )

Mag_g = [norm(flow_growth(x, y)) for x in x, y in y]
max_mag_g = maximum(Mag_g)
Mag_g = Mag_g/maximum(Mag_g)

# begin
fig_g = Figure(size = (800, 400), fontsize = 22, fonts = (;regular="CMU Serif"))
ax_g = Axis(fig_g[1, 1] , xlabel = L"x", ylabel = L"z",
            xgridvisible = false, ygridvisible = false,
            xreversed = false,
            xticksvisible = false,
            backgroundcolor = :white)
ax_g_2 = Axis(fig_g[1, 1] , ylabel = L"z",
            xgridvisible = false, ygridvisible = false,
            xreversed = true,
            xticksvisible = false,
            xminorticksvisible = false,
            backgroundcolor = :white)
hidedecorations!(ax_g_2)
fs = heatmap!(ax_g, x, y , Mag_g, colormap = colormap )
fs2 = heatmap!(ax_g_2, x, y , Mag_g*NaN)
# fs.colormap = Reverse(colormap)

# (0.05, x[end]), (y[5], y[end])
# st 
# starts = [Point2f(0,0), Point2f(1,1), Point2f(2,0)]
st1 = streamplot!(ax_g, X -> Point2f((flow_growth(X[1], X[2])/max_mag_g)...), (0, 2), (y[10], y[end]), colormap = [:white, :white],
    gridsize = (17, 24), arrow_size = 10,  maxsteps = 800, density = 0.9, stepsize = 0.01) 
st2 = streamplot!(ax_g_2, X -> Point2f((flow_growth(X[1], X[2])/max_mag_g)...), (0, 2), (y[10], y[end]), colormap = [:white, :white],
    gridsize = (17, 24), arrow_size = 10,  maxsteps = 800, density = 0.9, stepsize = 0.01) 
        

lines!(ax_g_2, [-1,1], [0,0], color = RGBf(209/255, 231/255, 206/255), linewidth = 15)
# end
# save("flow_buoyancy.png", fig_b)

save("flow_growth.pdf", fig_g)
save("flow_buoyancy.pdf", fig_b)
fig_g
fig_b