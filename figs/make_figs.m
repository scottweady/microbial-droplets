
mkdir("output")
cd scripts

fig_axisymmetric
print("../output/fig-axisymmetric", "-dpdf")

fig_concentration
print("../output/fig-concentration", "-dpdf")

sol = 'g';
fig_sigma
print("../output/fig-sigma_g", "-dpdf")
sol = 'b';
fig_sigma
print("../output/fig-sigma_b", "-dpdf")

fig_stability
print("../output/fig-stability", "-dpdf")

fig_transition
print("../output/fig-transition", "-dpdf")

fig_convergence
print("../output/fig-convergence", "-dpdf")

close all

cd ..
