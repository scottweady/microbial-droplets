
mkdir("output")
cd scripts

fig_axisymmetric
print("../output/fig_axisymmetric", "-dpdf")

fig_concentration
print("../output/fig_concentration", "-dpdf")

sol = 'g';
fig_sigma
print("../output/fig_sigma_g", "-dpdf")

sol = 'b';
fig_sigma
print("../output/fig_sigma_b", "-dpdf")

fig_stability
print("../output/fig_stability", "-dpdf")

fig_transition
print("../output/fig_transition", "-dpdf")

cd ..