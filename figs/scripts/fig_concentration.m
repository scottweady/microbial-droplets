
close all
fig = journal_figure([6.5 2.5], 2, 'tex');

% Degree of eigenfunction expansion
M = 64;

% Discretize
[theta, dtheta] = legpts(M + 1, [0 pi/2]);
[phi, dphi] = trigpts(2 * M + 1, [0 2 * pi]);

% Reshape
dtheta = dtheta(:);
dphi = dphi(:);

% Disk points (include 1 / sqrt(1 - r^2) weight)
zeta = sin(theta) .* exp(1i * phi'); zeta = zeta(:);
wdzeta = sin(theta) .* dtheta .* dphi'; wdzeta = wdzeta(:);

% Density
sigma = (4 / lambda_PSH(0, 0)) * ones(length(zeta), 1);

% Number of evaluation points
N = 64;
x = linspace(-4, 4, N);
z = linspace(-4, -1e-2, N);

% Meshgrid
[zz, xx] = meshgrid(z, x);
xx = xx(:); zz = zz(:);

% Single layer potential
S = wdzeta.' ./ sqrt( (xx - real(zeta)').^2 + (imag(zeta)').^2 + zz.^2);

% Concentration
c = 1 - (1 / (4 * pi)) * S * sigma;

% Reshape
xx = reshape(xx, [N N]);
zz = reshape(zz, [N N]);
c = reshape(c, [N N]);

% Plot
contourf(xx, zz, c, 0 : 0.1 : 1, 'EdgeColor', 'none')

axis equal tight

xlim([-4 4])
ylim([-4 0])
xticks(-4 : 2 : 4)
yticks(-4 : 0)

xlabel('{\it x}')
ylabel('{\it z}')

colormap(flipud(cmocean('matter')))
clim([0 1])

% Format
hb = colorbar;
hb.Title.String = '{\it c}';
hb.FontSize = 18;

function val = lambda_PSH(l, m)
  id = abs(m) <= l;
  val = (gamma((l + m + 1) / 2) .* gamma((l - m + 1) / 2)) ./ (gamma((l + m + 2) / 2) .* gamma((l - m + 2) / 2));
  val = val .* id;
end