
close all
fig = journal_figure([6.5 2.5], 2, 'tex');

% Domain size
R = 4;

% Load solution
sol = load('../../data/sol_b.dat');
zeta = sol(:, 1) + 1i * sol(:, 2);
dzeta = sol(:, 3);
sigma = sol(:, 5);

% Number of evaluation points
N = 64;
x = linspace(-R, R, N);
z = linspace(-R, -1e-2, N);

% Meshgrid
[zz, xx] = meshgrid(z, x);
xx = xx(:); zz = zz(:);

% Single layer potential
S = dzeta.' ./ sqrt( (xx - real(zeta)').^2 + (imag(zeta)').^2 + zz.^2) / (4 * pi);

% Concentration
c = 1 - S * sigma;

% Reshape
xx = reshape(xx, [N N]);
zz = reshape(zz, [N N]);
c = reshape(c, [N N]);

% Plot
contourf(xx, zz, c, 0 : 0.1 : 1, 'EdgeColor', 'none')

axis equal tight

xlim([-R R])
ylim([-R 0])
xticks(-R : 2 : R)
yticks(-R : 0)

xlabel('{\it x}')
ylabel('{\it z}')

colormap(flipud(cmocean('matter')))
clim([0 1])

% Format
hb = colorbar;
hb.Title.String = '{\it c}';
hb.FontSize = 18;
