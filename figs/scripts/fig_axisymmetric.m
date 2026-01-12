
close all
fig = journal_figure([6.5 2.5], 2, 'tex');

% Degree of eigenfunction expansion
M = 64;

% Load growth solution
sol_g = load('../../data/sol_g.dat');

% Get variables
x = reshape(sol_g(:, 1), [M + 1, 2 * M + 1]);
y = reshape(sol_g(:, 2), [M + 1, 2 * M + 1]);
dx = reshape(sol_g(:, 3), [M + 1, 2 * M + 1]);
p_g = reshape(sol_g(:, 4), [M + 1, 2 * M + 1]);

% Load buoyancy solution
sol_b = load('../../data/sol_b.dat');
p_b = reshape(sol_b(:, 4), [M + 1, 2 * M + 1]);
sigma_b = reshape(sol_b(:, 5), [M + 1, 2 * M + 1]);

% Number of evaluation points
N = 1024;

% Exterior points
r_ext = linspace(1.01, 2, N)';

% Preallocate
Sp_g = zeros(N, 1);
Sp_b = zeros(N, 1);
Bs_b = zeros(N, 1);

% Loop over exterior points
for n = 1 : N

  dr = sqrt( (r_ext(n) - x).^2 + y.^2);
  arg_g = (1 /(4 * pi)) * 1 ./ dr;
  Sp_g(n) = sum(arg_g(:) .* p_g(:) .* dx(:));

  arg_b = (1 /(4 * pi)) * 1 ./ dr;
  Sp_b(n) = sum(arg_b(:) .* p_b(:) .* dx(:));
  arg_b = (1 /(8 * pi)) * dr.^2 .* (log(dr) - 1);
  Bs_b(n) = sum(arg_b(:) .* sigma_b(:) .* dx(:));

end

% Finite difference to get velocity
dr = diff(r_ext);
r_ext = 0.5 * (r_ext(1 : end - 1) + r_ext(2 : end));
u_g = -diff(Sp_g) ./ dr;
u_b = -diff(Sp_b + (1 / 16) * Bs_b) ./ dr;

% Get theta = 0 slice
r = x(:, 1);
p_g = p_g(:, 1);
p_b = p_b(:, 1);

% Plot pressure
sp1 = subplot(1, 2, 1);
plot(r, p_g / max(abs(p_g)), 'Color', [0.4 0.8 0.4], 'DisplayName', 'growth, {\it p_g}'), hold on
plot(r, p_b / max(abs(p_b)), 'Color', [0.4 0.4 0.8], 'DisplayName', 'buoyancy, {\it p_b}')
plot(r, 0 * r, 'k--', 'LineWidth', 1.5, 'HandleVisibility','off')

ylim([-1.25 1.25])
yticks(-1 : 1)
yticklabels({'\it p_{b,0}', '0', '\it p_{g,0}'})
legend('location', 'southeast', 'Color', 'none')

xlabel('{\it r}')
ylabel('pressure, {\it p}')

% Plot velocity
sp2 = subplot(1, 2, 2);
plot(r_ext, u_g, 'Color', [0.4 0.8 0.4], 'DisplayName', 'growth, {\it u_g}'), hold on
plot(r_ext, 50 * u_b, 'Color', [0.4 0.4 0.8], 'DisplayName', 'buoyancy, {\it u_b}'), hold on
plot(r_ext, 0 * r_ext, 'k--', 'LineWidth', 1.5, 'HandleVisibility','off')

ylim([-0.5 0.5])
yticks(-0.5 : 0.25 : 0.5)
legend('location', 'northeast', 'Color', 'none')

xlabel('{\it r}')
ylabel('radial velocity, {\it u}')

% Format
sp1.Units = 'inches';
sp2.Units = 'inches';

sp1.Position(1) = 1;
sp1.Position(2) = 0.75;
sp1.Position(3) = 0.4 * fig.PaperSize(1);
sp1.Position(4) = 0.8 * fig.PaperSize(2);

sp2.Position(1) = sum(sp1.Position([1 3])) + 1.25;
sp2.Position(2 : 4) = sp1.Position(2 : 4);

a = annotation('textbox', 'string', '\it (a)', 'FontName', 'Times', 'edgecolor', 'none', 'FontSize', 18, 'Units', 'inches');
a.Position(1) = sp1.Position(1);
a.Position(2) = sum(sp1.Position([2 4])) - a.Position(4);

b = annotation('textbox', 'string', '\it (b)', 'FontName', 'Times', 'edgecolor', 'none', 'FontSize', 18, 'Units', 'inches');
b.Position(1) = sp2.Position(1);
b.Position(2) = sum(sp2.Position([2 4])) - b.Position(4);
