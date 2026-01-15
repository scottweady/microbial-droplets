

close all
fig = journal_figure([6.5 2.5], 2, 'tex');

% Degree of eigenfunction expansion
M = 64;

% Load growth solution
sol_g = load('../../data/sol_g.dat');
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
p_g = p_g(:, 1);
p_b = p_b(:, 1);
r = x(:, 1);

% Range of Rayleigh numbers
Ra = 95 + 5 * (-2 : 2);
Ra(3) = 96;
  
% Color map
cmap = [0.4 0.8 0.4; ...
        0.7 0.95 0.7; ...
        0.1 0.1 0.1; ...
        0.7 0.7 0.95; ...
        0.4 0.4 0.8];

% Plot pressure for different Rayleigh numbers
sp1 = subplot(1, 2, 1);

  for n = 1 : length(Ra)

    if n == (length(Ra)+1)/2
      ls = '-';
    else
      ls = '--';
    end
  
    plot(r, p_g + Ra(n) * p_b, 'Color', cmap(n, :), 'LineStyle', ls, 'DisplayName', strcat("{\it Ra} = ", num2str(Ra(n)))), hold on
  
  end

  xlim([0.8 1]);

  xlabel('{\it r}')
  ylabel('pressure, {\it p}')

  legend('location', 'southeast', 'FontSize', 18)
  set(gca, 'FontName', 'Times', 'FontSize', 18)

% Plot velocity for different Rayleigh numbers
sp2 = subplot(1, 2, 2);
  
  for n = 1 : length(Ra)
    if n == (length(Ra)+1)/2
      ls = '-';
    else
      ls = '--';
    end
  
    plot(r_ext, u_g + Ra(n) * u_b, 'Color', cmap(n, :), 'LineStyle', ls, 'DisplayName', strcat("{\it Ra} = ", num2str(Ra(n)))), hold on
  
  end
  
  xlim([1 1.2])

  xlabel('{\it r}')
  ylabel('radial velocity, {\it u}')

  legend('location', 'southwest', 'FontSize', 18)
  set(gca, 'FontName', 'Times', 'FontSize', 18)

% Format
sp1.Units = 'inches';
sp2.Units = 'inches';

sp1.Position(1) = 1;
sp1.Position(2) = 0.75;

sp1.Position(3) = 0.4 * fig.PaperSize(1);
sp1.Position(4) = 0.8 * fig.PaperSize(2);

sp2.Position(1) = sum(sp1.Position([1 3])) + 1.25;
sp2.Position(2 : 4) = sp1.Position(2 : 4);

a = annotation('textbox', 'string', '\it (a)', 'fontname', 'times', 'edgecolor', 'none', 'FontSize', sp1.FontSize, 'Units', 'inches');
a.Position(1) = sp1.Position(1) - 0.05;
a.Position(2) = sum(sp1.Position([2 4])) - a.Position(4);

b = annotation('textbox', 'string', '\it (b)', 'fontname', 'times', 'edgecolor', 'none', 'FontSize', sp1.FontSize, 'Units', 'inches');
b.Position(1) = sp2.Position(1) - 0.05;
b.Position(2) = sum(sp2.Position([2 4])) - b.Position(4);
