
close all
fig = journal_figure([6.5 2.5], 2, 'tex');

% Choose solution: 'g' for growth, 'b' for buoyancy
% sol = 'g';

% Load and set formatting parameters
if strcmp(sol, 'g')
  sigma = load('../../data/sigma_g.dat');
  name = '{\it \sigma_g}'; 
  ys = [-4 1]; 
  yt = -4 : 1; 
  loc = 'southwest'; 
  c = [0.4 0.8 0.4];
elseif strcmp(sol, 'b')
  sigma = load('../../data/sigma_b.dat');
  name = '{\it \sigma_b}'; 
  ys = [-0.01 0.05]; 
  yt = -0.01 : 0.01 : 0.05; 
  loc = 'northwest'; 
  c = [0.4 0.4 0.8];
end

m = sigma(:, 1);
sigma = sigma(:, 2);

% Plot
plot(m, sigma, 'o-', 'Color', c, 'DisplayName', name), hold on
plot(m, 0 * m, '--', 'Color', 'k', 'DisplayName', 'none', 'LineWidth', 1, 'HandleVisibility','off')
plot(m, sigma(end) * (m / m(end)).^0.5, ':', 'Color', [0.2 0.2 0.2], 'DisplayName', '\propto {\it m^{1/2}}')
xlabel('mode, {\it m}');
ylabel(strcat('growth rate, ', name))

xlim([0 40])
ylim(ys)
xticks(0 : 10 : 40)
yticks(yt)

legend('location', loc)

% Format
ax = gca;
ax.Units = 'inches';

ax.Position(3) = 0.45 * fig.PaperSize(1);
ax.Position(1) = 0.5 * fig.PaperSize(1) - 0.5 * ax.Position(3);
ax.Position(4) = 0.7 * fig.PaperSize(2);
ax.Position(2) = 0.5 * fig.PaperSize(2) - 0.5 * ax.Position(4);
ax.FontSize = 18;
