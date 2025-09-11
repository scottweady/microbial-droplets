
close all
fig = journal_figure([6.5 2.5], 2, 'tex');

% Load data
sigma_g = load('../../data/sigma_g.dat');
sigma_b = load('../../data/sigma_b.dat');

% Range of wavenumbers
m = 0 : 40;
Ra_c = (0.5 - sigma_g) ./ sigma_b;

% Trim wavenumber range
m = m(4 : end);
Ra_c = Ra_c(4 : end, 2)';

% Plot
mm = [m fliplr(m)];
RRa_c = [Ra_c Ra_c(1) * ones(size(m))];
fill(mm, RRa_c, [0.4 0.4 0.8].^0.125, 'LineWidth', 2, 'EdgeColor', 'none'), hold on
plot(m, Ra_c, 'ko-', 'LineWidth', 2)

xlim([3 40])
ylim([0 Ra_c(1)])
xticks([3 10 : 10 : 40])

xlabel('mode, {\it m}')
ylabel('Rayleigh number, {\it Ra}')

% Format
ax = gca;
ax.Units = 'inches';

ax.Position(3) = 0.45 * fig.PaperSize(1);
ax.Position(1) = 0.5 * fig.PaperSize(1) - 0.5 * ax.Position(3);
ax.FontSize = 18;

stable = annotation('textbox', 'string', 'stable', 'FontName', 'times', 'EdgeColor', 'none', 'BackgroundColor', 'none', 'FontSize', 18, 'Units', 'inches', 'HorizontalAlignment', 'center');
unstable = annotation('textbox', 'string', 'unstable', 'FontName', 'times', 'EdgeColor', 'none', 'BackgroundColor', 'none', 'FontSize', 18, 'Units', 'inches', 'HorizontalAlignment', 'center');

stable.Position(1) = 0.5 * fig.PaperSize(1) - 0.5 * stable.Position(3);
stable.Position(2) = 0.9;

unstable.Position(1) = 0.5 * fig.PaperSize(1) - 0.5 * unstable.Position(3);
unstable.Position(2) = 2.9;

