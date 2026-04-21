
close all
fig = journal_figure([6.5 2.25], 2, 'tex');

% Load data
data = load('../../data/kappa_H.dat');

H = data(:,1);
kappa = data(:,2);

% Rayleigh number
Ra0 = 96;
Ra = Ra0 ./ (1 - 6*kappa);
c = [0.4 0.4 0.8];

Hspan = 7 ./ [8 1];

% Plot kappa
sp1 = subplot(1, 2, 1);
plot(H, kappa, 'Color', c, 'DisplayName', '{\it \kappa}');
hold on;
plot(H, 1 - log(2 * H), '--', 'Color', [0.2 0.2 0.2], 'DisplayName', '1 - log 2|{\it H}|');
plot(H, 0*H, ':', 'Color', [0.2 0.2 0.2], 'HandleVisibility', 'off')
xlabel('|{\it H}|'); ylabel('{\it \kappa}');
legend('Location', 'northeast');
xlim(Hspan);
set(gca, 'FontName', 'Times', 'FontSize', 18)

% Plot Rayleigh number
sp2 = subplot(1, 2, 2);
plot(H, Ra, 'Color', c, 'DisplayName', '{\it Ra_{*,\kappa}}')
xlabel('|{\it H}|')
ylabel('{\it Ra_{*,\kappa}}')
xlim(Hspan)
ylim([0 100])
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
a.Position(1) = sp1.Position(1);
a.Position(2) = sum(sp1.Position([2 4])) - a.Position(4);

b = annotation('textbox', 'string', '\it (b)', 'fontname', 'times', 'edgecolor', 'none', 'FontSize', sp1.FontSize, 'Units', 'inches');
b.Position(1) = sp2.Position(1);
b.Position(2) = sum(sp2.Position([2 4])) - b.Position(4);
