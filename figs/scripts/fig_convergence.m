
close all
fig = journal_figure([6.5 2.5], 2, 'tex');

% Mode numbers
m = (0 : 8)';

% Load and set formatting parameters
sg8 = load('../../data/convergence/sigma_g8.dat');
sg16 = load('../../data/convergence/sigma_g16.dat');
sg32 = load('../../data/convergence/sigma_g32.dat');
sg64 = load('../../data/convergence/sigma_g64.dat');

cg = [0.4 0.8 0.4];

sg = 1/2 - m .* lambda_PSH(1, 0) .* lambda_PSH(m, m) / 4;
sg8(end-1:end, 2) = nan;

sb8 = load('../../data/convergence/sigma_b8.dat');
sb16 = load('../../data/convergence/sigma_b16.dat');
sb32 = load('../../data/convergence/sigma_b32.dat');
sb64 = load('../../data/convergence/sigma_b64.dat');

cb = [0.4 0.4 0.8];

sb = (1 / 96) * (m .* lambda_PSH(1, 0) .* lambda_PSH(m, m) / 4 - 9 ./ (2 * (m - 1) .* (2 * m + 1)));
sb(1:2) = 0;

i = 6;
cmapg = 1 - linearrgbmap(1 - cg, 5);
cmapg = cmapg(2:end, :);
cmapb = 1 - linearrgbmap(1 - cb, 5);
cmapb = cmapb(2:end, :);

sp1 = subplot(1, 2, 1);

plot(m, sg, 'o-','Color', cg, 'DisplayName', '{\it \sigma_g} exact'), hold on
plot(m, sg64(:, 2), 's--', 'Color', cg/2, 'DisplayName', '{\it \sigma_g} numerical'), hold on


plot(m, 96 * sb, 'o-','Color', cb, 'DisplayName', '{\it Ra_* \sigma_b} exact'), hold on
plot(m, 96 * sb64(:, 2), 's--', 'Color', cb/2, 'DisplayName', '{\it Ra_* \sigma_b} numerical'), hold on
xlabel('mode, {\it m}');
ylabel('stability coefficient, {\it \sigma^m}')
xlim([0 8])
legend('Location', 'east', 'Color', 'none')

sp2 = subplot(1, 2, 2);
res = [8 16 32 64];
errg = [abs(sg(i) - sg8(i, 2)), abs(sg(i) - sg16(i, 2)), abs(sg(i) - sg32(i, 2)), abs(sg(i) - sg64(i, 2))];
errb = [abs(sb(i) - sb8(i, 2)), abs(sb(i) - sb16(i, 2)), abs(sb(i) - sb32(i, 2)), abs(sb(i) - sb64(i, 2))];
loglog(res, errg, '-o', 'Color', cg, 'DisplayName', '{\it \sigma_g}'), hold on
loglog(res, errb, '-o', 'Color', cb, 'DisplayName', '{\it \sigma_b}'), hold on
loglog(res, 4 * errb(1)*(res / res(1)).^(-2.5), '--', 'Color', [0.4 0.4 0.4]', 'HandleVisibility','off')
legend('location', 'west')
annotation('textbox', [0.8 0.5 0.3 0.3], 'String', '\propto{\it M^{-2.5}}', 'FontName', 'Times', 'FontSize', gca().FontSize, 'FitBoxToText','on', 'EdgeColor','none');
xlabel('degree of expansion, {\it M}');
ylabel('absolute error, |{\it \sigma_{M}}^{5} - {\it \sigma}^5|')
xlim([0.8*8 1.2*64])
xticks(res)
xticklabels({'8', '16', '32', '64'})
ylim([1e-15 1e0])

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

function val = lambda_PSH(l, m)

  id = abs(m) <= l;
  val = (gamma((l + m + 1) / 2) .* gamma((l - m + 1) / 2)) ./ (gamma((l + m + 2) / 2) .* gamma((l - m + 2) / 2));
  val = val .* id;

end
