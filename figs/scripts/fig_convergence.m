
close all
fig = journal_figure([6.5 2.5], 2, 'tex');

% Mode numbers
m = (0 : 8)';

% Numerical values from simulations
sg_num = [0.5000000000000019, -2.55351295663786e-15, -0.25000000000000444, -0.4375000000000049, -0.5937500000000036, -0.7304687500000102, -0.8535156250000009, -0.966308593750024, -1.0710449218750018];
sb_num = [-6.841443289541467e-8, 7.774365682571746e-7, -0.001560232528819519, 0.006420550427962096, 0.009661147579125777, 0.011756958761655648, 0.013383751645312816, 0.014759831915538264, 0.01597725303247551];

% Exact values
sg = 1/2 - m .* lambda_PSH(1, 0) .* lambda_PSH(m, m) / 4;
sb = (1 / 96) * (m .* lambda_PSH(1, 0) .* lambda_PSH(m, m) / 4 - 9 ./ (2 * (m - 1) .* (2 * m + 1)));
sb(1:2) = 0;

cb = [0.4 0.4 0.8];
cg = [0.4 0.8 0.4];


cmapg = 1 - linearrgbmap(1 - cg, 5);
cmapg = cmapg(2:end, :);
cmapb = 1 - linearrgbmap(1 - cb, 5);
cmapb = cmapb(2:end, :);

sp1 = subplot(1, 2, 1);

plot(m, sg, 'o-','Color', cg, 'DisplayName', '{\it \sigma_g} exact'), hold on
plot(m, sg_num, 's--', 'Color', cg/2, 'DisplayName', '{\it \sigma_g} numerical'), hold on


plot(m, 96 * sb, 'o-','Color', cb, 'DisplayName', '{\it Ra_* \sigma_b} exact'), hold on
plot(m, 96 * sb_num, 's--', 'Color', cb/2, 'DisplayName', '{\it Ra_* \sigma_b} numerical'), hold on
xlabel('mode, {\it m}');
ylabel('stability coefficient, {\it \sigma^m}')
xlim([0 8])
legend('Location', 'east', 'Color', 'none')


sp2 = subplot(1, 2, 2);

% Errors
errg = abs([-0.7304687500000027, -0.7304687500000129, -0.7304687500000102, -0.7304687500000129, -0.7304687500000013] - sg(6));
errb = abs([0.01147270745963821, 0.01177776684299821, 0.011756958761655648, 0.0117519910113929, 0.011752038480886295] - sb(6));
res = [8 16 32 64 128];

loglog(res, errg, '-o', 'Color', cg, 'DisplayName', '{\it \sigma_g}'), hold on
loglog(res, errb, '-o', 'Color', cb, 'DisplayName', '{\it \sigma_b}'), hold on
loglog(res, 32 * errb(end)*(res / res(end)).^(-4), '--', 'Color', [0.4 0.4 0.4]', 'HandleVisibility','off')
legend('location', 'west')
annotation('textbox', [0.85 0.5 0.3 0.3], 'String', '\propto{\it M^{-4}}', 'FontName', 'Times', 'FontSize', gca().FontSize, 'FitBoxToText','on', 'EdgeColor','none');
xlabel('degree of expansion, {\it M}');
ylabel('absolute error, |{\it \sigma_{M}}^{5} - {\it \sigma}^5|')
xlim([0.8*8 1.2*128])
xticks(res)
xticklabels({'8', '16', '32', '64', '128'})
ylim([1e-16 1e-1])

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
